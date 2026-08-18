/****************************************************************************
* Triangle-mesh Fischer-Burmeister contact implementation.
****************************************************************************/
#pragma once

#include <sofa/ncp/contact/MeshNCPContactForceField.h>
#include <sofa/ncp/contact/FischerBurmeisterContactForceField.inl>

#include <sofa/core/objectmodel/ComponentState.h>
#include <sofa/helper/logging/Messaging.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <utility>

namespace sofa::ncp
{

template<class T1, class T2>
MeshNCPContactForceField<T1, T2>::MeshNCPContactForceField()
    : Inherit()
    , d_meshFilename(initData(&d_meshFilename, "meshFilename",
          "Triangle OBJ used directly as the vessel contact surface."))
    , d_flipNormals(initData(&d_flipNormals, false, "flipNormals",
          "Flip the OBJ orientation. The final normal must point toward the admissible lumen."))
    , d_normalMode(initData(&d_normalMode, static_cast<unsigned int>(FaceNormal), "normalMode",
          "0=piecewise-exact face normal, 1=barycentrically interpolated vertex normal."))
    , d_contactOffset(initData(&d_contactOffset, Real(0), "contactOffset",
          "Distance subtracted from the lumen-side gap, typically the catheter radius."))
    , d_maxSearchDistance(initData(&d_maxSearchDistance, Real(0), "maxSearchDistance",
          "Maximum closest-point distance. Zero disables the limit and queries the whole mesh."))
    , d_bvhLeafSize(initData(&d_bvhLeafSize, 8u, "bvhLeafSize",
          "Maximum number of triangles stored in an AABB-tree leaf."))
    , d_ignoreBoundaryEdges(initData(&d_ignoreBoundaryEdges, true, "ignoreBoundaryEdges",
          "Pin a contact row when its closest feature lies on an open mesh boundary edge."))
    , d_boundaryBarycentricTolerance(initData(&d_boundaryBarycentricTolerance, Real(1e-8),
          "boundaryBarycentricTolerance", "Tolerance used to detect a closest point on a boundary edge."))
    , d_pinnedIndices(initData(&d_pinnedIndices, sofa::type::vector<unsigned int>{ 0u },
          "pinnedIndices", "Object1 point indices excluded from mesh contact; index 0 is pinned by default."))
    , d_meshVertexCount(initData(&d_meshVertexCount, sofa::Size(0), "meshVertexCount",
          "Number of loaded OBJ vertices."))
    , d_meshTriangleCount(initData(&d_meshTriangleCount, sofa::Size(0), "meshTriangleCount",
          "Number of nondegenerate loaded triangles."))
    , d_meshBoundaryEdgeCount(initData(&d_meshBoundaryEdgeCount, sofa::Size(0), "meshBoundaryEdgeCount",
          "Number of open boundary edges in the contact mesh."))
{
    d_meshVertexCount.setReadOnly(true);
    d_meshTriangleCount.setReadOnly(true);
    d_meshBoundaryEdgeCount.setReadOnly(true);
}

template<class T1, class T2>
void MeshNCPContactForceField<T1, T2>::init()
{
    Inherit::init();

    if (!loadMesh())
    {
        this->d_componentState.setValue(core::objectmodel::ComponentState::Invalid);
        return;
    }
}

template<class T1, class T2>
void MeshNCPContactForceField<T1, T2>::reinit()
{
    Inherit::reinit();

    if (!loadMesh())
    {
        this->d_componentState.setValue(core::objectmodel::ComponentState::Invalid);
        return;
    }
}

template<class T1, class T2>
bool MeshNCPContactForceField<T1, T2>::loadMesh()
{
    m_vertices.clear();
    m_objNormals.clear();
    m_generatedVertexNormals.clear();
    m_triangles.clear();
    m_triangleOrder.clear();
    m_bvhNodes.clear();

    d_meshVertexCount.setValue(0);
    d_meshTriangleCount.setValue(0);
    d_meshBoundaryEdgeCount.setValue(0);

    if (!d_meshFilename.isSet() || d_meshFilename.getValue().empty())
    {
        msg_error() << "meshFilename is required.";
        return false;
    }

    if (d_normalMode.getValue() > static_cast<unsigned int>(InterpolatedVertexNormal))
    {
        msg_error() << "normalMode must be 0=face or 1=interpolated vertex normal.";
        return false;
    }

    if (!(d_contactOffset.getValue() >= Real(0))
        || !std::isfinite(static_cast<double>(d_contactOffset.getValue())))
    {
        msg_error() << "contactOffset must be finite and nonnegative.";
        return false;
    }

    if (!(d_maxSearchDistance.getValue() >= Real(0))
        || !std::isfinite(static_cast<double>(d_maxSearchDistance.getValue())))
    {
        msg_error() << "maxSearchDistance must be finite and nonnegative.";
        return false;
    }

    if (!(d_boundaryBarycentricTolerance.getValue() >= Real(0))
        || !std::isfinite(static_cast<double>(d_boundaryBarycentricTolerance.getValue())))
    {
        msg_error() << "boundaryBarycentricTolerance must be finite and nonnegative.";
        return false;
    }

    if (d_bvhLeafSize.getValue() == 0)
    {
        msg_error() << "bvhLeafSize must be at least one.";
        return false;
    }

    const std::string filename = d_meshFilename.getFullPath();
    if (!loadObj(filename) || !finalizeMesh() || !buildBvh())
        return false;

    d_meshVertexCount.setValue(m_vertices.size());
    d_meshTriangleCount.setValue(m_triangles.size());

    msg_info() << "Loaded mesh contact surface '" << filename << "': vertices="
               << m_vertices.size() << ", triangles=" << m_triangles.size()
               << ", boundaryEdges=" << d_meshBoundaryEdgeCount.getValue()
               << ", bvhNodes=" << m_bvhNodes.size() << ".";
    return true;
}

template<class T1, class T2>
bool MeshNCPContactForceField<T1, T2>::loadObj(const std::string& filename)
{
    std::ifstream stream(filename);
    if (!stream)
    {
        msg_error() << "Cannot open OBJ contact mesh: " << filename;
        return false;
    }

    std::string line;
    sofa::Size lineNumber = 0;

    while (std::getline(stream, line))
    {
        ++lineNumber;
        std::istringstream input(line);
        std::string keyword;
        input >> keyword;

        if (keyword.empty() || keyword[0] == '#')
            continue;

        if (keyword == "v")
        {
            Real x = Real(0), y = Real(0), z = Real(0);
            if (!(input >> x >> y >> z))
            {
                msg_error() << "Invalid OBJ vertex at line " << lineNumber << ".";
                return false;
            }
            m_vertices.emplace_back(x, y, z);
            continue;
        }

        if (keyword == "vn")
        {
            Real x = Real(0), y = Real(0), z = Real(0);
            if (!(input >> x >> y >> z))
            {
                msg_error() << "Invalid OBJ normal at line " << lineNumber << ".";
                return false;
            }

            Vec3 normal(x, y, z);
            if (!normalize(normal))
                normal.clear();
            m_objNormals.push_back(normal);
            continue;
        }

        if (keyword != "f")
            continue;

        sofa::type::vector<FaceCorner> polygon;
        std::string token;
        while (input >> token)
        {
            FaceCorner corner;
            if (!parseFaceCorner(token, m_vertices.size(), m_objNormals.size(), corner))
            {
                msg_error() << "Invalid OBJ face corner '" << token << "' at line " << lineNumber << ".";
                return false;
            }
            polygon.push_back(corner);
        }

        if (polygon.size() < 3)
        {
            msg_error() << "OBJ face at line " << lineNumber << " has fewer than three vertices.";
            return false;
        }

        for (sofa::Size corner = 1; corner + 1 < polygon.size(); ++corner)
        {
            Triangle triangle;
            const std::array<FaceCorner, 3> fan { polygon[0], polygon[corner], polygon[corner + 1] };
            for (sofa::Size local = 0; local < 3; ++local)
            {
                triangle.vertex[local] = fan[local].vertex;
                triangle.normal[local] = fan[local].normal;
            }
            m_triangles.push_back(triangle);
        }
    }

    if (m_vertices.empty() || m_triangles.empty())
    {
        msg_error() << "OBJ contact mesh contains no usable vertices or faces: " << filename;
        return false;
    }

    return true;
}

template<class T1, class T2>
bool MeshNCPContactForceField<T1, T2>::finalizeMesh()
{
    sofa::type::vector<Triangle> validTriangles;
    validTriangles.reserve(m_triangles.size());
    m_generatedVertexNormals.assign(m_vertices.size(), Vec3(Real(0), Real(0), Real(0)));

    for (Triangle triangle : m_triangles)
    {
        const Vec3& a = m_vertices[triangle.vertex[0]];
        const Vec3& b = m_vertices[triangle.vertex[1]];
        const Vec3& c = m_vertices[triangle.vertex[2]];
        const Vec3 areaNormal = cross(b - a, c - a);

        triangle.faceNormal = areaNormal;
        if (!normalize(triangle.faceNormal))
            continue;

        for (sofa::Size axis = 0; axis < 3; ++axis)
        {
            triangle.minimum[axis] = std::min({ a[axis], b[axis], c[axis] });
            triangle.maximum[axis] = std::max({ a[axis], b[axis], c[axis] });
            triangle.centroid[axis] = (a[axis] + b[axis] + c[axis]) / Real(3);
        }

        m_generatedVertexNormals[triangle.vertex[0]] += areaNormal;
        m_generatedVertexNormals[triangle.vertex[1]] += areaNormal;
        m_generatedVertexNormals[triangle.vertex[2]] += areaNormal;
        validTriangles.push_back(triangle);
    }

    m_triangles = std::move(validTriangles);
    if (m_triangles.empty())
    {
        msg_error() << "All OBJ triangles are degenerate.";
        return false;
    }

    for (Vec3& normal : m_generatedVertexNormals)
        normalize(normal);

    using Edge = std::pair<sofa::Index, sofa::Index>;
    std::map<Edge, unsigned int> edgeCount;

    auto orderedEdge = [](sofa::Index a, sofa::Index b)
    {
        return a < b ? Edge(a, b) : Edge(b, a);
    };

    for (const Triangle& triangle : m_triangles)
    {
        ++edgeCount[orderedEdge(triangle.vertex[1], triangle.vertex[2])];
        ++edgeCount[orderedEdge(triangle.vertex[2], triangle.vertex[0])];
        ++edgeCount[orderedEdge(triangle.vertex[0], triangle.vertex[1])];
    }

    sofa::Size boundaryCount = 0;
    for (const auto& [edge, count] : edgeCount)
    {
        (void)edge;
        if (count == 1)
            ++boundaryCount;
    }

    for (Triangle& triangle : m_triangles)
    {
        triangle.boundaryEdge[0] = edgeCount[orderedEdge(triangle.vertex[1], triangle.vertex[2])] == 1;
        triangle.boundaryEdge[1] = edgeCount[orderedEdge(triangle.vertex[2], triangle.vertex[0])] == 1;
        triangle.boundaryEdge[2] = edgeCount[orderedEdge(triangle.vertex[0], triangle.vertex[1])] == 1;
    }

    d_meshBoundaryEdgeCount.setValue(boundaryCount);
    return true;
}

template<class T1, class T2>
bool MeshNCPContactForceField<T1, T2>::buildBvh()
{
    m_triangleOrder.resize(m_triangles.size());
    for (sofa::Index triangle = 0; triangle < m_triangles.size(); ++triangle)
        m_triangleOrder[triangle] = triangle;

    m_bvhNodes.clear();
    m_bvhNodes.reserve(2 * m_triangles.size());
    buildBvhNode(0, m_triangleOrder.size());
    return !m_bvhNodes.empty();
}

template<class T1, class T2>
sofa::Index MeshNCPContactForceField<T1, T2>::buildBvhNode(sofa::Size begin, sofa::Size end)
{
    const Real largest = std::numeric_limits<Real>::max();
    const Real lowest = std::numeric_limits<Real>::lowest();
    Vec3 minimum(largest, largest, largest);
    Vec3 maximum(lowest, lowest, lowest);
    Vec3 centroidMinimum(largest, largest, largest);
    Vec3 centroidMaximum(lowest, lowest, lowest);

    for (sofa::Size orderIndex = begin; orderIndex < end; ++orderIndex)
    {
        const Triangle& triangle = m_triangles[m_triangleOrder[orderIndex]];
        for (sofa::Size axis = 0; axis < 3; ++axis)
        {
            minimum[axis] = std::min(minimum[axis], triangle.minimum[axis]);
            maximum[axis] = std::max(maximum[axis], triangle.maximum[axis]);
            centroidMinimum[axis] = std::min(centroidMinimum[axis], triangle.centroid[axis]);
            centroidMaximum[axis] = std::max(centroidMaximum[axis], triangle.centroid[axis]);
        }
    }

    const sofa::Index nodeIndex = m_bvhNodes.size();
    m_bvhNodes.emplace_back();
    m_bvhNodes[nodeIndex].minimum = minimum;
    m_bvhNodes[nodeIndex].maximum = maximum;

    const sofa::Size count = end - begin;
    if (count <= static_cast<sofa::Size>(d_bvhLeafSize.getValue()))
    {
        m_bvhNodes[nodeIndex].begin = begin;
        m_bvhNodes[nodeIndex].count = count;
        return nodeIndex;
    }

    sofa::Size splitAxis = 0;
    Real splitExtent = centroidMaximum[0] - centroidMinimum[0];
    for (sofa::Size axis = 1; axis < 3; ++axis)
    {
        const Real extent = centroidMaximum[axis] - centroidMinimum[axis];
        if (extent > splitExtent)
        {
            splitExtent = extent;
            splitAxis = axis;
        }
    }

    if (!(splitExtent > Real(1e-15)))
    {
        m_bvhNodes[nodeIndex].begin = begin;
        m_bvhNodes[nodeIndex].count = count;
        return nodeIndex;
    }

    const sofa::Size middle = begin + count / 2;
    std::nth_element(
        m_triangleOrder.begin() + begin,
        m_triangleOrder.begin() + middle,
        m_triangleOrder.begin() + end,
        [this, splitAxis](sofa::Index lhs, sofa::Index rhs)
        {
            return m_triangles[lhs].centroid[splitAxis] < m_triangles[rhs].centroid[splitAxis];
        });

    const sofa::Index left = buildBvhNode(begin, middle);
    const sofa::Index right = buildBvhNode(middle, end);
    m_bvhNodes[nodeIndex].left = left;
    m_bvhNodes[nodeIndex].right = right;
    return nodeIndex;
}

template<class T1, class T2>
typename MeshNCPContactForceField<T1, T2>::ContactStatus
MeshNCPContactForceField<T1, T2>::computeContactKinematics(const Vec3& position, Contact& contact) const
{
    if (isPinned(contact.pointIndex))
        return ContactStatus::Pinned;

    if (m_bvhNodes.empty())
        return ContactStatus::InvalidGeometry;

    ClosestPointResult result;
    const Real maximumDistance = d_maxSearchDistance.getValue();
    if (maximumDistance > Real(0))
        result.squaredDistance = maximumDistance * maximumDistance;

    if (!findClosestPoint(position, result))
        return ContactStatus::Pinned;

    const Triangle& triangle = m_triangles[result.triangle];
    if (liesOnIgnoredBoundary(triangle, result.barycentric))
        return ContactStatus::Pinned;

    const Vec3 normal = contactNormal(triangle, result.barycentric);
    if (!(normal.norm2() > Real(1e-30)))
        return ContactStatus::InvalidGeometry;

    contact.gapGradient = normal;
    contact.gap = dot(position - result.point, normal) - d_contactOffset.getValue();
    return ContactStatus::Active;
}

template<class T1, class T2>
bool MeshNCPContactForceField<T1, T2>::findClosestPoint(
    const Vec3& position, ClosestPointResult& result) const
{
    if (m_bvhNodes.empty())
        return false;

    queryBvhNode(0, position, result);
    return result.triangle != InvalidIndex;
}

template<class T1, class T2>
void MeshNCPContactForceField<T1, T2>::queryBvhNode(
    sofa::Index nodeIndex, const Vec3& position, ClosestPointResult& result) const
{
    const BvhNode& node = m_bvhNodes[nodeIndex];
    if (squaredDistanceToBox(position, node.minimum, node.maximum) > result.squaredDistance)
        return;

    if (node.isLeaf())
    {
        for (sofa::Size local = 0; local < node.count; ++local)
        {
            const sofa::Index triangleIndex = m_triangleOrder[node.begin + local];
            const Triangle& triangle = m_triangles[triangleIndex];
            Vec3 barycentric;
            const Vec3 closest = closestPointOnTriangle(
                position,
                m_vertices[triangle.vertex[0]],
                m_vertices[triangle.vertex[1]],
                m_vertices[triangle.vertex[2]],
                barycentric);
            const Real distance2 = (position - closest).norm2();

            if (distance2 < result.squaredDistance)
            {
                result.triangle = triangleIndex;
                result.point = closest;
                result.barycentric = barycentric;
                result.squaredDistance = distance2;
            }
        }
        return;
    }

    const Real leftDistance = squaredDistanceToBox(
        position, m_bvhNodes[node.left].minimum, m_bvhNodes[node.left].maximum);
    const Real rightDistance = squaredDistanceToBox(
        position, m_bvhNodes[node.right].minimum, m_bvhNodes[node.right].maximum);

    if (leftDistance <= rightDistance)
    {
        queryBvhNode(node.left, position, result);
        queryBvhNode(node.right, position, result);
    }
    else
    {
        queryBvhNode(node.right, position, result);
        queryBvhNode(node.left, position, result);
    }
}

template<class T1, class T2>
typename MeshNCPContactForceField<T1, T2>::Vec3
MeshNCPContactForceField<T1, T2>::contactNormal(
    const Triangle& triangle, const Vec3& barycentric) const
{
    Vec3 normal = triangle.faceNormal;

    if (d_normalMode.getValue() == static_cast<unsigned int>(InterpolatedVertexNormal))
    {
        Vec3 interpolated(Real(0), Real(0), Real(0));
        for (sofa::Size corner = 0; corner < 3; ++corner)
        {
            Vec3 cornerNormal;
            if (triangle.normal[corner] != InvalidIndex
                && triangle.normal[corner] < m_objNormals.size()
                && m_objNormals[triangle.normal[corner]].norm2() > Real(1e-30))
            {
                cornerNormal = m_objNormals[triangle.normal[corner]];
            }
            else
            {
                cornerNormal = m_generatedVertexNormals[triangle.vertex[corner]];
            }

            interpolated += cornerNormal * barycentric[corner];
        }

        if (normalize(interpolated))
        {
            if (dot(interpolated, triangle.faceNormal) < Real(0))
                interpolated = -interpolated;
            normal = interpolated;
        }
    }

    if (d_flipNormals.getValue())
        normal = -normal;

    return normal;
}

template<class T1, class T2>
bool MeshNCPContactForceField<T1, T2>::liesOnIgnoredBoundary(
    const Triangle& triangle, const Vec3& barycentric) const
{
    if (!d_ignoreBoundaryEdges.getValue())
        return false;

    const Real tolerance = d_boundaryBarycentricTolerance.getValue();
    for (sofa::Size oppositeCorner = 0; oppositeCorner < 3; ++oppositeCorner)
    {
        if (triangle.boundaryEdge[oppositeCorner]
            && barycentric[oppositeCorner] <= tolerance)
        {
            return true;
        }
    }

    return false;
}

template<class T1, class T2>
bool MeshNCPContactForceField<T1, T2>::isPinned(sofa::Index pointIndex) const
{
    const auto& pinned = d_pinnedIndices.getValue();
    return std::find(pinned.begin(), pinned.end(), static_cast<unsigned int>(pointIndex)) != pinned.end();
}

template<class T1, class T2>
bool MeshNCPContactForceField<T1, T2>::parseFaceCorner(
    const std::string& token, sofa::Size vertexCount, sofa::Size normalCount, FaceCorner& corner)
{
    const std::size_t firstSlash = token.find('/');
    const std::size_t secondSlash = firstSlash == std::string::npos
        ? std::string::npos
        : token.find('/', firstSlash + 1);

    const std::string vertexToken = firstSlash == std::string::npos
        ? token
        : token.substr(0, firstSlash);
    const std::string normalToken = secondSlash == std::string::npos
        ? std::string()
        : token.substr(secondSlash + 1);

    try
    {
        if (vertexToken.empty()
            || !resolveObjIndex(std::stoll(vertexToken), vertexCount, corner.vertex))
        {
            return false;
        }

        if (!normalToken.empty())
        {
            sofa::Index normalIndex = InvalidIndex;
            if (!resolveObjIndex(std::stoll(normalToken), normalCount, normalIndex))
                return false;
            corner.normal = normalIndex;
        }
    }
    catch (...)
    {
        return false;
    }

    return true;
}

template<class T1, class T2>
bool MeshNCPContactForceField<T1, T2>::resolveObjIndex(
    long long rawIndex, sofa::Size count, sofa::Index& resolvedIndex)
{
    if (rawIndex == 0 || count == 0)
        return false;

    const long long resolved = rawIndex > 0
        ? rawIndex - 1
        : static_cast<long long>(count) + rawIndex;

    if (resolved < 0 || resolved >= static_cast<long long>(count))
        return false;

    resolvedIndex = static_cast<sofa::Index>(resolved);
    return true;
}

template<class T1, class T2>
typename MeshNCPContactForceField<T1, T2>::Vec3
MeshNCPContactForceField<T1, T2>::closestPointOnTriangle(
    const Vec3& position, const Vec3& a, const Vec3& b, const Vec3& c, Vec3& barycentric)
{
    const Vec3 ab = b - a;
    const Vec3 ac = c - a;
    const Vec3 ap = position - a;
    const Real d1 = dot(ab, ap);
    const Real d2 = dot(ac, ap);

    if (d1 <= Real(0) && d2 <= Real(0))
    {
        barycentric = Vec3(Real(1), Real(0), Real(0));
        return a;
    }

    const Vec3 bp = position - b;
    const Real d3 = dot(ab, bp);
    const Real d4 = dot(ac, bp);
    if (d3 >= Real(0) && d4 <= d3)
    {
        barycentric = Vec3(Real(0), Real(1), Real(0));
        return b;
    }

    const Real vc = d1 * d4 - d3 * d2;
    if (vc <= Real(0) && d1 >= Real(0) && d3 <= Real(0))
    {
        const Real v = d1 / (d1 - d3);
        barycentric = Vec3(Real(1) - v, v, Real(0));
        return a + ab * v;
    }

    const Vec3 cp = position - c;
    const Real d5 = dot(ab, cp);
    const Real d6 = dot(ac, cp);
    if (d6 >= Real(0) && d5 <= d6)
    {
        barycentric = Vec3(Real(0), Real(0), Real(1));
        return c;
    }

    const Real vb = d5 * d2 - d1 * d6;
    if (vb <= Real(0) && d2 >= Real(0) && d6 <= Real(0))
    {
        const Real w = d2 / (d2 - d6);
        barycentric = Vec3(Real(1) - w, Real(0), w);
        return a + ac * w;
    }

    const Real va = d3 * d6 - d5 * d4;
    if (va <= Real(0) && (d4 - d3) >= Real(0) && (d5 - d6) >= Real(0))
    {
        const Real w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        barycentric = Vec3(Real(0), Real(1) - w, w);
        return b + (c - b) * w;
    }

    const Real denominator = Real(1) / (va + vb + vc);
    const Real v = vb * denominator;
    const Real w = vc * denominator;
    barycentric = Vec3(Real(1) - v - w, v, w);
    return a + ab * v + ac * w;
}

template<class T1, class T2>
typename MeshNCPContactForceField<T1, T2>::Real
MeshNCPContactForceField<T1, T2>::squaredDistanceToBox(
    const Vec3& position, const Vec3& minimum, const Vec3& maximum)
{
    Real distance2 = Real(0);
    for (sofa::Size axis = 0; axis < 3; ++axis)
    {
        if (position[axis] < minimum[axis])
        {
            const Real delta = minimum[axis] - position[axis];
            distance2 += delta * delta;
        }
        else if (position[axis] > maximum[axis])
        {
            const Real delta = position[axis] - maximum[axis];
            distance2 += delta * delta;
        }
    }
    return distance2;
}

template<class T1, class T2>
typename MeshNCPContactForceField<T1, T2>::Real
MeshNCPContactForceField<T1, T2>::dot(const Vec3& a, const Vec3& b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template<class T1, class T2>
typename MeshNCPContactForceField<T1, T2>::Vec3
MeshNCPContactForceField<T1, T2>::cross(const Vec3& a, const Vec3& b)
{
    return Vec3(
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]);
}

template<class T1, class T2>
bool MeshNCPContactForceField<T1, T2>::normalize(Vec3& value)
{
    const Real norm2 = dot(value, value);
    if (!std::isfinite(static_cast<double>(norm2)) || norm2 <= Real(1e-30))
        return false;

    value *= Real(1) / std::sqrt(norm2);
    return true;
}

} // namespace sofa::ncp
