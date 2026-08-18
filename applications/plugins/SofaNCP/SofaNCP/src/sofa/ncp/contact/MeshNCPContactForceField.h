/****************************************************************************
* Triangle-mesh specialization of FischerBurmeisterContactForceField.
*
* The mesh is loaded directly from an OBJ file and queried through a static
* AABB tree. The contact gap is the oriented tangent-plane distance
*
*     g(x) = (x - q) . n - contactOffset,
*
* where q is the closest point on the selected triangle and n points toward
* the admissible lumen. Set flipNormals=true when the OBJ winding points in
* the opposite direction.
****************************************************************************/
#pragma once

#include <sofa/ncp/config.h>
#include <sofa/ncp/contact/FischerBurmeisterContactForceField.h>

#include <sofa/core/objectmodel/DataFileName.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/type/Vec.h>
#include <sofa/type/vector.h>

#include <array>
#include <limits>
#include <string>

namespace sofa::core { class ObjectFactory; }

namespace sofa::ncp
{

template<class TDataTypes1, class TDataTypes2>
class MeshNCPContactForceField final
    : public FischerBurmeisterContactForceField<TDataTypes1, TDataTypes2>
{
public:
    SOFA_CLASS(
        SOFA_TEMPLATE2(MeshNCPContactForceField, TDataTypes1, TDataTypes2),
        SOFA_TEMPLATE2(FischerBurmeisterContactForceField, TDataTypes1, TDataTypes2));

    using Inherit = FischerBurmeisterContactForceField<TDataTypes1, TDataTypes2>;
    using Real = typename Inherit::Real;
    using Vec3 = typename Inherit::Vec3;
    using Contact = typename Inherit::Contact;
    using ContactStatus = typename Inherit::ContactStatus;

    enum NormalMode : unsigned int
    {
        FaceNormal = 0,
        InterpolatedVertexNormal = 1
    };

    sofa::core::objectmodel::DataFileName d_meshFilename;
    Data<bool> d_flipNormals;
    Data<unsigned int> d_normalMode;
    Data<Real> d_contactOffset;
    Data<Real> d_maxSearchDistance;
    Data<unsigned int> d_bvhLeafSize;
    Data<bool> d_ignoreBoundaryEdges;
    Data<Real> d_boundaryBarycentricTolerance;
    Data<sofa::type::vector<unsigned int>> d_pinnedIndices;

    Data<sofa::Size> d_meshVertexCount;
    Data<sofa::Size> d_meshTriangleCount;
    Data<sofa::Size> d_meshBoundaryEdgeCount;

protected:
    MeshNCPContactForceField();

public:
    void init() override;
    void reinit() override;

protected:
    ContactStatus computeContactKinematics(const Vec3& position, Contact& contact) const override;

private:
    static constexpr sofa::Index InvalidIndex = std::numeric_limits<sofa::Index>::max();

    struct FaceCorner
    {
        sofa::Index vertex = InvalidIndex;
        sofa::Index normal = InvalidIndex;
    };

    struct Triangle
    {
        std::array<sofa::Index, 3> vertex { InvalidIndex, InvalidIndex, InvalidIndex };
        std::array<sofa::Index, 3> normal { InvalidIndex, InvalidIndex, InvalidIndex };
        std::array<bool, 3> boundaryEdge { false, false, false };
        Vec3 faceNormal = Vec3(Real(0), Real(0), Real(0));
        Vec3 minimum = Vec3(Real(0), Real(0), Real(0));
        Vec3 maximum = Vec3(Real(0), Real(0), Real(0));
        Vec3 centroid = Vec3(Real(0), Real(0), Real(0));
    };

    struct BvhNode
    {
        Vec3 minimum = Vec3(Real(0), Real(0), Real(0));
        Vec3 maximum = Vec3(Real(0), Real(0), Real(0));
        sofa::Index left = InvalidIndex;
        sofa::Index right = InvalidIndex;
        sofa::Size begin = 0;
        sofa::Size count = 0;

        bool isLeaf() const { return count != 0; }
    };

    struct ClosestPointResult
    {
        sofa::Index triangle = InvalidIndex;
        Vec3 point = Vec3(Real(0), Real(0), Real(0));
        Vec3 barycentric = Vec3(Real(0), Real(0), Real(0));
        Real squaredDistance = std::numeric_limits<Real>::max();
    };

    sofa::type::vector<Vec3> m_vertices;
    sofa::type::vector<Vec3> m_objNormals;
    sofa::type::vector<Vec3> m_generatedVertexNormals;
    sofa::type::vector<Triangle> m_triangles;
    sofa::type::vector<sofa::Index> m_triangleOrder;
    sofa::type::vector<BvhNode> m_bvhNodes;

    bool loadMesh();
    bool loadObj(const std::string& filename);
    bool finalizeMesh();
    bool buildBvh();
    sofa::Index buildBvhNode(sofa::Size begin, sofa::Size end);

    bool findClosestPoint(const Vec3& position, ClosestPointResult& result) const;
    void queryBvhNode(sofa::Index nodeIndex, const Vec3& position, ClosestPointResult& result) const;
    Vec3 contactNormal(const Triangle& triangle, const Vec3& barycentric) const;
    bool liesOnIgnoredBoundary(const Triangle& triangle, const Vec3& barycentric) const;
    bool isPinned(sofa::Index pointIndex) const;

    static bool parseFaceCorner(const std::string& token, sofa::Size vertexCount,
        sofa::Size normalCount, FaceCorner& corner);
    static bool resolveObjIndex(long long rawIndex, sofa::Size count, sofa::Index& resolvedIndex);
    static Vec3 closestPointOnTriangle(const Vec3& position, const Vec3& a,
        const Vec3& b, const Vec3& c, Vec3& barycentric);
    static Real squaredDistanceToBox(const Vec3& position, const Vec3& minimum, const Vec3& maximum);
    static Real dot(const Vec3& a, const Vec3& b);
    static Vec3 cross(const Vec3& a, const Vec3& b);
    static bool normalize(Vec3& value);
};

#if !defined(SOFANCP_MESH_NCP_CONTACT_FORCE_FIELD_CPP)
extern template class SOFANCP_API MeshNCPContactForceField<defaulttype::Rigid3Types, defaulttype::Vec1Types>;
extern template class SOFANCP_API MeshNCPContactForceField<defaulttype::Vec3Types, defaulttype::Vec1Types>;
#endif

SOFANCP_API void registerMeshNCPContactForceField(sofa::core::ObjectFactory* factory);

} // namespace sofa::ncp
