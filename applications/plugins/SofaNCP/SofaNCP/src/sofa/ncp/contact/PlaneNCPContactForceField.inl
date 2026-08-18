/****************************************************************************
* Infinite-plane geometry for the Fischer-Burmeister NCP contact base.
****************************************************************************/
#pragma once

#include <sofa/ncp/contact/PlaneNCPContactForceField.h>
#include <sofa/ncp/contact/FischerBurmeisterContactForceField.inl>
#include <sofa/core/visual/VisualParams.h>

#include <cmath>
#include <vector>

namespace sofa::ncp
{

template<class TDataTypes1, class TDataTypes2>
PlaneNCPContactForceField<TDataTypes1, TDataTypes2>::
PlaneNCPContactForceField()
    : Base()
    , d_point(this->initData(
          &d_point,
          Vec3(Real(0), Real(0), Real(0)),
          "point",
          "Any point on the infinite plane."))
    , d_normal(this->initData(
          &d_normal,
          Vec3(Real(0), Real(1), Real(0)),
          "normal",
          "Plane normal pointing into the admissible half-space. Normalized internally."))
    , d_showPlane(this->initData(
          &d_showPlane,
          true,
          "showPlane",
          "Draw a finite square representing the infinite plane."))
    , d_drawSize(this->initData(
          &d_drawSize,
          Real(50),
          "drawSize",
          "Half-size of the square used to draw the plane."))
    , d_color(this->initData(
          &d_color,
          sofa::type::RGBAColor(0.20f, 0.80f, 0.30f, 0.45f),
          "color",
          "Plane wireframe and normal color."))
{
}

template<class TDataTypes1, class TDataTypes2>
PlaneNCPContactForceField<TDataTypes1, TDataTypes2>::
PlaneNCPContactForceField(
    core::behavior::MechanicalState<DataTypes1>* object1,
    core::behavior::MechanicalState<DataTypes2>* object2)
    : PlaneNCPContactForceField()
{
    this->mstate1.set(object1);
    this->mstate2.set(object2);
}

template<class TDataTypes1, class TDataTypes2>
bool PlaneNCPContactForceField<TDataTypes1, TDataTypes2>::
tryGetUnitNormal(Vec3& normal) const
{
    normal = d_normal.getValue();

    const Real norm2 = normal.norm2();
    constexpr Real epsilon = Real(1e-12);

    if (!std::isfinite(static_cast<double>(normal[0])) ||
        !std::isfinite(static_cast<double>(normal[1])) ||
        !std::isfinite(static_cast<double>(normal[2])) ||
        !std::isfinite(static_cast<double>(norm2)) ||
        norm2 <= epsilon * epsilon)
    {
        return false;
    }

    normal *= Real(1) / std::sqrt(norm2);
    return true;
}

template<class TDataTypes1, class TDataTypes2>
typename PlaneNCPContactForceField<TDataTypes1, TDataTypes2>::ContactStatus
PlaneNCPContactForceField<TDataTypes1, TDataTypes2>::
computeContactKinematics(
    const Vec3& p,
    Contact& contact) const
{
    const Vec3 point = d_point.getValue();

    const auto finiteVec = [](const Vec3& value)
    {
        return
            std::isfinite(static_cast<double>(value[0])) &&
            std::isfinite(static_cast<double>(value[1])) &&
            std::isfinite(static_cast<double>(value[2]));
    };

    Vec3 normal;

    if (!finiteVec(p) ||
        !finiteVec(point) ||
        !tryGetUnitNormal(normal))
    {
        return ContactStatus::InvalidGeometry;
    }

    /*
     * The plane is infinite, so every valid point has a well-defined NCP row.
     * Far points do not need to be pinned: for g>0 and lambda=0, the exact FB
     * residual is zero and the mechanical coupling coefficient dPhiDgap is zero.
     *
     * Because the plane gap is affine:
     *   g(p)       = (p - point) . normal,
     *   grad g     = normal,
     *   Hessian(g) = 0.
     */
    contact.gap = (p - point) * normal;
    contact.gapGradient = normal;

    if (!std::isfinite(static_cast<double>(contact.gap)))
        return ContactStatus::InvalidGeometry;

    return ContactStatus::Active;
}

template<class TDataTypes1, class TDataTypes2>
void PlaneNCPContactForceField<TDataTypes1, TDataTypes2>::draw(
    const core::visual::VisualParams* vparams)
{
    // Preserve visualization of contact gradients supplied by the FB base.
    Base::draw(vparams);

    if (!vparams ||
        !vparams->drawTool() ||
        !d_showPlane.getValue())
    {
        return;
    }

    const Vec3 point = d_point.getValue();
    const Real size = d_drawSize.getValue();
    Vec3 normal;

    if (!tryGetUnitNormal(normal) ||
        !std::isfinite(static_cast<double>(size)) ||
        size <= Real(0))
    {
        return;
    }

    // Construct a stable orthonormal basis spanning the plane.
    const Vec3 reference =
        std::abs(normal[0]) < Real(0.9)
            ? Vec3(Real(1), Real(0), Real(0))
            : Vec3(Real(0), Real(1), Real(0));

    Vec3 tangent1 = reference - normal * (reference * normal);
    const Real tangentNorm2 = tangent1.norm2();
    constexpr Real epsilon = Real(1e-12);

    if (!std::isfinite(static_cast<double>(tangentNorm2)) ||
        tangentNorm2 <= epsilon * epsilon)
    {
        return;
    }

    tangent1 *= Real(1) / std::sqrt(tangentNorm2);

    const Vec3 tangent2(
        normal[1] * tangent1[2] - normal[2] * tangent1[1],
        normal[2] * tangent1[0] - normal[0] * tangent1[2],
        normal[0] * tangent1[1] - normal[1] * tangent1[0]);

    const Vec3 c0 = point - tangent1 * size - tangent2 * size;
    const Vec3 c1 = point + tangent1 * size - tangent2 * size;
    const Vec3 c2 = point + tangent1 * size + tangent2 * size;
    const Vec3 c3 = point - tangent1 * size + tangent2 * size;

    const auto toDrawVec = [](const Vec3& value)
    {
        return sofa::type::Vec3(
            static_cast<float>(value[0]),
            static_cast<float>(value[1]),
            static_cast<float>(value[2]));
    };

    std::vector<sofa::type::Vec3> lines;
    lines.reserve(12);

    // Square boundary.
    lines.emplace_back(toDrawVec(c0));
    lines.emplace_back(toDrawVec(c1));
    lines.emplace_back(toDrawVec(c1));
    lines.emplace_back(toDrawVec(c2));
    lines.emplace_back(toDrawVec(c2));
    lines.emplace_back(toDrawVec(c3));
    lines.emplace_back(toDrawVec(c3));
    lines.emplace_back(toDrawVec(c0));

    // Two center axes plus a normal arrow showing the admissible side.
    lines.emplace_back(toDrawVec(point - tangent1 * size));
    lines.emplace_back(toDrawVec(point + tangent1 * size));
    lines.emplace_back(toDrawVec(point));
    lines.emplace_back(toDrawVec(point + normal * (Real(0.25) * size)));

    vparams->drawTool()->drawLines(
        lines,
        1.5f,
        d_color.getValue());
}

} // namespace sofa::ncp
