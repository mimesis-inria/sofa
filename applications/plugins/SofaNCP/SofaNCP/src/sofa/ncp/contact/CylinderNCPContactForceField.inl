/****************************************************************************
* Finite-cylinder geometry for the Fischer-Burmeister NCP contact base.
****************************************************************************/
#pragma once

#include <sofa/ncp/contact/CylinderNCPContactForceField.h>
#include <sofa/ncp/contact/FischerBurmeisterContactForceField.inl>
#include <sofa/core/visual/VisualParams.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace sofa::ncp
{

template<class TDataTypes1, class TDataTypes2>
CylinderNCPContactForceField<TDataTypes1, TDataTypes2>::
CylinderNCPContactForceField()
    : Base()
    , d_center(this->initData( &d_center, Vec3(Real(0), Real(0), Real(0)), "center", "Cylinder center."))
    , d_axis(this->initData( &d_axis, Vec3(Real(1), Real(0), Real(0)), "axis", "Cylinder axis direction."))
    , d_radius(this->initData( &d_radius, Real(10), "radius", "Cylinder radius."))
    , d_length(this->initData( &d_length, Real(100), "length", "Cylinder length."))
    , d_endCaps(this->initData( &d_endCaps, static_cast<unsigned int>(NoCaps), "endCaps", "0=none, 1=negative-axis end, 2=positive-axis end, 3=both."))
    , d_showCylinder(this->initData(&d_showCylinder, true, "showCylinder", "Draw the cylinder wireframe."))
    , d_color(this->initData( &d_color, sofa::type::RGBAColor(0.05f, 0.55f, 1.0f, 0.70f), "color", "Cylinder wireframe color."))
{
}

template<class TDataTypes1, class TDataTypes2>
CylinderNCPContactForceField<TDataTypes1, TDataTypes2>::
CylinderNCPContactForceField(
    core::behavior::MechanicalState<DataTypes1>* object1,
    core::behavior::MechanicalState<DataTypes2>* object2)
    : CylinderNCPContactForceField()
{
    this->mstate1.set(object1);
    this->mstate2.set(object2);
}

template<class TDataTypes1, class TDataTypes2>
bool CylinderNCPContactForceField<TDataTypes1, TDataTypes2>::tryGetUnitAxis(Vec3& axis) const
{
    axis = d_axis.getValue();

    const Real norm2 = axis.norm2();
    constexpr Real epsilon = Real(1e-12);

    if (!std::isfinite(static_cast<double>(norm2)) ||
        norm2 <= epsilon * epsilon)
    {
        return false;
    }

    axis *= Real(1) / std::sqrt(norm2);
    return true;
}

template<class TDataTypes1, class TDataTypes2>
typename CylinderNCPContactForceField<TDataTypes1, TDataTypes2>::ContactStatus
CylinderNCPContactForceField<TDataTypes1, TDataTypes2>::computeContactKinematics(const Vec3& p,Contact& contact) const
{
    constexpr Real epsilon = Real(1e-12);

    const Vec3 center = d_center.getValue();
    const Real radius = d_radius.getValue();
    const Real length = d_length.getValue();
    const unsigned int capMode = d_endCaps.getValue();

    auto finiteVec = [](const Vec3& value)
    {
        return
            std::isfinite(static_cast<double>(value[0])) &&
            std::isfinite(static_cast<double>(value[1])) &&
            std::isfinite(static_cast<double>(value[2]));
    };

    Vec3 axis;

    if (!finiteVec(p) ||
        !finiteVec(center) ||
        !tryGetUnitAxis(axis) ||
        !std::isfinite(static_cast<double>(radius)) ||
        !std::isfinite(static_cast<double>(length)) ||
        radius <= epsilon ||
        length <= epsilon ||
        capMode > static_cast<unsigned int>(BothEndCaps))
    {
        return ContactStatus::InvalidGeometry;
    }

    const Real halfLength = Real(0.5) * length;
    const Vec3 relative = p - center;
    const Real axial = relative * axis;

    const bool negativeCap =
        (capMode & static_cast<unsigned int>(NegativeEndCap)) != 0;

    const bool positiveCap =
        (capMode & static_cast<unsigned int>(PositiveEndCap)) != 0;

    // Beyond an open end, this cylinder specialization defines no contact row.
    if ((axial < -halfLength && !negativeCap) ||
        (axial >  halfLength && !positiveCap))
    {
        return ContactStatus::Pinned;
    }

    const Vec3 radial = relative - axis * axial;
    const Real radialNorm2 = radial.norm2();

    if (!std::isfinite(static_cast<double>(radialNorm2)))
        return ContactStatus::InvalidGeometry;

    // The side-wall normal is undefined on the cylinder axis.
    if (radialNorm2 <= epsilon * epsilon)
        return ContactStatus::Pinned;

    const Real radialNorm = std::sqrt(radialNorm2);

    Real gap = radius - radialNorm;
    Vec3 gradient = radial * (Real(-1) / radialNorm);

    bool ambiguous = false;

    auto considerBoundary =
        [&](Real candidateGap, const Vec3& candidateGradient)
        {
            const Real scale = std::max(
                Real(1),
                std::max(std::abs(gap), std::abs(candidateGap)));

            const Real tolerance = epsilon * scale;

            if (candidateGap < gap - tolerance)
            {
                gap = candidateGap;
                gradient = candidateGradient;
                ambiguous = false;
            }
            else if (std::abs(candidateGap - gap) <= tolerance)
            {
                ambiguous = true;
            }
        };

    if (negativeCap)
        considerBoundary(halfLength + axial, axis);

    if (positiveCap)
        considerBoundary(halfLength - axial, axis * Real(-1));

    if (ambiguous ||
        !std::isfinite(static_cast<double>(gap)) ||
        !finiteVec(gradient))
    {
        return ContactStatus::InvalidGeometry;
    }

    contact.gap = gap;
    contact.gapGradient = gradient;
    return ContactStatus::Active;
}

// template<class TDataTypes1, class TDataTypes2>
// bool CylinderNCPContactForceField<TDataTypes1, TDataTypes2>::computeGapHessian(const Vec3& p,Mat3& hessian) const
// {
//     constexpr Real epsilon = Real(1e-12);

//     hessian.clear();

//     const Vec3 center = d_center.getValue();
//     const Real radius = d_radius.getValue();
//     const Real length = d_length.getValue();
//     const unsigned int capMode = d_endCaps.getValue();

//     auto finiteVec = [](const Vec3& value)
//     {
//         return
//             std::isfinite(static_cast<double>(value[0])) &&
//             std::isfinite(static_cast<double>(value[1])) &&
//             std::isfinite(static_cast<double>(value[2]));
//     };

//     Vec3 axis;

//     if (!finiteVec(p) ||
//         !finiteVec(center) ||
//         !tryGetUnitAxis(axis) ||
//         !std::isfinite(static_cast<double>(radius)) ||
//         !std::isfinite(static_cast<double>(length)) ||
//         radius <= epsilon ||
//         length <= epsilon ||
//         capMode > static_cast<unsigned int>(BothEndCaps))
//     {
//         return false;
//     }

//     const Real halfLength = Real(0.5) * length;
//     const Vec3 relative = p - center;
//     const Real axial = relative * axis;

//     const bool negativeCap =
//         (capMode & static_cast<unsigned int>(NegativeEndCap)) != 0;

//     const bool positiveCap =
//         (capMode & static_cast<unsigned int>(PositiveEndCap)) != 0;

//     // Same semantics as computeContactKinematics():
//     // beyond an open end there is no cylinder contact row.
//     if ((axial < -halfLength && !negativeCap) ||
//         (axial >  halfLength && !positiveCap))
//     {
//         return false;
//     }

//     const Vec3 radial = relative - axis * axial;
//     const Real radialNorm2 = radial.norm2();

//     if (!std::isfinite(static_cast<double>(radialNorm2)) ||
//         radialNorm2 <= epsilon * epsilon)
//     {
//         return false;
//     }

//     const Real radialNorm = std::sqrt(radialNorm2);

//     // Initially assume the side wall is the active boundary:
//     //
//     //     g_side = R - rho
//     //
//     // with
//     //
//     //     rho = ||P (p-c)||
//     //     P   = I - u u^T.
//     Real selectedGap = radius - radialNorm;

//     enum class Boundary
//     {
//         Side,
//         NegativeCap,
//         PositiveCap
//     };

//     Boundary selectedBoundary = Boundary::Side;
//     bool ambiguous = false;

//     auto considerBoundary =
//         [&](Real candidateGap, Boundary candidateBoundary)
//         {
//             const Real scale = std::max(
//                 Real(1),
//                 std::max(
//                     std::abs(selectedGap),
//                     std::abs(candidateGap)));

//             const Real tolerance = epsilon * scale;

//             if (candidateGap < selectedGap - tolerance)
//             {
//                 selectedGap = candidateGap;
//                 selectedBoundary = candidateBoundary;
//                 ambiguous = false;
//             }
//             else if (std::abs(candidateGap - selectedGap) <= tolerance)
//             {
//                 ambiguous = true;
//             }
//         };

//     if (negativeCap)
//     {
//         considerBoundary(
//             halfLength + axial,
//             Boundary::NegativeCap);
//     }

//     if (positiveCap)
//     {
//         considerBoundary(
//             halfLength - axial,
//             Boundary::PositiveCap);
//     }

//     // At an edge where two gap functions are equally active, the min-gap
//     // function is nonsmooth and there is no unique classical Hessian.
//     if (ambiguous)
//         return false;

//     // End-cap gaps are affine:
//     //
//     //     g_- = L/2 + axial
//     //     g_+ = L/2 - axial
//     //
//     // therefore their Hessian is exactly zero.
//     if (selectedBoundary != Boundary::Side)
//     {
//         hessian.clear();
//         return true;
//     }

//     // Side-wall gap:
//     //
//     //     g = R - rho
//     //     rho = ||P(p-c)||
//     //     P = I - u u^T
//     //     n = radial / rho
//     //
//     // therefore:
//     //
//     //     Hess(g) = -(P - n n^T) / rho.
//     const Vec3 n = radial * (Real(1) / radialNorm);
//     const Real inverseRadius = Real(1) / radialNorm;

//     for (sofa::Size i = 0; i < 3; ++i)
//     {
//         for (sofa::Size j = 0; j < 3; ++j)
//         {
//             const Real identity = i == j ? Real(1) : Real(0);

//             const Real projection =
//                 identity - axis[i] * axis[j];

//             hessian(i, j) =
//                 -(projection - n[i] * n[j])
//                 * inverseRadius;
//         }
//     }

//     return true;
// }

template<class TDataTypes1, class TDataTypes2>
void CylinderNCPContactForceField<TDataTypes1, TDataTypes2>::draw(const core::visual::VisualParams* vparams)
{
    Base::draw(vparams);

    if (!vparams || !vparams->drawTool() || !d_showCylinder.getValue())
        return;

    constexpr Real epsilon = Real(1e-12);
    constexpr Real twoPi =
        Real(6.283185307179586476925286766559);

    constexpr unsigned int radialSegments = 48;
    constexpr unsigned int longitudinalGuides = 8;

    const Vec3 center = d_center.getValue();
    const Real radius = d_radius.getValue();
    const Real length = d_length.getValue();
    const unsigned int capMode = d_endCaps.getValue();

    Vec3 axis;

    if (!tryGetUnitAxis(axis) ||
        !std::isfinite(static_cast<double>(radius)) ||
        !std::isfinite(static_cast<double>(length)) ||
        radius <= epsilon ||
        length <= epsilon ||
        capMode > static_cast<unsigned int>(BothEndCaps))
    {
        return;
    }

    const Real halfLength = Real(0.5) * length;

    const Vec3 negativeCenter = center - axis * halfLength;
    const Vec3 middleCenter = center;
    const Vec3 positiveCenter = center + axis * halfLength;

    const Vec3 reference =
        std::abs(axis[0]) < Real(0.9)
            ? Vec3(Real(1), Real(0), Real(0))
            : Vec3(Real(0), Real(1), Real(0));

    Vec3 e1 = reference - axis * (reference * axis);

    const Real e1Norm2 = e1.norm2();

    if (!std::isfinite(static_cast<double>(e1Norm2)) ||
        e1Norm2 <= epsilon * epsilon)
    {
        return;
    }

    e1 *= Real(1) / std::sqrt(e1Norm2);

    const Vec3 e2(
        axis[1] * e1[2] - axis[2] * e1[1],
        axis[2] * e1[0] - axis[0] * e1[2],
        axis[0] * e1[1] - axis[1] * e1[0]);

    const auto toDrawPoint = [](const Vec3& point)
    {
        return sofa::type::Vec3(
            static_cast<float>(point[0]),
            static_cast<float>(point[1]),
            static_cast<float>(point[2]));
    };

    const auto radialDirection =
        [&](unsigned int segment)
        {
            const Real angle =
                twoPi * Real(segment) / Real(radialSegments);

            return
                e1 * std::cos(angle) +
                e2 * std::sin(angle);
        };

    std::vector<sofa::type::Vec3> lines;
    lines.reserve(512);

    const auto appendCircle =
        [&](const Vec3& circleCenter, Real circleRadius)
        {
            for (unsigned int i = 0; i < radialSegments; ++i)
            {
                const Vec3 p0 =
                    circleCenter +
                    radialDirection(i) * circleRadius;

                const Vec3 p1 =
                    circleCenter +
                    radialDirection((i + 1) % radialSegments) *
                        circleRadius;

                lines.emplace_back(toDrawPoint(p0));
                lines.emplace_back(toDrawPoint(p1));
            }
        };

    // Outer shape.
    appendCircle(negativeCenter, radius);
    appendCircle(middleCenter, radius);
    appendCircle(positiveCenter, radius);

    // Sparse longitudinal guides preserve transparency.
    for (unsigned int i = 0; i < longitudinalGuides; ++i)
    {
        const unsigned int segment =
            i * radialSegments / longitudinalGuides;

        const Vec3 radial = radialDirection(segment);

        lines.emplace_back(
            toDrawPoint(negativeCenter + radial * radius));

        lines.emplace_back(
            toDrawPoint(positiveCenter + radial * radius));
    }

    // A second, smaller ring marks a closed end without obscuring the view.
    constexpr Real capIndicatorScale = Real(0.72);

    if ((capMode &
         static_cast<unsigned int>(NegativeEndCap)) != 0)
    {
        appendCircle(
            negativeCenter,
            radius * capIndicatorScale);
    }

    if ((capMode &
         static_cast<unsigned int>(PositiveEndCap)) != 0)
    {
        appendCircle(
            positiveCenter,
            radius * capIndicatorScale);
    }

    if (!lines.empty())
    {
        vparams->drawTool()->drawLines(
            lines,
            1.0f,
            d_color.getValue());
    }
}

} // namespace sofa::ncp