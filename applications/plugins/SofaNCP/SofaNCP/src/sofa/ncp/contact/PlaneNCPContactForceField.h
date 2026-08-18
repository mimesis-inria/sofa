/****************************************************************************
* Infinite-plane geometry for the Fischer-Burmeister NCP contact base.
*
* Gap convention:
*   g(p) = (p - point) . n
*
* where n is the normalized plane normal pointing into the admissible
* half-space. Therefore:
*   g >= 0 : admissible side,
*   g <  0 : penetration,
*   dg/dp  = n,
*   d2g/dp2 = 0.
*
* This component only defines analytical contact geometry. FB residuals,
* multipliers, compliance scaling, force assembly, symmetric tangent assembly,
* line-search refreshes, and diagnostics remain in the generic base class.
****************************************************************************/
#pragma once

#include <sofa/ncp/config.h>
#include <sofa/ncp/contact/FischerBurmeisterContactForceField.h>

#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/type/RGBAColor.h>

namespace sofa::core { class ObjectFactory; }

namespace sofa::ncp
{

template<class TDataTypes1, class TDataTypes2>
class PlaneNCPContactForceField
    : public FischerBurmeisterContactForceField<TDataTypes1, TDataTypes2>
{
public:
    SOFA_CLASS(
        SOFA_TEMPLATE2(PlaneNCPContactForceField, TDataTypes1, TDataTypes2),
        SOFA_TEMPLATE2(FischerBurmeisterContactForceField, TDataTypes1, TDataTypes2)
    );

    using Base = FischerBurmeisterContactForceField<TDataTypes1, TDataTypes2>;
    using DataTypes1 = TDataTypes1;
    using DataTypes2 = TDataTypes2;
    using Real = typename Base::Real;
    using Vec3 = typename Base::Vec3;
    using Contact = typename Base::Contact;
    using ContactStatus = typename Base::ContactStatus;

    /// Any point lying on the plane.
    Data<Vec3> d_point;

    /// Normal pointing into the admissible half-space. Normalized internally.
    Data<Vec3> d_normal;

    /// Visualization only: draw a finite square representing the infinite plane.
    Data<bool> d_showPlane;
    Data<Real> d_drawSize;
    Data<sofa::type::RGBAColor> d_color;

protected:
    PlaneNCPContactForceField();

    PlaneNCPContactForceField(
        core::behavior::MechanicalState<DataTypes1>* object1,
        core::behavior::MechanicalState<DataTypes2>* object2);

public:
    void draw(const core::visual::VisualParams* vparams) override;

protected:
    ContactStatus computeContactKinematics(
        const Vec3& p,
        Contact& contact) const override;

    bool tryGetUnitNormal(Vec3& normal) const;
};

#if !defined(SOFANCP_PLANE_CONTACT_FORCE_FIELD_CPP)
extern template class SOFANCP_API
    PlaneNCPContactForceField<defaulttype::Vec3Types, defaulttype::Vec1Types>;

extern template class SOFANCP_API
    PlaneNCPContactForceField<defaulttype::Rigid3Types, defaulttype::Vec1Types>;
#endif

SOFANCP_API void registerPlaneNCPContactForceField(sofa::core::ObjectFactory* factory);

} // namespace sofa::ncp
