/****************************************************************************
* Finite-cylinder geometry for the Fischer-Burmeister NCP contact base.
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
class CylinderNCPContactForceField
    : public FischerBurmeisterContactForceField<TDataTypes1, TDataTypes2>
{
public:
    SOFA_CLASS(
        SOFA_TEMPLATE2(CylinderNCPContactForceField, TDataTypes1, TDataTypes2),
        SOFA_TEMPLATE2(FischerBurmeisterContactForceField, TDataTypes1, TDataTypes2)
    );

    using Base = FischerBurmeisterContactForceField<TDataTypes1, TDataTypes2>;
    using DataTypes1 = TDataTypes1;
    using DataTypes2 = TDataTypes2;
    using Real = typename Base::Real;
    using Vec3 = typename Base::Vec3;
    using Mat3 = typename Base::Mat3;
    using Contact = typename Base::Contact;
    using ContactStatus = typename Base::ContactStatus;

    enum EndCapMode : unsigned int
    {
        NoCaps = 0,
        NegativeEndCap = 1,
        PositiveEndCap = 2,
        BothEndCaps = 3
    };

    Data<Vec3> d_center;
    Data<Vec3> d_axis;
    Data<Real> d_radius;
    Data<Real> d_length;
    Data<unsigned int> d_endCaps;
    Data<bool> d_showCylinder;
    Data<sofa::type::RGBAColor> d_color;

protected:
    CylinderNCPContactForceField();

    CylinderNCPContactForceField(
        core::behavior::MechanicalState<DataTypes1>* object1,
        core::behavior::MechanicalState<DataTypes2>* object2);

public:
    void draw(const core::visual::VisualParams* vparams) override;

protected:
    ContactStatus computeContactKinematics(
        const Vec3& p,
        Contact& contact) const override;

    // bool computeGapHessian(
    //     const Vec3& p,
    //     Mat3& hessian) const override;

    bool tryGetUnitAxis(Vec3& axis) const;
};

#if !defined(SOFANCP_CYLINDER_CONTACT_FORCE_FIELD_CPP)
extern template class SOFANCP_API
    CylinderNCPContactForceField<defaulttype::Vec3Types, defaulttype::Vec1Types>;

extern template class SOFANCP_API
    CylinderNCPContactForceField<defaulttype::Rigid3Types, defaulttype::Vec1Types>;
#endif

SOFANCP_API void registerCylinderNCPContactForceField(sofa::core::ObjectFactory* factory);

} // namespace sofa::ncp
