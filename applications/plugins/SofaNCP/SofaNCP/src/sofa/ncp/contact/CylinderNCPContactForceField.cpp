/****************************************************************************
* Open-cylinder Fischer-Burmeister contact forcefield registration.
****************************************************************************/
#define SOFANCP_CYLINDER_CONTACT_FORCE_FIELD_CPP

#include <sofa/ncp/contact/CylinderNCPContactForceField.inl>
#include <sofa/core/behavior/MixedInteractionForceField.inl>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa::ncp
{

using namespace sofa::defaulttype;

void registerCylinderNCPContactForceField(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(core::ObjectRegistrationData(
        "Open-cylinder Fischer-Burmeister NCP contact forcefield")
        .add< CylinderNCPContactForceField<Vec3Types, Vec1Types> >()
        .add< CylinderNCPContactForceField<Rigid3Types, Vec1Types> >());
}

template class SOFANCP_API CylinderNCPContactForceField<Vec3Types, Vec1Types>;
template class SOFANCP_API CylinderNCPContactForceField<Rigid3Types, Vec1Types>;

} // namespace sofa::ncp
