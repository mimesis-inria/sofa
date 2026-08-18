/****************************************************************************
* Infinite-plane Fischer-Burmeister NCP contact forcefield registration.
****************************************************************************/
#define SOFANCP_PLANE_CONTACT_FORCE_FIELD_CPP

#include <sofa/ncp/contact/PlaneNCPContactForceField.inl>
#include <sofa/core/behavior/MixedInteractionForceField.inl>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa::ncp
{

using namespace sofa::defaulttype;

void registerPlaneNCPContactForceField(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(core::ObjectRegistrationData(
        "Infinite-plane Fischer-Burmeister NCP contact forcefield. "
        "The configured normal points into the admissible half-space.")
        .add< PlaneNCPContactForceField<Vec3Types, Vec1Types> >()
        .add< PlaneNCPContactForceField<Rigid3Types, Vec1Types> >());
}

template class SOFANCP_API
    PlaneNCPContactForceField<Vec3Types, Vec1Types>;

template class SOFANCP_API
    PlaneNCPContactForceField<Rigid3Types, Vec1Types>;

} // namespace sofa::ncp
