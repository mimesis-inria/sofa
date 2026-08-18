/****************************************************************************
* SDF Fischer-Burmeister NCP contact forcefield registration.
****************************************************************************/
#define SOFANCP_SDF_CONTACT_FORCE_FIELD_CPP

#include <sofa/ncp/contact/SignedDistanceFieldNCPContactForceField.inl>
#include <sofa/core/behavior/MixedInteractionForceField.inl>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa::ncp
{

using namespace sofa::defaulttype;

void registerSignedDistanceFieldNCPContactForceField(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(core::ObjectRegistrationData(
        "Signed-distance-field Fischer-Burmeister NCP contact forcefield")
        .add< SignedDistanceFieldNCPContactForceField<Vec3Types, Vec1Types> >()
        .add< SignedDistanceFieldNCPContactForceField<Rigid3Types, Vec1Types> >());
}

template class SOFANCP_API SignedDistanceFieldNCPContactForceField<Vec3Types, Vec1Types>;
template class SOFANCP_API SignedDistanceFieldNCPContactForceField<Rigid3Types, Vec1Types>;

} // namespace sofa::ncp
