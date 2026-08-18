#define SOFANCP_MESH_NCP_CONTACT_FORCE_FIELD_CPP

#include <sofa/ncp/contact/MeshNCPContactForceField.inl>
#include <sofa/core/ObjectFactory.h>

namespace sofa::ncp
{

void registerMeshNCPContactForceField(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(core::ObjectRegistrationData(
        "Fischer-Burmeister contact against a directly queried OBJ triangle mesh.")
        .add<MeshNCPContactForceField<defaulttype::Rigid3Types, defaulttype::Vec1Types>>()
        .add<MeshNCPContactForceField<defaulttype::Vec3Types, defaulttype::Vec1Types>>());
}

template class SOFANCP_API MeshNCPContactForceField<defaulttype::Rigid3Types, defaulttype::Vec1Types>;
template class SOFANCP_API MeshNCPContactForceField<defaulttype::Vec3Types, defaulttype::Vec1Types>;

} // namespace sofa::ncp
