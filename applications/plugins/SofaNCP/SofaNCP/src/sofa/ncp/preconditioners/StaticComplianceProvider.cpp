#define SOFANCP_STATIC_COMPLIANCE_PROVIDER_CPP

#include <sofa/ncp/preconditioners/StaticComplianceProvider.inl>
#include <sofa/core/ObjectFactory.h>

namespace sofa::ncp
{

void registerStaticComplianceProvider(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(core::ObjectRegistrationData(
        "Transactionally captures lagged per-point translational compliance.")
        .add<StaticComplianceProvider<defaulttype::Rigid3Types>>()
        .add<StaticComplianceProvider<defaulttype::Vec3Types>>());
}

template class SOFANCP_API StaticComplianceProvider<defaulttype::Rigid3Types>;
template class SOFANCP_API StaticComplianceProvider<defaulttype::Vec3Types>;

} // namespace sofa::ncp
