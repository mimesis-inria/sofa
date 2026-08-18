#pragma once

#include <sofa/ncp/config.h>

#include <sofa/core/behavior/BaseMatrixLinearSystem.h>
#include <sofa/core/behavior/LinearSolverAccessor.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/type/Mat.h>
#include <sofa/type/Vec.h>
#include <sofa/type/vector.h>

namespace sofa::core { class ObjectFactory; }

namespace sofa::ncp
{

template<class TDataTypes>
class StaticComplianceProvider final
    : public core::behavior::LinearSolverAccessor
{
public:
    SOFA_CLASS(
        SOFA_TEMPLATE(StaticComplianceProvider, TDataTypes),
        core::behavior::LinearSolverAccessor);

    using Inherit = core::behavior::LinearSolverAccessor;
    using DataTypes = TDataTypes;
    using MechanicalState = core::behavior::MechanicalState<DataTypes>;
    using Deriv = typename DataTypes::Deriv;
    using Real = typename DataTypes::Real;
    using Vec3 = sofa::type::Vec<3, Real>;
    using Mat3 = sofa::type::Mat<3, 3, Real>;

    static constexpr sofa::Size DerivDim = Deriv::total_size;
    static constexpr sofa::Size TranslationalDim = 3;
    static_assert(DerivDim >= TranslationalDim,
        "StaticComplianceProvider requires three translational DOFs.");

    SingleLink<StaticComplianceProvider<DataTypes>,
        MechanicalState,
        BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_mechanicalState;

    SingleLink<StaticComplianceProvider<DataTypes>,
        core::behavior::BaseMatrixLinearSystem,
        BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_linearSystem;

    Data<bool> d_debug;

protected:
    StaticComplianceProvider();

public:
    void init() override;
    void reinit() override;

    /** Transactional capture: a failed capture preserves the last valid data. */
    bool captureFromCurrentLinearSystem();
    bool hasCompliance(sofa::Size expectedPointCount) const;
    void invalidate();
    Real scalarCompliance(sofa::Index pointIndex, const Vec3& gapGradient) const;

private:
    bool m_valid = false;
    sofa::Size m_generation = 0;
    sofa::type::vector<Mat3> m_pointCompliance;

    core::behavior::BaseMatrixLinearSystem* resolveLinearSystem() const;
    bool extractPointCompliance(const sofa::linearalgebra::BaseMatrix& matrix);

    sofa::Size localDof(sofa::Index pointIndex, sofa::Size component) const
    {
        return sofa::Size(pointIndex) * DerivDim + component;
    }
};

#if !defined(SOFANCP_STATIC_COMPLIANCE_PROVIDER_CPP)
extern template class SOFANCP_API StaticComplianceProvider<defaulttype::Rigid3Types>;
extern template class SOFANCP_API StaticComplianceProvider<defaulttype::Vec3Types>;
#endif

SOFANCP_API void registerStaticComplianceProvider(sofa::core::ObjectFactory* factory);

} // namespace sofa::ncp
