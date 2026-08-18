#pragma once

#include <sofa/ncp/config.h>
#include <sofa/ncp/contact/FischerBurmeisterContactForceField.h>

#include <sofa/component/odesolver/backward/StaticSolver.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>

namespace sofa::core
{
class ObjectFactory;
}

namespace sofa::ncp
{

/**
 * Static mixed mechanical/Fischer-Burmeister backend.
 *
 * Responsibilities:
 *   - project the solve-entry state once;
 *   - evaluate the complete current residual;
 *   - assemble one mixed residual Jacobian per Newton iteration;
 *   - construct exact line-search trials from the accepted checkpoint;
 *   - capture lagged compliance only after a successful nonlinear solve.
 */
class SOFANCP_API NCPStaticSolver
    : public sofa::component::odesolver::backward::StaticSolver
{
public:
    SOFA_CLASS(NCPStaticSolver,
               sofa::component::odesolver::backward::StaticSolver);

    using Inherit = sofa::component::odesolver::backward::StaticSolver;
    using ContactForceField = FischerBurmeisterContactForceField <defaulttype::Rigid3Types,defaulttype::Vec1Types>;

    NCPStaticSolver();

    void init() override;
    void solve(const core::ExecParams* params,
        SReal dt,
        core::MultiVecCoordId xResult,
        core::MultiVecDerivId vResult) override;

    SingleLink<NCPStaticSolver, ContactForceField, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_contactForceField;

    Data<SReal> d_mechanicalResidualReference;
    Data<SReal> d_complementarityResidualReference;
    Data<bool> d_debug;
    Data<bool> d_finiteDifferenceCheck;

    Data<SReal> d_lastMechanicalResidualNorm;
    Data<SReal> d_lastComplementarityResidualNorm;
    Data<bool> d_lastComplianceCaptureSucceeded;

protected:
    core::MultiVecCoordId m_solveStateId;
    core::MultiVecCoordId m_newtonStateId;
    core::MultiVecCoordId m_continuationStateId;
    core::MultiVecCoordId m_mechanicalStateId;
};

SOFANCP_API void registerNCPStaticSolver(sofa::core::ObjectFactory* factory);

} // namespace sofa::ncp
