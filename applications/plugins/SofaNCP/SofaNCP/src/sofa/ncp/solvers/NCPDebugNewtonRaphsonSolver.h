#pragma once

#include <sofa/ncp/config.h>
#include <sofa/component/odesolver/backward/NewtonRaphsonSolver.h>

#include <string>

namespace sofa::core
{
class ObjectFactory;
}

namespace sofa::ncp
{

struct NCPDebugResidualSummary
{
    SReal mechanicalResidualNorm = 0_sreal;
    SReal complementarityResidualNorm = 0_sreal;
    SReal scaledMechanicalSquaredNorm = 0_sreal;
    SReal scaledComplementaritySquaredNorm = 0_sreal;
    SReal minimumActiveGap = 0_sreal;
    SReal maximumPenetration = 0_sreal;
    SReal minimumLambda = 0_sreal;
    SReal maximumLambda = 0_sreal;
    sofa::Size activeContacts = 0;
    sofa::Size pinnedContacts = 0;
    sofa::Size invalidContacts = 0;
    bool valid = false;
};

struct NCPDebugDirectionSummary
{
    SReal totalNorm = 0_sreal;

    SReal translationNorm = 0_sreal;
    SReal rotationNorm = 0_sreal;
    SReal lambdaNorm = 0_sreal;

    SReal maximumNodeTranslation = 0_sreal;
    SReal maximumNodeRotation = 0_sreal;
    SReal maximumAbsLambda = 0_sreal;

    sofa::Size rigidNodes = 0;
    sofa::Size lambdaDofs = 0;

    bool valid = false;
};

/**
 * Minimal state transaction required by the NCP merit line search.
 *
 * The residual function owns two persistent checkpoints:
 *   - solve entry: restored when updateStateWhenDiverged=false;
 *   - Newton base: every trial is reconstructed from this accepted state.
 */
class SOFANCP_API NCPDebugNewtonFunctionInterface
{
public:
    virtual ~NCPDebugNewtonFunctionInterface() = default;

    virtual void storeSolveState() = 0;
    virtual void restoreSolveState() = 0;
    virtual void storeNewtonState() = 0;
    virtual void restoreNewtonState() = 0;
    virtual void setTrialFromNewtonState(SReal alpha) = 0;
    virtual NCPDebugResidualSummary currentNCPDebugSummary() const = 0;
    virtual NCPDebugDirectionSummary currentDirectionSummary() const = 0;
    virtual bool correctionIsFinite() const = 0;

    // Optional directional finite-difference diagnostic. The backend snapshots
    // the accepted Newton base after the final search direction is selected,
    // then compares each true line-search trial against the assembled J*d.
    virtual void beginFiniteDifferenceCheck() {}
    virtual void evaluateFiniteDifferenceTrial(SReal, unsigned int, unsigned int) {}

    /**
     * Apply Algorithm 4.1 search-direction safeguard to the current Newton
     * correction.
     *
     * The backend computes the weighted merit gradient
     *
     *     gradPsi = J^T W^2 R
     *
     * from the assembled residual Jacobian. If the current Newton direction
     * satisfies
     *
     *     gradPsi^T d_N <= -rho ||d_N||^p,
     *
     * it is retained. Otherwise the current correction is replaced by
     *
     *     d_G = -gradPsi.
     *
     * On success, meritSlope is gradPsi^T d for the retained direction.
     */
    virtual bool selectMeritSearchDirection(
        SReal rho,
        SReal exponent,
        SReal& meritSlope,
        bool& usedGradientFallback) = 0;
};

/**
 * Newton solver globalized by the weighted residual merit
 *
 *     Psi = 0.5 || [ R_mechanical / Mref ; Phi_FB / Cref ] ||_2^2.
 *
 * A Newton direction is retained only when it satisfies the sufficient-descent
 * test from Algorithm 4.1. Otherwise the backend replaces it by -grad(Psi).
 *
 * One matrix and one direction are built per Newton iteration. Every trial
 * evaluates the true current residual. A successful trial is retained exactly
 * as evaluated; there is no post-acceptance contact refresh.
 */
class SOFANCP_API NCPDebugNewtonRaphsonSolver
    : public sofa::component::odesolver::backward::NewtonRaphsonSolver
{
public:
    SOFA_CLASS(NCPDebugNewtonRaphsonSolver,
               sofa::component::odesolver::backward::NewtonRaphsonSolver);

    using Inherit = sofa::component::odesolver::backward::NewtonRaphsonSolver;

    NCPDebugNewtonRaphsonSolver();

    void init() override;
    void reset() override;
    void solveNCP(sofa::component::odesolver::backward::newton_raphson::BaseNonLinearFunction& function);

    bool lastSolveConverged() const { return d_lastFailureReason.getValue().empty(); }

    Data<bool> d_logIterationSummary;
    Data<bool> d_logLineSearchTrials;
    Data<unsigned int> d_nonMonotoneWindow;
    Data<unsigned int> d_lastAcceptedNewtonUpdates;
    Data<SReal> d_lastAcceptedAlpha;
    Data<std::string> d_lastFailureReason;
};

SOFANCP_API void registerNCPDebugNewtonRaphsonSolver(sofa::core::ObjectFactory* factory);

} // namespace sofa::ncp
