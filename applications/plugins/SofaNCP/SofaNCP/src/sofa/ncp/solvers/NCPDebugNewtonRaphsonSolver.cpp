#include <sofa/ncp/solvers/NCPDebugNewtonRaphsonSolver.h>

#include <sofa/component/odesolver/backward/NonLinearFunction.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/ScopedAdvancedTimer.h>
#include <sofa/helper/logging/Messaging.h>

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <string>

namespace sofa::ncp
{

using sofa::component::odesolver::backward::NewtonStatus;
using BaseNonLinearFunction = sofa::component::odesolver::backward::newton_raphson::BaseNonLinearFunction;

void registerNCPDebugNewtonRaphsonSolver(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(core::ObjectRegistrationData(
        "NCP Newton solver with Algorithm-4.1 merit descent safeguard and non-monotone Armijo backtracking.")
        .add<NCPDebugNewtonRaphsonSolver>());
}

NCPDebugNewtonRaphsonSolver::NCPDebugNewtonRaphsonSolver()
    : Inherit()
    , d_logIterationSummary(initData(&d_logIterationSummary, true,
        "logIterationSummary", "Print one base/linear/terminal line per Newton iteration."))
    , d_logLineSearchTrials(initData(&d_logLineSearchTrials, true,
        "logLineSearchTrials", "Print one line per merit-function backtracking trial."))
    , d_nonMonotoneWindow(initData(&d_nonMonotoneWindow, 1u,
        "nonMonotoneWindow", "Number of recent accepted merit values used by the non-monotone Armijo condition."))
    , d_lastAcceptedNewtonUpdates(initData(&d_lastAcceptedNewtonUpdates, 0u,
        "lastAcceptedNewtonUpdates", "Number of Newton updates retained by the last solve."))
    , d_lastAcceptedAlpha(initData(&d_lastAcceptedAlpha, 0_sreal,
        "lastAcceptedAlpha", "Step length of the most recently retained Newton update."))
    , d_lastFailureReason(initData(&d_lastFailureReason, std::string(),
        "lastFailureReason", "Failure stage and reason; empty after convergence."))
{
    static std::string groupDiagnostics{"NCP Debug"};
    d_logIterationSummary.setGroup(groupDiagnostics);
    d_logLineSearchTrials.setGroup(groupDiagnostics);
    d_nonMonotoneWindow.setGroup(groupDiagnostics);
    d_lastAcceptedNewtonUpdates.setGroup(groupDiagnostics);
    d_lastAcceptedAlpha.setGroup(groupDiagnostics);
    d_lastFailureReason.setGroup(groupDiagnostics);

    d_lastAcceptedNewtonUpdates.setReadOnly(true);
    d_lastAcceptedAlpha.setReadOnly(true);
    d_lastFailureReason.setReadOnly(true);
}

void NCPDebugNewtonRaphsonSolver::init()
{
    Inherit::init();

    const SReal coefficient = d_lineSearchCoefficient.getValue();
    if (!(coefficient > 0_sreal && coefficient < 1_sreal))
    {
        msg_error() << "lineSearchCoefficient must be strictly between 0 and 1.";
        d_componentState.setValue(core::objectmodel::ComponentState::Invalid);
    }
}

void NCPDebugNewtonRaphsonSolver::reset()
{
    Inherit::reset();
    d_lastAcceptedNewtonUpdates.setValue(0u);
    d_lastAcceptedAlpha.setValue(0_sreal);
    d_lastFailureReason.setValue(std::string());
}

namespace
{

constexpr SReal armijoCoefficient = 1e-4_sreal;   // sigma
constexpr SReal descentRho = 1e-8_sreal;           // rho
constexpr SReal descentExponent = 2_sreal;        // p

bool finiteNonNegative(const SReal value)
{
    return std::isfinite(static_cast<double>(value)) && value >= 0_sreal;
}

SReal normFromSquared(const SReal value)
{
    return std::sqrt(std::max(value, 0_sreal));
}

SReal normRatio(const SReal trialNorm, const SReal baseNorm)
{
    return trialNorm / std::max(baseNorm,1e-10_sreal);
}

bool convergedAtAbsoluteThreshold(const SReal squaredResidual, const SReal absoluteThreshold)
{
    return absoluteThreshold > 0_sreal && squaredResidual <= absoluteThreshold * absoluteThreshold;
}

bool convergedAtRelativeInitialThreshold(const SReal squaredResidual, const SReal initialSquaredResidual, const SReal relativeInitialThreshold)
{
    return relativeInitialThreshold > 0_sreal && initialSquaredResidual > std::numeric_limits<SReal>::epsilon() && squaredResidual <= relativeInitialThreshold * relativeInitialThreshold * initialSquaredResidual;
}

} // namespace

void NCPDebugNewtonRaphsonSolver::solveNCP(BaseNonLinearFunction& function)
{
    if (!isComponentStateValid())
        return;

    auto* ncpFunction = dynamic_cast<NCPDebugNewtonFunctionInterface*>(&function);
    if (!ncpFunction)
    {
        msg_error() << "NCPDebugNewtonRaphsonSolver requires NCPDebugNewtonFunctionInterface.";
        d_lastFailureReason.setValue("configuration: residual function does not implement NCP interface");
        return;
    }

    start();
    d_lastAcceptedNewtonUpdates.setValue(0u);
    d_lastAcceptedAlpha.setValue(0_sreal);
    d_lastFailureReason.setValue(std::string());

    const bool keepLastAcceptedState = d_updateStateWhenDiverged.getValue();
    const SReal absoluteThreshold = d_absoluteResidualStoppingThreshold.getValue();
    const SReal relativeInitialThreshold = d_relativeInitialStoppingThreshold.getValue();
    const unsigned int maxNewtonIterations = d_maxNbIterationsNewton.getValue();
    const unsigned int maxLineSearchTrials = std::max(d_maxNbIterationsLineSearch.getValue(), 1u);
    const SReal reductionFactor = d_lineSearchCoefficient.getValue();

    ncpFunction->storeSolveState();
    function.evaluateCurrentGuess();

    SReal squaredResidual = function.squaredNormLastEvaluation();
    NCPDebugResidualSummary summary = ncpFunction->currentNCPDebugSummary();
    const SReal initialSquaredResidual = squaredResidual;

    if (!finiteNonNegative(squaredResidual))
    {
        d_lastFailureReason.setValue("residual-evaluation: initial residual is non-finite");

        msg_error()
            << "[NCP FAIL] stage=residual-evaluation"
            << " reason=non-finite-initial-residual"
            << " R2=" << squaredResidual << ".";

        return;
    }

    if (!summary.valid)
    {
        d_lastFailureReason.setValue("residual-evaluation: invalid contact evaluation");

        msg_error()
            << "[NCP FAIL] stage=residual-evaluation"
            << " reason=invalid-contact-evaluation"
            << " contacts="
            << summary.activeContacts << "/"
            << summary.pinnedContacts << "/"
            << summary.invalidContacts
            << " M=" << summary.mechanicalResidualNorm
            << " C=" << summary.complementarityResidualNorm
            << ".";

        return;
    }

    std::deque<SReal> acceptedMeritHistory;
    acceptedMeritHistory.push_back(0.5_sreal * squaredResidual);

    unsigned int retainedUpdates = 0;

    for (unsigned int iteration = 0; iteration < maxNewtonIterations; ++iteration)
    {
        SCOPED_TIMER_VARNAME(stepTimer, "NCPMeritNewtonStep");

        const SReal baseSquaredResidual = squaredResidual;
        const SReal baseResidual = normFromSquared(baseSquaredResidual);
        const NCPDebugResidualSummary baseSummary = summary;

        const bool absoluteThresholdCriterion = convergedAtAbsoluteThreshold(baseSquaredResidual, absoluteThreshold);
        const bool relativeInitialThresholdCriterion = convergedAtRelativeInitialThreshold( baseSquaredResidual, initialSquaredResidual, relativeInitialThreshold);

        if (absoluteThresholdCriterion || relativeInitialThresholdCriterion)
        {
            if (absoluteThresholdCriterion)
            {
                static constexpr NewtonStatus convergedAbsolute("ConvergedAbsoluteResidual");
                d_status.setValue(convergedAbsolute);
            }
            else
            {
                static constexpr NewtonStatus convergedRelativeInitial("ConvergedResidualInitialRatio");
                d_status.setValue(convergedRelativeInitial);
            }

            d_lastFailureReason.setValue(std::string());

            if (f_printLog.getValue() && d_logIterationSummary.getValue())
            {
                msg_info() << "[NCP CONVERGED] iterations=" << retainedUpdates
                        << " R=" << baseResidual
                        << " M=" << baseSummary.mechanicalResidualNorm
                        << " C=" << baseSummary.complementarityResidualNorm
                        << " criterion="
                        << (absoluteThresholdCriterion ? "absolute" : "relative-initial");
            }

            sofa::helper::AdvancedTimer::valSet("nb_iterations", retainedUpdates);
            sofa::helper::AdvancedTimer::valSet("residual", baseResidual);
            return;
        }

        if (f_printLog.getValue() && d_logIterationSummary.getValue())
        {
            msg_info() << "[NCP BASE " << retainedUpdates << "]"
                    << " R=" << baseResidual
                    << " M=" << baseSummary.mechanicalResidualNorm
                    << " C=" << baseSummary.complementarityResidualNorm
                    << " contacts=" << baseSummary.activeContacts << "/"
                    << baseSummary.pinnedContacts << "/"
                    << baseSummary.invalidContacts
                    << " minGap=" << baseSummary.minimumActiveGap
                    << " penetration=" << baseSummary.maximumPenetration
                    << " lambda=[" << baseSummary.minimumLambda << ","
                    << baseSummary.maximumLambda << "]";
        }

        // Store the exact base state for this Newton iteration.
        ncpFunction->storeNewtonState();

        function.computeGradientFromCurrentGuess();
        function.solveLinearEquation();

        const SReal squaredCorrection = function.squaredNormDx();

        if (!ncpFunction->correctionIsFinite() || !finiteNonNegative(squaredCorrection))
        {
            ncpFunction->restoreNewtonState();
            function.evaluateCurrentGuess();

            squaredResidual = function.squaredNormLastEvaluation();
            summary = ncpFunction->currentNCPDebugSummary();

            if (!keepLastAcceptedState)
            {
                ncpFunction->restoreSolveState();
                function.evaluateCurrentGuess();

                squaredResidual = function.squaredNormLastEvaluation();
                summary = ncpFunction->currentNCPDebugSummary();

                retainedUpdates = 0;
                d_lastAcceptedNewtonUpdates.setValue(0u);
                d_lastAcceptedAlpha.setValue(0_sreal);
            }

            d_lastFailureReason.setValue("linear-solve: non-finite Newton correction");

            msg_error() << "[NCP FAIL] stage=linear-solve"
                        << " reason=non-finite-correction"
                        << " retainedUpdates=" << retainedUpdates << ".";

            return;
        }

        SReal meritSlope = std::numeric_limits<SReal>::quiet_NaN();
        bool usedGradientFallback = false;

        if (!ncpFunction->selectMeritSearchDirection(
                descentRho,
                descentExponent,
                meritSlope,
                usedGradientFallback)
            || !std::isfinite(static_cast<double>(meritSlope))
            || !(meritSlope < 0_sreal))
        {
            ncpFunction->restoreNewtonState();
            function.evaluateCurrentGuess();

            squaredResidual = function.squaredNormLastEvaluation();
            summary = ncpFunction->currentNCPDebugSummary();

            if (!keepLastAcceptedState)
            {
                ncpFunction->restoreSolveState();
                function.evaluateCurrentGuess();

                squaredResidual = function.squaredNormLastEvaluation();
                summary = ncpFunction->currentNCPDebugSummary();

                retainedUpdates = 0;
                d_lastAcceptedNewtonUpdates.setValue(0u);
                d_lastAcceptedAlpha.setValue(0_sreal);
            }

            d_lastFailureReason.setValue(
                "direction-selection: neither Newton nor negative merit gradient produced a valid descent direction");

            msg_warning() << "[NCP FAIL] stage=direction-selection"
                        << " reason=no-valid-descent-direction"
                        << " retainedUpdates=" << retainedUpdates << ".";

            return;
        }
        
        const auto direction = ncpFunction->currentDirectionSummary();

        if (f_printLog.getValue() && d_logIterationSummary.getValue())
        {
            msg_info() << "[NCP DIR " << retainedUpdates << "]"
                    << " type=" << (usedGradientFallback ? "gradient" : "newton")
                    << " |dz|=" << normFromSquared(function.squaredNormDx())
                    << " dPsi=" << meritSlope
                    << " |dx|=" << direction.translationNorm
                    << " |dtheta|=" << direction.rotationNorm
                    << " |dlambda|=" << direction.lambdaNorm
                    << " maxDx=" << direction.maximumNodeTranslation
                    << " maxDtheta=" << direction.maximumNodeRotation
                    << " maxDlambda=" << direction.maximumAbsLambda;
        }

        ncpFunction->beginFiniteDifferenceCheck();

        const SReal baseMerit = 0.5_sreal * baseSquaredResidual;
        SReal referenceMerit = baseMerit;

        // for (const SReal merit : acceptedMeritHistory)
        //     referenceMerit = std::max(referenceMerit, merit);

        bool accepted = false;
        bool sawValidTrial = false;

        SReal acceptedAlpha = 0_sreal;
        SReal alpha = 1_sreal;

        // Keep diagnostics of the genuinely best evaluated trial.
        SReal bestAlpha = 0_sreal;
        SReal bestRatio = std::numeric_limits<SReal>::infinity();
        SReal bestQM = std::numeric_limits<SReal>::infinity();
        SReal bestQC = std::numeric_limits<SReal>::infinity();
        SReal bestTrialResidual = std::numeric_limits<SReal>::infinity();
        SReal bestFdMeritSlope = std::numeric_limits<SReal>::quiet_NaN();
        SReal bestFdSlopeRatio = std::numeric_limits<SReal>::quiet_NaN();
        bool bestValidTrial = false;
        bool bestArmijoAccepted = false;

        for (unsigned int trial = 0; trial < maxLineSearchTrials; ++trial)
        {
            ncpFunction->setTrialFromNewtonState(alpha);
            function.evaluateCurrentGuess();

            const SReal trialSquaredResidual = function.squaredNormLastEvaluation();
            const NCPDebugResidualSummary trialSummary = ncpFunction->currentNCPDebugSummary();

            const bool validTrial = finiteNonNegative(trialSquaredResidual) && trialSummary.valid;
            const SReal trialResidual = validTrial ? normFromSquared(trialSquaredResidual) : std::numeric_limits<SReal>::infinity();

            const SReal qR = normRatio(trialResidual, baseResidual);
            const SReal qM = normRatio(trialSummary.mechanicalResidualNorm, baseSummary.mechanicalResidualNorm);
            const SReal qC = normRatio(trialSummary.complementarityResidualNorm, baseSummary.complementarityResidualNorm);

            ncpFunction->evaluateFiniteDifferenceTrial(alpha, retainedUpdates, trial);

            const SReal fdMeritSlope = validTrial && alpha > 0_sreal
                ? (trialSquaredResidual - baseSquaredResidual) / (2_sreal * alpha)
                : std::numeric_limits<SReal>::quiet_NaN();

            const SReal fdSlopeRatio =
                std::abs(meritSlope) > std::numeric_limits<SReal>::epsilon()
                    ? fdMeritSlope / meritSlope
                    : std::numeric_limits<SReal>::quiet_NaN();

            // GLL non-monotone Armijo condition on Psi = 0.5 ||R||^2.
            // The trial is compared against the largest accepted merit in the
            // recent window, not necessarily against the current base merit.
            const SReal requiredMerit = referenceMerit + armijoCoefficient * alpha * meritSlope;
            const bool armijoAccepted = validTrial && 0.5_sreal * trialSquaredResidual < 1.0 * referenceMerit;
            const bool trialAccepted = armijoAccepted;

            if (validTrial)
            {
                sawValidTrial = true;

                if (qR < bestRatio)
                {
                    bestValidTrial = true;
                    bestAlpha = alpha;
                    bestRatio = qR;
                    bestQM = qM;
                    bestQC = qC;
                    bestTrialResidual = trialResidual;
                    bestFdMeritSlope = fdMeritSlope;
                    bestFdSlopeRatio = fdSlopeRatio;
                    bestArmijoAccepted = armijoAccepted;
                }
            }

            if (f_printLog.getValue() && d_logLineSearchTrials.getValue() && validTrial)
            {
                msg_info() << "[NCP LS " << retainedUpdates << "." << trial << "]"
                        << " direction=" << (usedGradientFallback ? "gradient" : "newton")
                        << " alpha=" << alpha
                        << " qR=" << qR
                        << " qM=" << qM
                        << " qC=" << qC
                        << " fdSlope=" << fdMeritSlope
                        << " predictedSlope=" << meritSlope
                        << " fdSlopeRatio=" << fdSlopeRatio
                        << " referenceR=" << std::sqrt(2_sreal * referenceMerit)
                        << " armijo=" << armijoAccepted
                        << " stepMaxDx=" << alpha * direction.maximumNodeTranslation
                        << " stepMaxDtheta=" << alpha * direction.maximumNodeRotation
                        << " stepMaxDlambda=" << alpha * direction.maximumAbsLambda;
            }

            if (trialAccepted)
            {
                accepted = true;
                acceptedAlpha = alpha;

                squaredResidual = trialSquaredResidual;
                summary = trialSummary;

                if (f_printLog.getValue() && d_logIterationSummary.getValue())
                {
                    msg_info() << "[NCP ACCEPT " << retainedUpdates << "]"
                            << " type=nonmonotone-armijo"
                            << " alpha=" << alpha
                            << " qR=" << qR
                            << " qM=" << qM
                            << " qC=" << qC;
                }

                break;
            }

            alpha *= reductionFactor;
        }

        if (accepted)
        {
            ++retainedUpdates;

            acceptedMeritHistory.push_back(0.5_sreal * squaredResidual);

            const unsigned int nonMonotoneWindow = std::max(d_nonMonotoneWindow.getValue(), 1u);
            while (acceptedMeritHistory.size() > nonMonotoneWindow)
                acceptedMeritHistory.pop_front();

            d_lastAcceptedNewtonUpdates.setValue(retainedUpdates);
            d_lastAcceptedAlpha.setValue(acceptedAlpha);

            // Current coordinates and residual already correspond exactly to the
            // accepted trial. Do not reevaluate or refresh here.
            continue;
        }

        // ---------------------------------------------------------------------
        // No trial accepted.
        //
        // Do NOT reapply bestAlpha for diagnostics. All diagnostics needed for
        // the best trial were stored while that trial was actually evaluated.
        // First restore the exact Newton base.
        // ---------------------------------------------------------------------

        ncpFunction->restoreNewtonState();
        function.evaluateCurrentGuess();

        squaredResidual = function.squaredNormLastEvaluation();
        summary = ncpFunction->currentNCPDebugSummary();

        // If requested, discard every Newton update accepted during this solve
        // and restore the solve-entry state.
        if (!keepLastAcceptedState)
        {
            ncpFunction->restoreSolveState();
            function.evaluateCurrentGuess();

            squaredResidual = function.squaredNormLastEvaluation();
            summary = ncpFunction->currentNCPDebugSummary();

            retainedUpdates = 0;
            d_lastAcceptedNewtonUpdates.setValue(0u);
            d_lastAcceptedAlpha.setValue(0_sreal);
        }

        static constexpr NewtonStatus divergedLineSearch("DivergedLineSearch");
        d_status.setValue(divergedLineSearch);

        d_lastFailureReason.setValue(
            sawValidTrial
                ? "line-search: no trial satisfied non-monotone Armijo"
                : "line-search: every trial had invalid residual/contact geometry");

        msg_warning_when(d_warnWhenLineSearchFails.getValue())
            << "[NCP FAIL] stage=line-search reason="
            << (sawValidTrial ? "no-acceptable-decrease" : "all-trials-invalid")
            << " bestAlpha=" << bestAlpha
            << " bestqR=" << bestRatio
            << " bestqM=" << bestQM
            << " bestqC=" << bestQC
            << " bestResidual=" << bestTrialResidual
            << " bestValid=" << bestValidTrial
            << " bestArmijo=" << bestArmijoAccepted
            << " bestFdSlope=" << bestFdMeritSlope
            << " predictedSlope=" << meritSlope
            << " bestFdSlopeRatio=" << bestFdSlopeRatio
            << " retainedUpdates=" << retainedUpdates << ".";

        sofa::helper::AdvancedTimer::valSet("nb_iterations", retainedUpdates);
        sofa::helper::AdvancedTimer::valSet("residual", normFromSquared(squaredResidual));
        return;
    }

    if (!keepLastAcceptedState)
    {
        ncpFunction->restoreSolveState();
        function.evaluateCurrentGuess();
        squaredResidual = function.squaredNormLastEvaluation();
        retainedUpdates = 0;
        d_lastAcceptedNewtonUpdates.setValue(0u);
        d_lastAcceptedAlpha.setValue(0_sreal);
    }

    static constexpr NewtonStatus divergedMaxIterations("DivergedMaxIterations");
    d_status.setValue(divergedMaxIterations);
    d_lastFailureReason.setValue("max-iterations: nonlinear solve did not converge");

    msg_warning_when(d_warnWhenDiverge.getValue())
        << "[NCP FAIL] stage=max-iterations retainedUpdates=" << retainedUpdates
        << " R=" << normFromSquared(squaredResidual) << ".";

    sofa::helper::AdvancedTimer::valSet("nb_iterations", retainedUpdates);
    sofa::helper::AdvancedTimer::valSet("residual", normFromSquared(squaredResidual));
}

} // namespace sofa::ncp
