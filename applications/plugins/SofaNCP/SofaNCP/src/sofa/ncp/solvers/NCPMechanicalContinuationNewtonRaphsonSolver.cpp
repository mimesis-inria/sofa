#include <sofa/ncp/solvers/NCPMechanicalContinuationNewtonRaphsonSolver.h>
#include <sofa/ncp/solvers/NCPMechanicalContinuationFunctionInterface.h>

#include <sofa/component/odesolver/backward/NonLinearFunction.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/ScopedAdvancedTimer.h>
#include <sofa/helper/logging/Messaging.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace sofa::ncp
{

using sofa::component::odesolver::backward::NewtonStatus;
using BaseNonLinearFunction = sofa::component::odesolver::backward::newton_raphson::BaseNonLinearFunction;

namespace
{

constexpr SReal armijoCoefficient = 1e-4_sreal;
constexpr SReal descentRho = 1e-8_sreal;
constexpr SReal descentExponent = 2.1_sreal;

bool finiteNonNegative(SReal value)
{
    return std::isfinite(static_cast<double>(value)) && value >= 0_sreal;
}

SReal normFromSquared(SReal value)
{
    return std::sqrt(std::max(value, 0_sreal));
}

SReal ratio(SReal value, SReal reference)
{
    return reference > std::numeric_limits<SReal>::epsilon()
        ? value / reference
        : (value <= std::numeric_limits<SReal>::epsilon() ? 1_sreal : std::numeric_limits<SReal>::infinity());
}

bool convergedAt(SReal current2, SReal initial2, SReal absoluteTolerance, SReal relativeTolerance)
{
    if (current2 <= std::numeric_limits<SReal>::epsilon())
        return true;

    if (absoluteTolerance > 0_sreal && current2 <= absoluteTolerance * absoluteTolerance)
        return true;

    return relativeTolerance > 0_sreal
        && initial2 > std::numeric_limits<SReal>::epsilon()
        && current2 <= relativeTolerance * relativeTolerance * initial2;
}

bool fullNCPConverged(const NCPDebugResidualSummary& current, const NCPDebugResidualSummary& initial,
    SReal absoluteTolerance, SReal relativeTolerance)
{
    return convergedAt(current.scaledMechanicalSquaredNorm, initial.scaledMechanicalSquaredNorm,
               absoluteTolerance, relativeTolerance)
        && convergedAt(current.scaledComplementaritySquaredNorm, initial.scaledComplementaritySquaredNorm,
               absoluteTolerance, relativeTolerance);
}

struct NewtonIterationScope
{
    explicit NewtonIterationScope(BaseNonLinearFunction& function) : m_function(function)
    {
        m_function.startNewtonIteration();
    }

    ~NewtonIterationScope()
    {
        m_function.endNewtonIteration();
    }

    BaseNonLinearFunction& m_function;
};

} // namespace

void registerNCPMechanicalContinuationNewtonRaphsonSolver(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(core::ObjectRegistrationData(
        "Mixed Fischer-Burmeister Newton solver with target-tracked fixed-lambda continuation.")
        .add<NCPMechanicalContinuationNewtonRaphsonSolver>());
}

NCPMechanicalContinuationNewtonRaphsonSolver::NCPMechanicalContinuationNewtonRaphsonSolver()
    : Inherit()
    , d_lambdaStepInitial(initData(&d_lambdaStepInitial, 1_sreal, "lambdaStepInitial",
        "Initial infinity-norm substep used to inject one outer Newton dLambda target."))
    , d_lambdaStepMin(initData(&d_lambdaStepMin, 1e-4_sreal, "lambdaStepMin",
        "Smallest lambda continuation substep before target failure."))
    , d_lambdaStepMax(initData(&d_lambdaStepMax, 1e3_sreal, "lambdaStepMax",
        "Largest adaptive lambda continuation substep."))
    , d_lambdaStepGrowth(initData(&d_lambdaStepGrowth, 1.5_sreal, "lambdaStepGrowth",
        "Substep growth after an easy fixed-lambda mechanical solve."))
    , d_lambdaStepShrink(initData(&d_lambdaStepShrink, 0.5_sreal, "lambdaStepShrink",
        "Substep shrink after a difficult or failed mechanical solve."))
    , d_maxMechanicalIterations(initData(&d_maxMechanicalIterations, 12u, "maxMechanicalIterations",
        "Maximum mechanical Newton iterations per lambda substep."))
    , d_maxMechanicalLineSearchTrials(initData(&d_maxMechanicalLineSearchTrials, 12u, "maxMechanicalLineSearchTrials",
        "Maximum residual-decrease backtracks per inner mechanical Newton iteration."))
    , d_maxLambdaStageRetries(initData(&d_maxLambdaStageRetries, 8u, "maxLambdaStageRetries",
        "Maximum substep reductions attempted from one accepted partial target state."))
    , d_maxContinuationStages(initData(&d_maxContinuationStages, 200u, "maxContinuationStages",
        "Maximum accepted lambda substeps inside one static solve."))
    , d_lastLambdaStage(initData(&d_lastLambdaStage, 0_sreal, "lastLambdaStage",
        "Infinity norm of the latest accepted lambda substep."))
    , d_lastMechanicalIterations(initData(&d_lastMechanicalIterations, 0u, "lastMechanicalIterations",
        "Mechanical Newton iterations used by the latest accepted lambda substep."))
    , d_lastContinuationStages(initData(&d_lastContinuationStages, 0u, "lastContinuationStages",
        "Accepted continuation substeps in the latest static solve."))
{
    d_lastLambdaStage.setReadOnly(true);
    d_lastMechanicalIterations.setReadOnly(true);
    d_lastContinuationStages.setReadOnly(true);
}

void NCPMechanicalContinuationNewtonRaphsonSolver::reset()
{
    Inherit::reset();
    d_lastLambdaStage.setValue(0_sreal);
    d_lastMechanicalIterations.setValue(0u);
    d_lastContinuationStages.setValue(0u);
}

void NCPMechanicalContinuationNewtonRaphsonSolver::solveNCPContinuation(BaseNonLinearFunction& function)
{
    if (!isComponentStateValid())
        return;

    auto* ncp = dynamic_cast<NCPDebugNewtonFunctionInterface*>(&function);
    auto* continuation = dynamic_cast<NCPMechanicalContinuationFunctionInterface*>(&function);

    if (!ncp || !continuation)
    {
        msg_error() << "NCPMechanicalContinuationNewtonRaphsonSolver requires NCP debug and continuation interfaces.";
        return;
    }

    start();

    d_lastAcceptedNewtonUpdates.setValue(0u);
    d_lastAcceptedAlpha.setValue(0_sreal);
    d_lastFailureReason.setValue(std::string());
    d_lastLambdaStage.setValue(0_sreal);
    d_lastMechanicalIterations.setValue(0u);
    d_lastContinuationStages.setValue(0u);

    const bool keepLastAcceptedState = d_updateStateWhenDiverged.getValue();
    const SReal absoluteTolerance = d_absoluteResidualStoppingThreshold.getValue();
    const SReal relativeTolerance = d_relativeInitialStoppingThreshold.getValue();
    const SReal reductionFactor = d_lineSearchCoefficient.getValue();

    const SReal lambdaStepMin = std::max(d_lambdaStepMin.getValue(), std::numeric_limits<SReal>::epsilon());
    const SReal lambdaStepMax = std::max(d_lambdaStepMax.getValue(), lambdaStepMin);
    SReal lambdaStep = std::clamp(d_lambdaStepInitial.getValue(), lambdaStepMin, lambdaStepMax);

    const unsigned int maxOuterIterations = d_maxNbIterationsNewton.getValue();
    const unsigned int maxMixedTrials = std::max(d_maxNbIterationsLineSearch.getValue(), 1u);
    const unsigned int maxMechanicalIterations = std::max(d_maxMechanicalIterations.getValue(), 1u);
    const unsigned int maxMechanicalTrials = std::max(d_maxMechanicalLineSearchTrials.getValue(), 1u);
    const unsigned int maxStageRetries = std::max(d_maxLambdaStageRetries.getValue(), 1u);
    const unsigned int maxContinuationStages = std::max(d_maxContinuationStages.getValue(), 1u);

    ncp->storeSolveState();
    function.evaluateCurrentGuess();

    SReal residual2 = function.squaredNormLastEvaluation();
    NCPDebugResidualSummary summary = ncp->currentNCPDebugSummary();
    const NCPDebugResidualSummary initialSummary = summary;

    if (!finiteNonNegative(residual2) || !summary.valid)
    {
        static constexpr NewtonStatus invalidResidual("DivergedInvalidResidual");
        d_status.setValue(invalidResidual);
        d_lastFailureReason.setValue("invalid initial residual/contact state");
        msg_error() << "[NCP FAIL] invalid initial residual/contact state.";
        return;
    }

    unsigned int acceptedOuterUpdates = 0;
    unsigned int continuationStages = 0;
    bool continuationMode = false;
    bool failed = false;

    auto reportConverged = [&](const NCPDebugResidualSummary& stateSummary, SReal stateResidual2)
    {
        static constexpr NewtonStatus converged("ConvergedEquilibrium");
        d_status.setValue(converged);
        d_lastFailureReason.setValue(std::string());

        if (f_printLog.getValue() && d_logIterationSummary.getValue())
            msg_info() << "[NCP CONVERGED]"
                       << " outer=" << acceptedOuterUpdates
                       << " cont=" << continuationStages
                       << " R=" << normFromSquared(stateResidual2)
                       << " M=" << stateSummary.mechanicalResidualNorm
                       << " C=" << stateSummary.complementarityResidualNorm;
    };

    auto solveMechanicalEquilibrium = [&](unsigned int stageIndex, unsigned int& iterations)
    {
        const SReal mechanicalInitial2 = continuation->scaledMechanicalSquaredNorm();

        auto mechanicalConverged = [&](SReal value)
        {
            return continuation->continuationStateValid()
                && finiteNonNegative(value)
                && convergedAt(value, mechanicalInitial2, absoluteTolerance, relativeTolerance);
        };

        iterations = 0;

        while (!mechanicalConverged(continuation->scaledMechanicalSquaredNorm())
            && iterations < maxMechanicalIterations)
        {
            const SReal mechanicalBase2 = continuation->scaledMechanicalSquaredNorm();
            continuation->storeMechanicalState();

            if (!continuation->solveMechanicalLinearEquation())
                return false;

            bool accepted = false;
            SReal alpha = 1_sreal;

            for (unsigned int trial = 0; trial < maxMechanicalTrials; ++trial)
            {
                continuation->setMechanicalTrialFromState(alpha);
                continuation->evaluateMechanicalTrial();

                const SReal trial2 = continuation->scaledMechanicalSquaredNorm();
                const bool valid = continuation->continuationStateValid() && finiteNonNegative(trial2);
                const bool decrease = valid && trial2 < mechanicalBase2;

                if (f_printLog.getValue() && d_logLineSearchTrials.getValue())
                    msg_info() << "[NCP MECH " << stageIndex << "." << iterations << "." << trial << "]"
                               << " alpha=" << alpha
                               << " qM=" << (valid ? std::sqrt(ratio(trial2, mechanicalBase2))
                                                  : std::numeric_limits<SReal>::infinity())
                               << " accepted=" << decrease;

                if (decrease)
                {
                    accepted = true;
                    break;
                }

                alpha *= reductionFactor;
            }

            if (!accepted)
            {
                continuation->restoreMechanicalState();
                return false;
            }

            ++iterations;
        }

        return mechanicalConverged(continuation->scaledMechanicalSquaredNorm());
    };

    for (unsigned int outer = 0; outer < maxOuterIterations && !failed; ++outer)
    {
        SCOPED_TIMER("NCPTargetTrackedOuter");
        NewtonIterationScope iterationScope(function);

        function.evaluateCurrentGuess();
        residual2 = function.squaredNormLastEvaluation();
        summary = ncp->currentNCPDebugSummary();

        if (!finiteNonNegative(residual2) || !summary.valid)
        {
            d_lastFailureReason.setValue("invalid outer residual/contact state");
            failed = true;
            break;
        }

        if (fullNCPConverged(summary, initialSummary, absoluteTolerance, relativeTolerance))
        {
            reportConverged(summary, residual2);
            return;
        }

        const SReal baseResidual2 = residual2;
        const NCPDebugResidualSummary baseSummary = summary;

        // One outer linearization. m_newtonState remains fixed until this entire
        // dLambda target has either been applied to 100% or the outer step fails.
        ncp->storeNewtonState();
        function.computeGradientFromCurrentGuess();
        function.solveLinearEquation();

        if (!ncp->correctionIsFinite())
        {
            d_lastFailureReason.setValue("non-finite mixed Newton correction");
            failed = true;
            break;
        }

        const SReal rawLambda = continuation->maxAbsLambdaCorrection();
        if (!std::isfinite(static_cast<double>(rawLambda)) || !continuation->cacheLambdaPredictor())
        {
            d_lastFailureReason.setValue("invalid lambda predictor");
            failed = true;
            break;
        }

        if (f_printLog.getValue() && d_logIterationSummary.getValue())
            msg_info() << "[NCP BASE " << outer << "]"
                       << " R=" << normFromSquared(baseResidual2)
                       << " M=" << baseSummary.mechanicalResidualNorm
                       << " C=" << baseSummary.complementarityResidualNorm
                       << " raw|dlambda|=" << rawLambda
                       << " step=" << lambdaStep
                       << " mode=" << (continuationMode ? "continuation" : "mixed");

        /**
         * Before continuation is activated, preserve the benchmark mixed solver.
         * A large dLambda activates continuation immediately. A failed mixed
         * globalization also activates continuation using the SAME cached Newton
         * dLambda predictor.
         */
        if (!continuationMode && rawLambda <= lambdaStep)
        {
            SReal meritSlope = std::numeric_limits<SReal>::quiet_NaN();
            bool usedGradientFallback = false;
            bool mixedAccepted = false;

            if (ncp->selectMeritSearchDirection(descentRho, descentExponent, meritSlope, usedGradientFallback))
            {
                SReal alpha = 1_sreal;

                for (unsigned int trial = 0; trial < maxMixedTrials; ++trial)
                {
                    ncp->setTrialFromNewtonState(alpha);
                    function.evaluateCurrentGuess();

                    const SReal trial2 = function.squaredNormLastEvaluation();
                    const NCPDebugResidualSummary trialSummary = ncp->currentNCPDebugSummary();
                    const bool valid = finiteNonNegative(trial2) && trialSummary.valid;
                    const SReal qR = valid ? normFromSquared(trial2) / normFromSquared(baseResidual2)
                                           : std::numeric_limits<SReal>::infinity();
                    const SReal armijoBound = baseResidual2 + 2_sreal * armijoCoefficient * alpha * meritSlope;
                    const bool accepted = valid && trial2 <= armijoBound;

                    if (f_printLog.getValue() && d_logLineSearchTrials.getValue())
                        msg_info() << "[NCP LS " << outer << "." << trial << "]"
                                   << " alpha=" << alpha
                                   << " qR=" << qR
                                   << " qM=" << ratio(trialSummary.mechanicalResidualNorm, baseSummary.mechanicalResidualNorm)
                                   << " qC=" << ratio(trialSummary.complementarityResidualNorm, baseSummary.complementarityResidualNorm)
                                   << " accepted=" << accepted;

                    if (accepted)
                    {
                        mixedAccepted = true;
                        residual2 = trial2;
                        summary = trialSummary;
                        d_lastAcceptedAlpha.setValue(alpha);
                        break;
                    }

                    alpha *= reductionFactor;
                }
            }

            if (mixedAccepted)
            {
                ++acceptedOuterUpdates;
                d_lastAcceptedNewtonUpdates.setValue(acceptedOuterUpdates);

                if (f_printLog.getValue() && d_logIterationSummary.getValue())
                    msg_info() << "[NCP MIXED " << acceptedOuterUpdates << "]"
                               << " R=" << normFromSquared(residual2)
                               << " M=" << summary.mechanicalResidualNorm
                               << " C=" << summary.complementarityResidualNorm
                               << " direction=" << (usedGradientFallback ? "gradient" : "newton");

                continue;
            }

            // The monolithic direction was not globally usable. Roll back the
            // state, but retain the cached Newton dLambda and globalize that target.
            ncp->restoreNewtonState();
            function.evaluateCurrentGuess();
            continuationMode = true;

            if (f_printLog.getValue() && d_logIterationSummary.getValue())
                msg_info() << "[NCP MODE] mixed -> continuation"
                           << " raw|dlambda|=" << rawLambda;
        }
        else if (!continuationMode)
        {
            continuationMode = true;

            if (f_printLog.getValue() && d_logIterationSummary.getValue())
                msg_info() << "[NCP MODE] mixed -> continuation"
                           << " raw|dlambda|=" << rawLambda
                           << " trigger=" << lambdaStep;
        }

        if (!(rawLambda > std::numeric_limits<SReal>::epsilon()))
        {
            d_lastFailureReason.setValue("continuation requested with zero lambda predictor");
            ncp->restoreNewtonState();
            failed = true;
            break;
        }

        /**
         * Target-tracked continuation of ONE outer Newton dLambda_k.
         *
         * eta=0 : lambda = lambda_k
         * eta=1 : lambda = lambda_k + dLambda_k
         *
         * The cached dLambda_k is never recomputed inside this loop. Each accepted
         * substep advances eta, then solves Rx=0 with lambda fixed. Only after eta
         * reaches one do we evaluate the full NCP and allow a new outer linearization.
         */
        SReal eta = 0_sreal;
        unsigned int targetStages = 0;
        continuation->storeContinuationState();

        while (eta < 1_sreal && !failed)
        {
            if (continuationStages >= maxContinuationStages)
            {
                d_lastFailureReason.setValue("maximum continuation stages reached");
                failed = true;
                break;
            }

            const SReal remainingFraction = 1_sreal - eta;
            const SReal remainingMagnitude = rawLambda * remainingFraction;
            bool stageAccepted = false;

            for (unsigned int retry = 0; retry < maxStageRetries && !stageAccepted; ++retry)
            {
                continuation->restoreContinuationState();

                const SReal requestedMagnitude = std::min(lambdaStep, remainingMagnitude);
                const SReal fraction = requestedMagnitude / rawLambda;
                SReal appliedMagnitude = 0_sreal;

                if (!continuation->applyLambdaPredictorFraction(fraction, appliedMagnitude))
                {
                    d_lastFailureReason.setValue("failed to apply lambda predictor fraction");
                    failed = true;
                    break;
                }

                unsigned int mechanicalIterations = 0;

                if (solveMechanicalEquilibrium(continuationStages, mechanicalIterations))
                {
                    eta = std::min(1_sreal, eta + fraction);
                    stageAccepted = true;
                    ++targetStages;
                    ++continuationStages;

                    // Commit this partial target only for the next substep/retry.
                    // m_newtonState remains the original outer base.
                    continuation->storeContinuationState();

                    d_lastLambdaStage.setValue(appliedMagnitude);
                    d_lastMechanicalIterations.setValue(mechanicalIterations);
                    d_lastContinuationStages.setValue(continuationStages);

                    if (mechanicalIterations <= 3)
                        lambdaStep = std::min(lambdaStep * d_lambdaStepGrowth.getValue(), lambdaStepMax);
                    else if (mechanicalIterations > 10)
                        lambdaStep = std::max(lambdaStep * d_lambdaStepShrink.getValue(), lambdaStepMin);

                    if (f_printLog.getValue() && d_logIterationSummary.getValue())
                        msg_info() << "[NCP CONT " << continuationStages << "]"
                                   << " target=" << rawLambda
                                   << " eta=" << eta
                                   << " applied=" << appliedMagnitude
                                   << " remaining=" << rawLambda * (1_sreal - eta)
                                   << " mechIt=" << mechanicalIterations
                                   << " nextStep=" << lambdaStep;

                    break;
                }

                continuation->restoreContinuationState();

                const SReal previousStep = lambdaStep;
                lambdaStep = std::max(lambdaStep * d_lambdaStepShrink.getValue(), lambdaStepMin);

                if (f_printLog.getValue() && d_logIterationSummary.getValue())
                    msg_info() << "[NCP RETRY " << continuationStages << "." << retry << "]"
                               << " eta=" << eta
                               << " step=" << previousStep << " -> " << lambdaStep;

                if (lambdaStep <= lambdaStepMin && previousStep <= lambdaStepMin)
                    break;
            }

            if (failed)
                break;

            if (!stageAccepted)
            {
                // Partial injection is not an accepted outer Newton update.
                ncp->restoreNewtonState();
                function.evaluateCurrentGuess();
                d_lastFailureReason.setValue("could not inject complete lambda target");
                failed = true;
                break;
            }
        }

        if (failed)
            break;

        // The full outer dLambda_k has now been injected and x has been
        // mechanically equilibrated at the final lambda target.
        ++acceptedOuterUpdates;
        d_lastAcceptedNewtonUpdates.setValue(acceptedOuterUpdates);
        d_lastAcceptedAlpha.setValue(1_sreal);

        function.evaluateCurrentGuess();
        residual2 = function.squaredNormLastEvaluation();
        summary = ncp->currentNCPDebugSummary();

        if (!finiteNonNegative(residual2) || !summary.valid)
        {
            d_lastFailureReason.setValue("invalid residual after complete lambda target");
            failed = true;
            break;
        }

        if (f_printLog.getValue() && d_logIterationSummary.getValue())
            msg_info() << "[NCP TARGET " << acceptedOuterUpdates << "]"
                       << " dLambda=" << rawLambda
                       << " stages=" << targetStages
                       << " R=" << normFromSquared(residual2)
                       << " M=" << summary.mechanicalResidualNorm
                       << " C=" << summary.complementarityResidualNorm;

        // This is the first full-NCP decision after continuation started.
        if (fullNCPConverged(summary, initialSummary, absoluteTolerance, relativeTolerance))
        {
            reportConverged(summary, residual2);
            return;
        }

        // Not converged: the next OUTER iteration rebuilds J and computes a new
        // dLambda target. continuationMode remains latched for this solve.
    }

    if (!failed && d_lastFailureReason.getValue().empty())
        d_lastFailureReason.setValue("maximum outer Newton iterations reached");

    if (!keepLastAcceptedState)
    {
        ncp->restoreSolveState();
        function.evaluateCurrentGuess();
        acceptedOuterUpdates = 0;
        continuationStages = 0;
        d_lastAcceptedNewtonUpdates.setValue(0u);
        d_lastAcceptedAlpha.setValue(0_sreal);
        d_lastContinuationStages.setValue(0u);
    }

    static constexpr NewtonStatus failedStatus("DivergedContinuation");
    d_status.setValue(failedStatus);

    msg_warning() << "[NCP FAIL] " << d_lastFailureReason.getValue()
                  << " outer=" << acceptedOuterUpdates
                  << " cont=" << continuationStages << ".";
}

} // namespace sofa::ncp
