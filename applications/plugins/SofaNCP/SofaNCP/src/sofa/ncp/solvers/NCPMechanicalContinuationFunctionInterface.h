#pragma once

#include <sofa/ncp/config.h>

namespace sofa::ncp
{

/**
 * Hooks implemented by the NCP static residual backend for target-tracked,
 * fixed-lambda mechanical continuation.
 */
class SOFANCP_API NCPMechanicalContinuationFunctionInterface
{
public:
    virtual ~NCPMechanicalContinuationFunctionInterface() = default;

    /// Infinity norm of the lambda block of the current mixed Newton correction.
    virtual SReal maxAbsLambdaCorrection() const = 0;

    /// Cache the current mixed Newton lambda predictor before any merit fallback can overwrite it.
    virtual bool cacheLambdaPredictor() = 0;

    /// Checkpoint the last accepted partial continuation state.
    virtual void storeContinuationState() = 0;
    virtual void restoreContinuationState() = 0;

    /**
     * Add fraction * dLambda_k to the CURRENT accepted continuation state.
     *
     * dLambda_k is the cached outer mixed predictor and remains fixed until the
     * whole target has been injected. fraction is incremental, not cumulative.
     */
    virtual bool applyLambdaPredictorFraction(SReal fraction, SReal& appliedInfinityNorm) = 0;

    /// Inner fixed-lambda mechanical Newton checkpoint.
    virtual void storeMechanicalState() = 0;
    virtual void restoreMechanicalState() = 0;

    /// Assemble the current tangent and solve dx while enforcing dLambda = 0.
    virtual bool solveMechanicalLinearEquation() = 0;

    /// Mechanical line-search trial from the same inner checkpoint.
    virtual void setMechanicalTrialFromState(SReal alpha) = 0;
    virtual void evaluateMechanicalTrial() = 0;

    /// Dimensionless ||Rx / mechanicalResidualReference||^2.
    virtual SReal scaledMechanicalSquaredNorm() const = 0;
    virtual bool continuationStateValid() const = 0;
};

} // namespace sofa::ncp
