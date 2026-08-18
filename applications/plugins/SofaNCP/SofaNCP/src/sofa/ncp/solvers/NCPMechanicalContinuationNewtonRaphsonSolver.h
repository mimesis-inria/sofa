#pragma once

#include <sofa/ncp/config.h>
#include <sofa/ncp/solvers/NCPDebugNewtonRaphsonSolver.h>

namespace sofa::core { class ObjectFactory; }

namespace sofa::ncp
{

/**
 * Full mixed FB Newton with a fixed-lambda mechanical corrector used only when
 * the mixed predictor requests an excessive multiplier increment.
 */
class SOFANCP_API NCPMechanicalContinuationNewtonRaphsonSolver : public NCPDebugNewtonRaphsonSolver
{
public:
    SOFA_CLASS(NCPMechanicalContinuationNewtonRaphsonSolver, NCPDebugNewtonRaphsonSolver);

    using Inherit = NCPDebugNewtonRaphsonSolver;
    using BaseNonLinearFunction = sofa::component::odesolver::backward::newton_raphson::BaseNonLinearFunction;

    NCPMechanicalContinuationNewtonRaphsonSolver();
    void reset() override;
    void solveNCPContinuation(BaseNonLinearFunction& function);

    Data<SReal> d_lambdaStepInitial;
    Data<SReal> d_lambdaStepMin;
    Data<SReal> d_lambdaStepMax;
    Data<SReal> d_lambdaStepGrowth;
    Data<SReal> d_lambdaStepShrink;
    Data<unsigned int> d_maxMechanicalIterations;
    Data<unsigned int> d_maxMechanicalLineSearchTrials;
    Data<unsigned int> d_maxLambdaStageRetries;
    Data<unsigned int> d_maxContinuationStages;

    Data<SReal> d_lastLambdaStage;
    Data<unsigned int> d_lastMechanicalIterations;
    Data<unsigned int> d_lastContinuationStages;
};

SOFANCP_API void registerNCPMechanicalContinuationNewtonRaphsonSolver(sofa::core::ObjectFactory* factory);

} // namespace sofa::ncp
