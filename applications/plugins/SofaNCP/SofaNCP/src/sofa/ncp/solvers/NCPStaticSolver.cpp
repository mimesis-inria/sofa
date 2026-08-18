#include <sofa/ncp/solvers/NCPStaticSolver.h>
#include <sofa/ncp/solvers/NCPDebugNewtonRaphsonSolver.h>
#include <sofa/ncp/solvers/NCPMechanicalContinuationFunctionInterface.h>
#include <sofa/ncp/solvers/NCPMechanicalContinuationNewtonRaphsonSolver.h>
#include <sofa/ncp/contact/FischerBurmeisterContactForceField.inl>

#include <sofa/component/odesolver/backward/NonLinearFunction.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/ScopedAdvancedTimer.h>
#include <sofa/helper/logging/Messaging.h>
#include <sofa/linearalgebra/BaseMatrix.h>
#include <sofa/linearalgebra/BaseVector.h>
#include <sofa/linearalgebra/FullVector.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/VectorOperations.h>

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace sofa::ncp
{

void registerNCPStaticSolver(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(core::ObjectRegistrationData(
        "Static Fischer-Burmeister backend with exact trial residuals and merit backtracking.")
        .add<NCPStaticSolver>());
}

NCPStaticSolver::NCPStaticSolver()
    : Inherit()
    , l_contactForceField(initLink("contactForceField", "Fischer-Burmeister contact force field."))
    , d_mechanicalResidualReference(initData(&d_mechanicalResidualReference, 1_sreal,
        "mechanicalResidualReference", "Positive mechanical residual scale."))
    , d_complementarityResidualReference(initData(&d_complementarityResidualReference, 1_sreal,
        "complementarityResidualReference", "Positive complementarity residual scale."))
    , d_debug(initData(&d_debug, false, "debug", "Print residual-evaluation summaries."))
    , d_finiteDifferenceCheck(initData(&d_finiteDifferenceCheck, false,
        "finiteDifferenceCheck", "Compare actual line-search residual derivatives against J*d and contact derivatives."))
    , d_lastMechanicalResidualNorm(initData(&d_lastMechanicalResidualNorm, 0_sreal,
        "lastMechanicalResidualNorm", "Mechanical residual norm from the retained state."))
    , d_lastComplementarityResidualNorm(initData(&d_lastComplementarityResidualNorm, 0_sreal,
        "lastComplementarityResidualNorm", "Complementarity residual norm from the retained state."))
    , d_lastComplianceCaptureSucceeded(initData(&d_lastComplianceCaptureSucceeded, false,
        "lastComplianceCaptureSucceeded", "Whether compliance was captured after the last converged solve."))
{
    d_lastMechanicalResidualNorm.setReadOnly(true);
    d_lastComplementarityResidualNorm.setReadOnly(true);
    d_lastComplianceCaptureSucceeded.setReadOnly(true);
}

void NCPStaticSolver::init()
{
    Inherit::init();

    if (!l_contactForceField)
    {
        l_contactForceField.set(getContext()->get<ContactForceField>(getContext()->getTags(),core::objectmodel::BaseContext::SearchDown));
    }

    auto* contact = l_contactForceField.get();
    auto* diagnosticNewton = dynamic_cast<NCPDebugNewtonRaphsonSolver*>(l_newtonSolver.get());

    if (!contact || !l_linearSolver.get() || !diagnosticNewton)
    {
        msg_error() << "NCPStaticSolver requires contactForceField, linearSolver and NCPDebugNewtonRaphsonSolver links.";
        d_componentState.setValue(core::objectmodel::ComponentState::Invalid);
        return;
    }

    if (!(d_mechanicalResidualReference.getValue() > 0_sreal)
        || !(d_complementarityResidualReference.getValue() > 0_sreal))
    {
        msg_error() << "mechanicalResidualReference and complementarityResidualReference must be positive.";
        d_componentState.setValue(core::objectmodel::ComponentState::Invalid);
        return;
    }

    d_componentState.setValue(core::objectmodel::ComponentState::Valid);
}

namespace
{

using BaseNonLinearFunction = sofa::component::odesolver::backward::newton_raphson::BaseNonLinearFunction;

class NCPStaticResidualFunction final
    : public BaseNonLinearFunction
    , public NCPDebugNewtonFunctionInterface
    , public NCPMechanicalContinuationFunctionInterface
{
public:
    using MechanicalOperations = sofa::simulation::common::MechanicalOperations;
    using ContactForceField = NCPStaticSolver::ContactForceField;

    NCPStaticResidualFunction(
        const NCPStaticSolver& owner,
        MechanicalOperations& mechanicalOperations,
        core::behavior::MultiVecCoord& position,
        core::behavior::MultiVecDeriv& residual,
        core::behavior::MultiVecDeriv& correction,
        core::behavior::LinearSolver* linearSolver,
        ContactForceField* contact,
        core::MultiVecCoordId solveStateId,
        core::MultiVecCoordId newtonStateId,
        core::MultiVecCoordId continuationStateId,
        core::MultiVecCoordId mechanicalStateId,
        SReal mechanicalReference,
        SReal complementarityReference,
        bool debug,
        bool finiteDifferenceCheck)
        : m_owner(owner)
        , m_mechanicalOperations(mechanicalOperations)
        , m_position(position)
        , m_residual(residual)
        , m_correction(correction)
        , m_linearSolver(linearSolver)
        , m_contact(contact)
        , m_solveState(position.ops(), solveStateId)
        , m_newtonState(position.ops(), newtonStateId)
        , m_continuationState(position.ops(), continuationStateId)
        , m_mechanicalState(position.ops(), mechanicalStateId)
        , m_mechanicalReference(mechanicalReference)
        , m_complementarityReference(complementarityReference)
        , m_debug(debug)
        , m_finiteDifferenceCheck(finiteDifferenceCheck)
    {
    }

    /** Enforce the current prescribed positions before the solve-entry checkpoint. */
    void projectCurrentState()
    {
        projectAndPropagate();
    }

    /** Evaluate the complete true residual at the current x and lambda. */
    void evaluateCurrentGuess() override
    {
        SCOPED_TIMER("NCPResidualEvaluation");
        m_mechanicalOperations.computeForce(m_residual, true, true);
        m_mechanicalOperations.projectResponse(m_residual);

        m_summary = NCPDebugResidualSummary{};

        if (!m_contact->isCurrentEvaluationValid())
        {
            m_scaledResidualSquaredNorm = std::numeric_limits<SReal>::infinity();
            m_summary.valid = false;
            return;
        }

        const auto blocks = m_contact->currentResidualBlockNorms();
        const auto diagnostics = m_contact->currentContactDiagnostics();
        const SReal rawSquaredNorm = m_residual.dot(m_residual);

        m_complementarityResidualSquaredNorm = blocks.complementaritySquaredNorm;
        m_mechanicalResidualSquaredNorm = std::max(rawSquaredNorm - m_complementarityResidualSquaredNorm,SReal(0));

        m_scaledResidualSquaredNorm =
            m_mechanicalResidualSquaredNorm / (m_mechanicalReference * m_mechanicalReference)
            + m_complementarityResidualSquaredNorm / (m_complementarityReference * m_complementarityReference);

        m_summary.mechanicalResidualNorm = std::sqrt(m_mechanicalResidualSquaredNorm);
        m_summary.complementarityResidualNorm = std::sqrt(m_complementarityResidualSquaredNorm);
        m_summary.scaledMechanicalSquaredNorm = m_mechanicalResidualSquaredNorm / (m_mechanicalReference * m_mechanicalReference);
        m_summary.scaledComplementaritySquaredNorm = m_complementarityResidualSquaredNorm / (m_complementarityReference * m_complementarityReference);
        m_summary.minimumActiveGap = diagnostics.minimumActiveGap;
        m_summary.maximumPenetration = diagnostics.maximumPenetration;
        m_summary.minimumLambda = diagnostics.minimumLambda;
        m_summary.maximumLambda = diagnostics.maximumLambda;
        m_summary.activeContacts = diagnostics.activeCount;
        m_summary.pinnedContacts = diagnostics.pinnedCount;
        m_summary.invalidContacts = diagnostics.invalidCount;
        m_summary.valid = diagnostics.invalidCount == 0 && std::isfinite(static_cast<double>(m_scaledResidualSquaredNorm));

    }

    SReal squaredNormLastEvaluation() override
    {
        return m_scaledResidualSquaredNorm;
    }

    /** Assemble A=-dR/dz at the current accepted Newton base. */
    void computeGradientFromCurrentGuess() override
    {
        assembleLinearizedSystem();
    }

    /** Solve the current mixed linear system and dispatch dz to all states. */
    void solveLinearEquation() override
    {
        SCOPED_TIMER("NCPLinearSolve");
        auto* system = m_linearSolver->getLinearSystem();
        system->setSystemSolution(m_correction);
        system->setRHS(m_residual);
        m_linearSolver->solveSystem();
        system->dispatchSystemSolution(m_correction);
    }

    // void solveLinearEquation() override
    // {
    //     SCOPED_TIMER("NCPLMSolve");

    //     auto* system = m_linearSolver->getLinearSystem();
    //     system->setSystemSolution(m_correction);
    //     system->setRHS(m_residual);

    //     auto* matrix = system->getSystemBaseMatrix();
    //     auto* rhs = system->getSystemRHSBaseVector();
    //     auto* solution = system->getSystemSolutionBaseVector();

    //     if (!matrix || !rhs || !solution)
    //         return;

    //     const Eigen::Index rows = matrix->rowSize();
    //     const Eigen::Index cols = matrix->colSize();

    //     Eigen::MatrixXd augmented = Eigen::MatrixXd::Zero(rows + cols, cols);
    //     Eigen::VectorXd augmentedRhs = Eigen::VectorXd::Zero(rows + cols);

    //     for (Eigen::Index row = 0; row < rows; ++row)
    //     {
    //         augmentedRhs[row] = rhs->element(row);

    //         for (Eigen::Index col = 0; col < cols; ++col)
    //             augmented(row, col) = matrix->element(row, col);
    //     }

    //     const Eigen::Index lambdaDofs = static_cast<Eigen::Index>(lambdaDofCount());
    //     const Eigen::Index lambdaBegin = cols - lambdaDofs;

    //     const SReal muMechanical = 1e-10_sreal;
    //     const SReal muLambda = 3e-5_sreal;

    //     for (Eigen::Index col = 0; col < cols; ++col)
    //     {
    //         const SReal mu = col < lambdaBegin ? muMechanical : muLambda;
    //         augmented(rows + col, col) = std::sqrt(mu);
    //     }

    //     const Eigen::VectorXd correction = augmented.colPivHouseholderQr().solve(augmentedRhs);

    //     for (Eigen::Index i = 0; i < cols; ++i)
    //         solution->set(i, static_cast<SReal>(correction[i]));

    //     system->dispatchSystemSolution(m_correction);
    // }

    /** Compatibility path; the merit solver reconstructs trials from m_newtonState. */
    void updateGuessFromLinearSolution(SReal alpha) override
    {
        m_position.peq(m_correction, alpha);
        projectAndPropagate();
    }

    SReal squaredNormDx() override
    {
        return m_correction.dot(m_correction);
    }

    SReal squaredLastEvaluation() override
    {
        return m_position.dot(m_position);
    }

    /** Store/restore the nonlinear solve-entry state. */
    void storeSolveState() override
    {
        m_solveState.eq(m_position.id());
    }

    void restoreSolveState() override
    {
        restoreState(m_solveState);
    }

    /** Store/restore the current accepted Newton base. */
    void storeNewtonState() override
    {
        m_newtonState.eq(m_position.id());
        m_cachedLambdaPredictor.clear();
    }

    void restoreNewtonState() override
    {
        restoreState(m_newtonState);
    }

    /** Construct z(alpha)=z_k+alpha*dz from the same accepted base every time. */
    void setTrialFromNewtonState(SReal alpha) override
    {
        m_position.eq(m_newtonState.id());
        m_position.peq(m_correction, alpha);
        projectAndPropagate();
    }

    NCPDebugResidualSummary currentNCPDebugSummary() const override
    {
        return m_summary;
    }

    NCPDebugDirectionSummary currentDirectionSummary() const override
    {
        NCPDebugDirectionSummary summary;

        auto* system = m_linearSolver ? m_linearSolver->getLinearSystem() : nullptr;
        auto* solution = system ? system->getSystemSolutionBaseVector() : nullptr;

        if (!solution)
            return summary;

        const sofa::SignedIndex totalDofs = solution->size();
        const sofa::SignedIndex lambdaDofs = static_cast<sofa::SignedIndex>(lambdaDofCount());
        const sofa::SignedIndex mechanicalDofs = totalDofs - lambdaDofs;

        if (totalDofs <= 0 || lambdaDofs < 0 || mechanicalDofs < 0)
            return summary;

        // Current NCP scene assumption:
        // one leading Rigid3 mechanical block: 6 scalar DOFs per node,
        // followed by one Vec1 multiplier DOF per contact row.
        if (mechanicalDofs % 6 != 0)
            return summary;

        SReal totalSquaredNorm = 0_sreal;
        SReal translationSquaredNorm = 0_sreal;
        SReal rotationSquaredNorm = 0_sreal;
        SReal lambdaSquaredNorm = 0_sreal;

        const sofa::SignedIndex rigidNodes = mechanicalDofs / 6;

        for (sofa::SignedIndex node = 0; node < rigidNodes; ++node)
        {
            const sofa::SignedIndex offset = 6 * node;

            SReal nodeTranslationSquaredNorm = 0_sreal;
            SReal nodeRotationSquaredNorm = 0_sreal;

            for (sofa::SignedIndex component = 0; component < 3; ++component)
            {
                const SReal value = solution->element(offset + component);

                if (!std::isfinite(static_cast<double>(value)))
                    return NCPDebugDirectionSummary{};

                nodeTranslationSquaredNorm += value * value;
            }

            for (sofa::SignedIndex component = 3; component < 6; ++component)
            {
                const SReal value = solution->element(offset + component);

                if (!std::isfinite(static_cast<double>(value)))
                    return NCPDebugDirectionSummary{};

                nodeRotationSquaredNorm += value * value;
            }

            translationSquaredNorm += nodeTranslationSquaredNorm;
            rotationSquaredNorm += nodeRotationSquaredNorm;

            summary.maximumNodeTranslation =
                std::max(summary.maximumNodeTranslation, std::sqrt(nodeTranslationSquaredNorm));

            summary.maximumNodeRotation =
                std::max(summary.maximumNodeRotation, std::sqrt(nodeRotationSquaredNorm));
        }

        for (sofa::SignedIndex i = mechanicalDofs; i < totalDofs; ++i)
        {
            const SReal value = solution->element(i);

            if (!std::isfinite(static_cast<double>(value)))
                return NCPDebugDirectionSummary{};

            lambdaSquaredNorm += value * value;
            summary.maximumAbsLambda = std::max(summary.maximumAbsLambda, std::abs(value));
        }

        totalSquaredNorm =
            translationSquaredNorm
            + rotationSquaredNorm
            + lambdaSquaredNorm;

        summary.totalNorm = std::sqrt(totalSquaredNorm);
        summary.translationNorm = std::sqrt(translationSquaredNorm);
        summary.rotationNorm = std::sqrt(rotationSquaredNorm);
        summary.lambdaNorm = std::sqrt(lambdaSquaredNorm);

        summary.rigidNodes = static_cast<sofa::Size>(rigidNodes);
        summary.lambdaDofs = static_cast<sofa::Size>(lambdaDofs);
        summary.valid = true;

        return summary;
    }

    bool correctionIsFinite() const override
    {
        const SReal norm2 = m_correction.dot(m_correction);
        return std::isfinite(static_cast<double>(norm2)) && norm2 >= 0_sreal;
    }

    void beginFiniteDifferenceCheck() override
    {
        m_fdReady = false;
        m_fdDetailedLogged = false;

        if (!m_finiteDifferenceCheck)
            return;

        auto* system = m_linearSolver ? m_linearSolver->getLinearSystem() : nullptr;
        auto* matrix = system ? system->getSystemBaseMatrix() : nullptr;
        auto* rhs = system ? system->getSystemRHSBaseVector() : nullptr;
        auto* solution = system ? system->getSystemSolutionBaseVector() : nullptr;

        if (!matrix || !rhs || !solution)
            return;

        const sofa::SignedIndex rows = matrix->rowSize();
        const sofa::SignedIndex cols = matrix->colSize();
        if (rows <= 0 || rows != cols || rhs->size() != rows || solution->size() != cols)
            return;

        const sofa::Size lambdaDofs = lambdaDofCount();
        if (lambdaDofs > static_cast<sofa::Size>(rows))
            return;

        m_fdMechanicalDofs = static_cast<sofa::Size>(rows) - lambdaDofs;
        m_fdBaseResidual.resize(static_cast<std::size_t>(rows));
        m_fdJd.resize(static_cast<std::size_t>(rows));

        sofa::linearalgebra::FullVector<double> direction(cols);
        sofa::linearalgebra::FullVector<double> Ad(rows);
        for (sofa::SignedIndex col = 0; col < cols; ++col)
            direction[col] = static_cast<double>(solution->element(col));

        matrix->opMulV(&Ad, &direction);

        SReal linearDefect2 = 0_sreal;
        SReal rhs2 = 0_sreal;

        for (sofa::SignedIndex row = 0; row < rows; ++row)
        {
            const SReal base = rhs->element(row);
            const SReal matrixDirection = static_cast<SReal>(Ad[row]);

            m_fdBaseResidual[static_cast<std::size_t>(row)] = base;
            m_fdJd[static_cast<std::size_t>(row)] = -matrixDirection; // A=-J.

            const SReal defect = matrixDirection - base; // linear system is A*d=R.
            linearDefect2 += defect * defect;
            rhs2 += base * base;
        }

        const SReal denom = std::max(std::sqrt(rhs2), std::numeric_limits<SReal>::epsilon());
        msg_info(&m_owner) << "[NCP FD BASE] rows=" << rows
                           << " mechanicalRows=" << m_fdMechanicalDofs
                           << " lambdaRows=" << lambdaDofs
                           << " linearRelErr=" << std::sqrt(linearDefect2) / denom;

        m_contact->storeFiniteDifferenceBase();
        m_fdReady = true;
    }

    void evaluateFiniteDifferenceTrial(SReal alpha, unsigned int iteration, unsigned int trial) override
    {
        if (!m_finiteDifferenceCheck || !m_fdReady || !(alpha > 0_sreal))
            return;

        auto* system = m_linearSolver ? m_linearSolver->getLinearSystem() : nullptr;
        if (!system)
            return;

        // Flatten the already evaluated true trial residual in the same global ordering.
        system->setRHS(m_residual);
        auto* rhs = system->getSystemRHSBaseVector();
        if (!rhs || static_cast<std::size_t>(rhs->size()) != m_fdBaseResidual.size())
            return;

        SReal error2 = 0_sreal;
        SReal fd2 = 0_sreal;
        SReal predicted2 = 0_sreal;
        SReal mechanicalError2 = 0_sreal;
        SReal mechanicalFD2 = 0_sreal;
        SReal mechanicalPredicted2 = 0_sreal;
        SReal complementarityError2 = 0_sreal;
        SReal complementarityFD2 = 0_sreal;
        SReal complementarityPredicted2 = 0_sreal;

        // One row/contact dump per Newton iteration, at a genuinely local step.
        // Global norms are still printed for every trial.
        const bool detailed = !m_fdDetailedLogged && alpha <= 1e-3_sreal;

        for (sofa::SignedIndex row = 0; row < rhs->size(); ++row)
        {
            const SReal base = m_fdBaseResidual[static_cast<std::size_t>(row)];
            const SReal current = rhs->element(row);
            const SReal fd = (current - base) / alpha;
            const SReal predicted = m_fdJd[static_cast<std::size_t>(row)];
            const SReal error = fd - predicted;
            const SReal predictedTrial = base + alpha * predicted;
            const bool mechanical = static_cast<sofa::Size>(row) < m_fdMechanicalDofs;

            if (detailed)
            {
                msg_info(&m_owner) << "[NCP FD ROW] it=" << iteration
                                   << " trial=" << trial
                                   << " alpha=" << alpha
                                   << " row=" << row
                                   << " block=" << (mechanical ? "M" : "FB")
                                   << " R0=" << base
                                   << " Rtrial=" << current
                                   << " Rpred=" << predictedTrial
                                   << " dRfd=" << fd
                                   << " Jd=" << predicted
                                   << " err=" << error;
            }

            error2 += error * error;
            fd2 += fd * fd;
            predicted2 += predicted * predicted;

            if (mechanical)
            {
                mechanicalError2 += error * error;
                mechanicalFD2 += fd * fd;
                mechanicalPredicted2 += predicted * predicted;
            }
            else
            {
                complementarityError2 += error * error;
                complementarityFD2 += fd * fd;
                complementarityPredicted2 += predicted * predicted;
            }
        }

        const SReal eps = std::numeric_limits<SReal>::epsilon();
        const SReal relScale = std::max({std::sqrt(fd2), std::sqrt(predicted2), eps});
        const SReal mechScale = std::max({std::sqrt(mechanicalFD2), std::sqrt(mechanicalPredicted2), eps});
        const SReal fbScale = std::max({std::sqrt(complementarityFD2), std::sqrt(complementarityPredicted2), eps});
        const SReal relErr = std::sqrt(error2) / relScale;
        const SReal mechRelErr = std::sqrt(mechanicalError2) / mechScale;
        const SReal fbRelErr = std::sqrt(complementarityError2) / fbScale;

        msg_info(&m_owner) << "[NCP FD GLOBAL] it=" << iteration
                           << " trial=" << trial
                           << " alpha=" << alpha
                           << " relErr=" << relErr
                           << " mechRelErr=" << mechRelErr
                           << " fbRelErr=" << fbRelErr
                           << " |dRfd|=" << std::sqrt(fd2)
                           << " |Jd|=" << std::sqrt(predicted2);

        m_contact->logFiniteDifferenceTrial(static_cast<ContactForceField::Real>(alpha));
        if (detailed)
        {
            m_fdDetailedLogged = true;
        }
    }

    /**
     * Select the line-search direction exactly in the spirit of Algorithm 4.1.
     *
     * The static system assembled by this class is
     *
     *     A = -J,
     *
     * because setSystemMBKMatrix() uses stiffnessFactor=-1 while the nonlinear
     * residual is R(z). For
     *
     *     Psi(z) = 0.5 R(z)^T W^2 R(z),
     *
     * we therefore have
     *
     *     -grad(Psi) = -J^T W^2 R = A^T W^2 R.
     *
     * SOFA already stores both A and the flattened RHS R in the linear system,
     * so the fallback requires only one transpose matrix-vector product.
     *
     * IMPORTANT:
     * This is the true merit gradient only when the assembled J is the true
     * derivative of evaluateCurrentGuess(). This is exact for the affine plane
     * contact with compliance frozen during the nonlinear solve. For curved
     * geometries, add the missing lambda*Hess(g) and any dr/dx term first.
     */
    bool selectMeritSearchDirection(SReal rho,SReal exponent,SReal& meritSlope,bool& usedGradientFallback) override
    {
        SOFA_UNUSED(rho);
        SOFA_UNUSED(exponent);

        meritSlope = std::numeric_limits<SReal>::quiet_NaN();
        usedGradientFallback = false;

        auto* system = m_linearSolver ? m_linearSolver->getLinearSystem() : nullptr;
        auto* matrix = system ? system->getSystemBaseMatrix() : nullptr;
        auto* rhs = system ? system->getSystemRHSBaseVector() : nullptr;
        auto* solution = system ? system->getSystemSolutionBaseVector() : nullptr;

        if (!system || !matrix || !rhs || !solution)
            return false;

        const auto rows = matrix->rowSize();
        const auto cols = matrix->colSize();

        if (rows <= 0
            || rows != cols
            || rhs->size() != rows
            || solution->size() != cols)
        {
            return false;
        }

        // One Vec1 multiplier DOF exists per contact row. In the present plugin
        // architecture the multiplier block is the trailing block of the global
        // system. The three status counts always sum to the number of rows.
        const auto diagnostics = m_contact->currentContactDiagnostics();
        const sofa::Size lambdaDofs = diagnostics.activeCount + diagnostics.pinnedCount + diagnostics.invalidCount;

        if (lambdaDofs > static_cast<sofa::Size>(rows))
            return false;

        const sofa::Size mechanicalDofs = static_cast<sofa::Size>(rows) - lambdaDofs;

        const SReal invMechanicalReference2 = 1_sreal / (m_mechanicalReference * m_mechanicalReference);
        const SReal invComplementarityReference2 = 1_sreal / (m_complementarityReference * m_complementarityReference);

        sofa::linearalgebra::FullVector<double> weightedResidual(rows);
        sofa::linearalgebra::FullVector<double> negativeGradient(cols);

        for (sofa::SignedIndex i = 0; i < rows; ++i)
        {
            const SReal weight = static_cast<sofa::Size>(i) < mechanicalDofs
                                    ? invMechanicalReference2
                                    : invComplementarityReference2;

            const SReal value = weight * rhs->element(i);
            if (!std::isfinite(static_cast<double>(value)))
                return false;

            weightedResidual[i] = static_cast<double>(value);
        }

        // A = -J, hence this is -grad(Psi).
        matrix->opMulTV(&negativeGradient, &weightedResidual);

        SReal gradientNorm2 = 0_sreal;
        SReal newtonNorm2 = 0_sreal;
        SReal negativeGradientDotNewton = 0_sreal;

        for (sofa::SignedIndex i = 0; i < cols; ++i)
        {
            const SReal minusGradient = static_cast<SReal>(negativeGradient[i]);
            const SReal newtonDirection = solution->element(i);

            if (!std::isfinite(static_cast<double>(minusGradient))
                || !std::isfinite(static_cast<double>(newtonDirection)))
            {
                return false;
            }

            gradientNorm2 += minusGradient * minusGradient;
            newtonNorm2 += newtonDirection * newtonDirection;
            negativeGradientDotNewton += minusGradient * newtonDirection;
        }

        if (!std::isfinite(static_cast<double>(gradientNorm2))
            || !std::isfinite(static_cast<double>(newtonNorm2)))
        {
            return false;
        }

        // Since negativeGradient=-grad(Psi):
        //
        //     grad(Psi)^T d_N = -negativeGradient^T d_N.
        const SReal newtonMeritSlope = -negativeGradientDotNewton;

        const SReal newtonNorm = std::sqrt(std::max(newtonNorm2, 0_sreal));
        // const SReal requiredSlope = -rho * std::pow(newtonNorm, exponent);
        const SReal requiredSlope = 0.0;

        if (m_debug)
        {
            msg_info(&m_owner) << "[NCP DESCENT]"
                               << " slope=" << newtonMeritSlope
                               << " |d|=" << newtonNorm
                               << " Criterion=" << requiredSlope
                               << " |gradPsi|=" << std::sqrt(gradientNorm2)
                               << " descent=" << (newtonMeritSlope < 0_sreal);
        }

        if (std::isfinite(static_cast<double>(newtonMeritSlope)) && newtonMeritSlope <= requiredSlope)
        {
            meritSlope = newtonMeritSlope;
            usedGradientFallback = false;

            // The current system solution already contains d_N and was already
            // dispatched to m_correction by solveLinearEquation().
            return true;
        }

        // Newton direction does not have enough descent. Use
        //
        //     d_G = -grad(Psi).
        if (!(gradientNorm2 > std::numeric_limits<SReal>::epsilon()))
            return false;

        for (sofa::SignedIndex i = 0; i < cols; ++i)
            solution->set(i, static_cast<SReal>(negativeGradient[i]));

        system->dispatchSystemSolution(m_correction);

        if (m_debug)
        {
            msg_info(&m_owner) << "[NCP FALLBACK]"
                               << " |gradPsi|=" << std::sqrt(gradientNorm2)
                               << " slope=" << -gradientNorm2
                               << " |d|=" << std::sqrt(m_correction.dot(m_correction));
        }

        meritSlope = -gradientNorm2;
        usedGradientFallback = correctionIsFinite() && std::isfinite(static_cast<double>(meritSlope)) && meritSlope < 0_sreal;

        return usedGradientFallback;
    }

    SReal maxAbsLambdaCorrection() const override
    {
        auto* system = m_linearSolver ? m_linearSolver->getLinearSystem() : nullptr;
        auto* solution = system ? system->getSystemSolutionBaseVector() : nullptr;
        if (!solution)
            return std::numeric_limits<SReal>::infinity();

        const sofa::SignedIndex lambdaDofs = static_cast<sofa::SignedIndex>(lambdaDofCount());
        const sofa::SignedIndex size = solution->size();
        if (lambdaDofs <= 0 || lambdaDofs > size)
            return std::numeric_limits<SReal>::infinity();

        SReal maxValue = 0_sreal;
        for (sofa::SignedIndex i = size - lambdaDofs; i < size; ++i)
            maxValue = std::max(maxValue, std::abs(solution->element(i)));

        return maxValue;
    }

    bool cacheLambdaPredictor() override
    {
        auto* system = m_linearSolver ? m_linearSolver->getLinearSystem() : nullptr;
        auto* solution = system ? system->getSystemSolutionBaseVector() : nullptr;
        if (!solution)
            return false;

        const sofa::SignedIndex lambdaDofs = static_cast<sofa::SignedIndex>(lambdaDofCount());
        const sofa::SignedIndex size = solution->size();
        const sofa::SignedIndex lambdaBegin = size - lambdaDofs;
        if (lambdaDofs <= 0 || lambdaBegin < 0)
            return false;

        m_cachedLambdaPredictor.resize(static_cast<std::size_t>(lambdaDofs));

        for (sofa::SignedIndex i = 0; i < lambdaDofs; ++i)
        {
            const SReal value = solution->element(lambdaBegin + i);
            if (!std::isfinite(static_cast<double>(value)))
                return false;

            m_cachedLambdaPredictor[static_cast<std::size_t>(i)] = value;
        }

        return true;
    }

    void storeContinuationState() override
    {
        m_continuationState.eq(m_position.id());
    }

    void restoreContinuationState() override
    {
        restoreState(m_continuationState);
        evaluateCurrentGuess();
    }

    bool applyLambdaPredictorFraction(SReal fraction, SReal& appliedInfinityNorm) override
    {
        appliedInfinityNorm = 0_sreal;

        if (!(fraction > 0_sreal) || fraction > 1_sreal || m_cachedLambdaPredictor.empty())
            return false;

        auto* system = m_linearSolver ? m_linearSolver->getLinearSystem() : nullptr;
        auto* solution = system ? system->getSystemSolutionBaseVector() : nullptr;
        if (!solution)
            return false;

        const sofa::SignedIndex lambdaDofs = static_cast<sofa::SignedIndex>(m_cachedLambdaPredictor.size());
        const sofa::SignedIndex size = solution->size();
        const sofa::SignedIndex lambdaBegin = size - lambdaDofs;
        if (lambdaDofs <= 0 || lambdaBegin < 0)
            return false;

        for (sofa::SignedIndex i = 0; i < lambdaBegin; ++i)
            solution->set(i, 0_sreal);

        for (sofa::SignedIndex i = 0; i < lambdaDofs; ++i)
        {
            const SReal increment = fraction * m_cachedLambdaPredictor[static_cast<std::size_t>(i)];
            solution->set(lambdaBegin + i, increment);
            appliedInfinityNorm = std::max(appliedInfinityNorm, std::abs(increment));
        }

        system->dispatchSystemSolution(m_correction);

        // Important: increment the CURRENT accepted continuation state. Do not
        // reconstruct from m_newtonState; m_newtonState is the fixed outer base.
        m_position.peq(m_correction, 1_sreal);
        projectAndPropagate();
        evaluateCurrentGuess();

        return m_summary.valid;
    }

    void storeMechanicalState() override
    {
        m_mechanicalState.eq(m_position.id());
    }

    void restoreMechanicalState() override
    {
        restoreState(m_mechanicalState);
        evaluateCurrentGuess();
    }

    bool solveMechanicalLinearEquation() override
    {
        SCOPED_TIMER("NCPMechanicalLinearSolve");

        assembleLinearizedSystem();

        auto* system = m_linearSolver ? m_linearSolver->getLinearSystem() : nullptr;
        if (!system)
            return false;

        system->setSystemSolution(m_correction);
        system->setRHS(m_residual);

        auto* matrix = system->getSystemBaseMatrix();
        auto* rhs = system->getSystemRHSBaseVector();
        auto* solution = system->getSystemSolutionBaseVector();
        if (!matrix || !rhs || !solution)
            return false;

        const sofa::SignedIndex lambdaDofs = static_cast<sofa::SignedIndex>(lambdaDofCount());
        const sofa::SignedIndex size = rhs->size();
        const sofa::SignedIndex lambdaBegin = size - lambdaDofs;

        if (lambdaDofs <= 0 || lambdaBegin < 0 || matrix->rowSize() != size || matrix->colSize() != size)
            return false;

        /**
         * Fixed-lambda mechanical Newton:
         *
         * [ A_xx  A_xλ ] [dx]   [R_x]
         * [   0     I  ] [dλ] = [ 0 ].
         *
         * Only lambda ROWS are replaced. A_xλ is deliberately kept: the first
         * residual equation still contains the current H^T lambda force, while
         * dλ=0 makes A_xλ*dλ vanish in this inner correction.
         */
        matrix->clearRows(lambdaBegin, size);
        for (sofa::SignedIndex i = lambdaBegin; i < size; ++i)
        {
            matrix->set(i, i, 1.0);
            rhs->set(i, 0_sreal);
        }

        matrix->compress();
        m_linearSolver->solveSystem();

        // Enforce the mathematical constraint exactly before dispatch, even if
        // the linear solver leaves tiny round-off values in the identity rows.
        for (sofa::SignedIndex i = lambdaBegin; i < size; ++i)
            solution->set(i, 0_sreal);

        system->dispatchSystemSolution(m_correction);
        return correctionIsFinite();
    }

    void setMechanicalTrialFromState(SReal alpha) override
    {
        m_position.eq(m_mechanicalState.id());
        m_position.peq(m_correction, alpha);
        projectAndPropagate();
    }

    void evaluateMechanicalTrial() override
    {
        evaluateCurrentGuess();
    }

    SReal scaledMechanicalSquaredNorm() const override
    {
        return m_summary.scaledMechanicalSquaredNorm;
    }

    bool continuationStateValid() const override
    {
        return m_summary.valid && std::isfinite(static_cast<double>(m_summary.scaledMechanicalSquaredNorm));
    }

    SReal mechanicalResidualNorm() const
    {
        return std::sqrt(std::max(m_mechanicalResidualSquaredNorm, SReal(0)));
    }

    SReal complementarityResidualNorm() const
    {
        return std::sqrt(std::max(m_complementarityResidualSquaredNorm, SReal(0)));
    }

private:
    const NCPStaticSolver& m_owner;
    MechanicalOperations& m_mechanicalOperations;
    core::behavior::MultiVecCoord& m_position;
    core::behavior::MultiVecDeriv& m_residual;
    core::behavior::MultiVecDeriv& m_correction;
    core::behavior::LinearSolver* m_linearSolver = nullptr;
    ContactForceField* m_contact = nullptr;
    core::behavior::MultiVecCoord m_solveState;
    core::behavior::MultiVecCoord m_newtonState;
    core::behavior::MultiVecCoord m_continuationState;
    core::behavior::MultiVecCoord m_mechanicalState;
    std::vector<SReal> m_cachedLambdaPredictor;
    const SReal m_mechanicalReference;
    const SReal m_complementarityReference;
    const bool m_debug;
    const bool m_finiteDifferenceCheck;

    std::vector<SReal> m_fdBaseResidual;
    std::vector<SReal> m_fdJd;
    sofa::Size m_fdMechanicalDofs = 0;
    bool m_fdReady = false;
    bool m_fdDetailedLogged = false;

    NCPDebugResidualSummary m_summary;
    SReal m_scaledResidualSquaredNorm = std::numeric_limits<SReal>::infinity();
    SReal m_mechanicalResidualSquaredNorm = std::numeric_limits<SReal>::infinity();
    SReal m_complementarityResidualSquaredNorm = std::numeric_limits<SReal>::infinity();

    sofa::Size lambdaDofCount() const
    {
        const auto diagnostics = m_contact->currentContactDiagnostics();
        return diagnostics.activeCount + diagnostics.pinnedCount + diagnostics.invalidCount;
    }

    void assembleLinearizedSystem()
    {
        SCOPED_TIMER("NCPJacobianAssembly");
        static constexpr core::MatricesFactors::M massFactor(0);
        static constexpr core::MatricesFactors::B dampingFactor(0);
        static constexpr core::MatricesFactors::K stiffnessFactor(-1);
        m_mechanicalOperations.setSystemMBKMatrix(massFactor,dampingFactor,stiffnessFactor,m_linearSolver);
    }

    /** Checkpoints are already projected; restoration only propagates mappings. */
    void restoreState(const core::behavior::MultiVecCoord& state)
    {
        m_position.eq(state.id());
        m_mechanicalOperations.propagateX(m_position);
    }

    void projectAndPropagate()
    {
        m_mechanicalOperations.solveConstraint(m_position,core::ConstraintOrder::POS);
        m_mechanicalOperations.propagateX(m_position);
    }
};

} // namespace

void NCPStaticSolver::solve(const core::ExecParams* params,SReal dt,core::MultiVecCoordId xResult,core::MultiVecDerivId vResult)
{
    if (!isComponentStateValid())
        return;

    SOFA_UNUSED(dt);
    SOFA_UNUSED(vResult);

    sofa::simulation::common::VectorOperations vectorOperations(params, getContext());
    sofa::simulation::common::MechanicalOperations mechanicalOperations(params, getContext());
    mechanicalOperations->setImplicit(true);

    core::behavior::MultiVecCoord position(&vectorOperations, xResult);
    core::behavior::MultiVecDeriv residual(&vectorOperations, core::vec_id::write_access::force);
    core::behavior::MultiVecDeriv correction(&vectorOperations, core::vec_id::write_access::dx);
    correction.realloc(&vectorOperations, true, true);

    auto prepareCheckpoint = [&](core::MultiVecCoordId& id, const char* name)
    {
        core::behavior::MultiVecCoord state(&vectorOperations, id);
        state.realloc( &vectorOperations, true, true, core::VecIdProperties(name, getClassName()));
        id = state.id();
    };

    prepareCheckpoint(m_solveStateId, "ncpSolveState");
    prepareCheckpoint(m_newtonStateId, "ncpNewtonState");
    prepareCheckpoint(m_continuationStateId, "ncpContinuationState");
    prepareCheckpoint(m_mechanicalStateId, "ncpMechanicalState");

    auto* contact = l_contactForceField.get();
    contact->beginNonlinearSolve();

    NCPStaticResidualFunction residualFunction(
        *this,
        mechanicalOperations,
        position,
        residual,
        correction,
        l_linearSolver.get(),
        l_contactForceField.get(),
        m_solveStateId,
        m_newtonStateId,
        m_continuationStateId,
        m_mechanicalStateId,
        d_mechanicalResidualReference.getValue(),
        d_complementarityResidualReference.getValue(),
        d_debug.getValue(),
        d_finiteDifferenceCheck.getValue());

    residualFunction.projectCurrentState();

    auto* newton = dynamic_cast<NCPDebugNewtonRaphsonSolver*>(l_newtonSolver.get());

    if (auto* continuationNewton = dynamic_cast<NCPMechanicalContinuationNewtonRaphsonSolver*>(newton))
        continuationNewton->solveNCPContinuation(residualFunction);
    else
        newton->solveNCP(residualFunction);

    d_lastMechanicalResidualNorm.setValue(residualFunction.mechanicalResidualNorm());
    d_lastComplementarityResidualNorm.setValue(residualFunction.complementarityResidualNorm());

    bool complianceReady = false;
    if (newton->lastSolveConverged())
        complianceReady = contact->commitLaggedComplianceSnapshot();
    else
        contact->discardLaggedComplianceSnapshot();

    d_lastComplianceCaptureSucceeded.setValue(complianceReady);
}

} // namespace sofa::ncp
