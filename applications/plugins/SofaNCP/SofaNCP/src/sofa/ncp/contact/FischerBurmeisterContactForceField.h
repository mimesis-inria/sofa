/****************************************************************************
* Generic Fischer-Burmeister NCP contact force field.
*
* Residual:
*     R_x      += H^T lambda
*     R_lambda += phi(g, r lambda)
*
* With r frozen during a nonlinear solve:
*
*         [ lambda Hess(g)    H^T ]
*     J = [                       ]
*         [      a H          b   ]
*
* Compliance modes:
*   0 Fixed   : use fixedComplianceScale.
*   1 Lagged  : use per-contact r = H C_ref H^T, frozen for the whole nonlinear solve.
*               C_ref is built once from a symmetric positive reference beam metric
*               K_ref = sum(T0 K0 T0^T), after eliminating fixed DOFs.
*   2 Current : reserved for same-linearization scaling. Not enabled yet because
*               changing r inside a Newton pass requires a matching residual policy.
****************************************************************************/
#pragma once

#include <sofa/ncp/config.h>

#include <sofa/component/constraint/projective/FixedProjectiveConstraint.h>
#include <sofa/component/solidmechanics/fem/elastic/BeamFEMForceField.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/behavior/MixedInteractionForceField.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/type/RGBAColor.h>
#include <sofa/type/Mat.h>
#include <sofa/type/Vec.h>
#include <sofa/type/vector.h>

#include <cstddef>
#include <limits>

namespace sofa::ncp
{

enum class ContactRowStatus : unsigned int
{
    Active = 0,
    Pinned = 1,
    InvalidGeometry = 2
};

template<class DataTypes1, class DataTypes2>
struct ContactRow
{
    using Real = typename DataTypes1::Real;
    using Vec3 = sofa::type::Vec<3, Real>;

    sofa::Index pointIndex = 0;
    sofa::Index lambdaIndex = 0;
    ContactRowStatus status = ContactRowStatus::Pinned;

    Vec3 gapGradient = Vec3(Real(0), Real(0), Real(0));
    Real gap = Real(0);
    Real lambda = Real(0);
    Real complianceScale = Real(1);
    Real scaledLambda = Real(0);
    Real phi = Real(0);
    Real dPhiDgap = Real(0);
    Real dPhiDlambda = Real(1);
};

template<class TDataTypes1, class TDataTypes2>
class FischerBurmeisterContactForceField
    : public core::behavior::MixedInteractionForceField<TDataTypes1, TDataTypes2>
{
public:
    SOFA_ABSTRACT_CLASS(
        SOFA_TEMPLATE2(FischerBurmeisterContactForceField, TDataTypes1, TDataTypes2),
        SOFA_TEMPLATE2(core::behavior::MixedInteractionForceField, TDataTypes1, TDataTypes2));

    using Inherit = core::behavior::MixedInteractionForceField<TDataTypes1, TDataTypes2>;
    using DataTypes1 = TDataTypes1;
    using DataTypes2 = TDataTypes2;
    using Real = typename DataTypes1::Real;
    using Vec3 = sofa::type::Vec<3, Real>;
    using Mat3 = sofa::type::Mat<3, 3, Real>;
    using Coord1 = typename DataTypes1::Coord;
    using Deriv1 = typename DataTypes1::Deriv;
    using Deriv2 = typename DataTypes2::Deriv;
    using VecCoord1 = typename DataTypes1::VecCoord;
    using VecCoord2 = typename DataTypes2::VecCoord;
    using VecDeriv1 = typename DataTypes1::VecDeriv;
    using VecDeriv2 = typename DataTypes2::VecDeriv;
    using DataVecCoord1 = core::objectmodel::Data<VecCoord1>;
    using DataVecCoord2 = core::objectmodel::Data<VecCoord2>;
    using DataVecDeriv1 = core::objectmodel::Data<VecDeriv1>;
    using DataVecDeriv2 = core::objectmodel::Data<VecDeriv2>;
    using Contact = ContactRow<DataTypes1, DataTypes2>;
    using ContactStatus = ContactRowStatus;
    using FixedConstraint = sofa::component::constraint::projective::FixedProjectiveConstraint<DataTypes1>;
    using BeamForceField = sofa::component::solidmechanics::fem::elastic::BeamFEMForceField<sofa::defaulttype::Rigid3Types>;

    struct ResidualBlockNorms
    {
        SReal mechanicalSquaredNorm { 0_sreal };
        SReal complementaritySquaredNorm { 0_sreal };
    };

    struct ContactDiagnostics
    {
        sofa::Size activeCount { 0 };
        sofa::Size pinnedCount { 0 };
        sofa::Size invalidCount { 0 };
        Real minimumActiveGap { Real(0) };
        Real maximumPenetration { Real(0) };
        Real minimumLambda { Real(0) };
        Real maximumLambda { Real(0) };
        Real phiSquaredNorm { Real(0) };
        bool hasActiveContact { false };
    };

    static constexpr sofa::Size DerivDim1 = Deriv1::total_size;
    static constexpr sofa::Size DerivDim2 = Deriv2::total_size;
    static constexpr sofa::Size TranslationalDim = 3;

    Data<Real> d_fixedComplianceScale;
    Data<unsigned int> d_complianceMode;
    Data<Real> d_fbEpsilon;

    SingleLink<FischerBurmeisterContactForceField<DataTypes1, DataTypes2>, BeamForceField,
        BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_beamForceField;

    SingleLink<FischerBurmeisterContactForceField<DataTypes1, DataTypes2>, FixedConstraint,
        BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_fixedConstraint;

    Data<bool> d_showContactGradients;
    Data<Real> d_drawGradientScale;
    Data<sofa::type::RGBAColor> d_contactColor;

    Data<bool> d_debug;
    Data<bool> d_publishContactData;
    Data<sofa::type::vector<unsigned int>> d_contactStatus;
    Data<sofa::type::vector<Vec3>> d_contactGapGradient;
    Data<sofa::type::vector<Real>> d_contactGap;
    Data<sofa::type::vector<Real>> d_contactLambda;
    Data<sofa::type::vector<Real>> d_contactComplianceScale;
    Data<sofa::type::vector<Real>> d_contactScaledLambda;
    Data<sofa::type::vector<Real>> d_contactPhi;
    Data<sofa::type::vector<Real>> d_contactDPhiDgap;
    Data<sofa::type::vector<Real>> d_contactDPhiDlambda;
    Data<sofa::Size> d_activeContactCount;
    Data<sofa::Size> d_pinnedContactCount;
    Data<sofa::Size> d_invalidContactCount;

protected:
    FischerBurmeisterContactForceField();
    FischerBurmeisterContactForceField(core::behavior::MechanicalState<DataTypes1>* object1,
        core::behavior::MechanicalState<DataTypes2>* object2);

public:
    void init() override;
    void reinit() override;

    void addForce(const sofa::core::MechanicalParams*, DataVecDeriv1&, DataVecDeriv2&,
        const DataVecCoord1&, const DataVecCoord2&, const DataVecDeriv1&, const DataVecDeriv2&) override;

    void addDForce(const sofa::core::MechanicalParams*, DataVecDeriv1&, DataVecDeriv2&,
        const DataVecDeriv1&, const DataVecDeriv2&) override;

    void addKToMatrix(const sofa::core::MechanicalParams*, const sofa::core::behavior::MultiMatrixAccessor*) override;
    void buildStiffnessMatrix(core::behavior::StiffnessMatrix*) override;
    void buildDampingMatrix(core::behavior::DampingMatrix*) final;

    SReal getPotentialEnergy(const sofa::core::MechanicalParams*, const DataVecCoord1&, const DataVecCoord2&) const override;
    void draw(const core::visual::VisualParams*) override;

    bool isCurrentEvaluationValid() const { return m_validState; }
    ResidualBlockNorms currentResidualBlockNorms() const;
    ContactDiagnostics currentContactDiagnostics() const;

    bool usesLaggedCompliance() const;
    bool usesCurrentCompliance() const;

    // Lagged-scaling transaction. current r stays frozen between begin and commit/discard.
    void beginNonlinearSolve();
    bool commitLaggedComplianceSnapshot();
    void discardLaggedComplianceSnapshot();

    // Directional FD diagnostics use the already evaluated Newton base/trials.
    void storeFiniteDifferenceBase();
    void logFiniteDifferenceTrial(Real alpha) const;

protected:
    enum ComplianceMode : unsigned int
    {
        Fixed = 0,
        Lagged = 1,
        Current = 2
    };

    struct FiniteDifferenceContactSnapshot
    {
        Contact contact;
        Vec3 position = Vec3(Real(0), Real(0), Real(0));
    };

    sofa::type::vector<Contact> m_contacts;
    sofa::type::vector<FiniteDifferenceContactSnapshot> m_fdBaseContacts;

    // Used by the current nonlinear solve and never changed during that solve.
    sofa::type::vector<Real> m_currentCompliance;
    bool m_hasCurrentCompliance = false;

    // Accepted values waiting for the next nonlinear solve.
    sofa::type::vector<Real> m_nextCompliance;
    bool m_hasNextCompliance = false;

    // Constant diagonal translational blocks of K_ref^{-1}, one 3x3 block per beam node.
    sofa::type::vector<Mat3> m_referencePointCompliance;
    sofa::Size m_cachedReferenceMetricVersion = std::numeric_limits<sofa::Size>::max();
    std::size_t m_cachedConstraintSignature = 0;
    bool m_referenceComplianceCacheValid = false;
    sofa::Size m_complianceGeneration = 0;

    bool m_validState = false;

    virtual ContactStatus computeContactKinematics(const Vec3&, Contact&) const = 0;
    virtual bool computeGapHessian(const Vec3& p, Mat3& hessian) const;

    bool initializeContactRows();
    bool rebuildCurrentContacts(const VecCoord1&, const VecCoord2&);
    void finalizeContactRow(Contact&, ContactStatus, Real fixedComplianceScale) const;
    void updateFischerBurmeisterTerms(Contact&) const;

    Real complianceForPoint(sofa::Index pointIndex, Real fallback) const;

    void invalidateReferenceComplianceCache();
    bool ensureReferenceComplianceCache();
    bool rebuildReferenceComplianceCache();
    std::size_t fixedConstraintSignature() const;
    bool computeLaggedComplianceFromCurrentContacts(sofa::type::vector<Real>& candidate) const;

    ContactDiagnostics summarizeContacts() const;
    void publishDebugData();
    void drawActiveNormals(const core::visual::VisualParams*) const;

    static bool hasValidKinematics(const Contact&);
    static Vec3 extractPosition(const Coord1&);
    static Real positiveOrFallback(Real value, Real fallback);
};

} // namespace sofa::ncp