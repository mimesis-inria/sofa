/****************************************************************************
* Generic Fischer-Burmeister NCP contact force field implementation.
****************************************************************************/
#pragma once

#include <sofa/ncp/contact/FischerBurmeisterContactForceField.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/core/behavior/BaseLocalForceFieldMatrix.h>
#include <sofa/core/behavior/MultiMatrixAccessor.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/logging/Messaging.h>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

namespace sofa::ncp
{

template<class T1, class T2>
FischerBurmeisterContactForceField<T1, T2>::FischerBurmeisterContactForceField()
    : Inherit()
    , d_fixedComplianceScale(initData(&d_fixedComplianceScale, Real(1.0),
          "fixedComplianceScale", "Positive fallback r used in phi(g,r*lambda)."))
    , d_complianceMode(initData(&d_complianceMode, static_cast<unsigned int>(Fixed),
          "complianceMode", "0=fixed, 1=lagged reference-elastic scale, 2=current scale (reserved)."))
    , d_fbEpsilon(initData(&d_fbEpsilon, Real(0),
          "fbEpsilon", "Nonnegative Fischer-Burmeister smoothing epsilon."))
    , l_beamForceField(initLink("beamForceField",
          "BeamFEMForceField providing the constant positive reference elastic metric."))
    , l_fixedConstraint(initLink("fixedConstraint",
          "Fixed projective constraint eliminated from the reference elastic metric."))
    , d_showContactGradients(initData(&d_showContactGradients, true,
          "showContactGradients", "Draw active near-contact gap gradients."))
    , d_drawGradientScale(initData(&d_drawGradientScale, Real(5),
          "drawGradientScale", "Contact-gradient drawing scale."))
    , d_contactColor(initData(&d_contactColor, sofa::type::RGBAColor(1.0f, 0.85f, 0.0f, 1.0f),
          "contactColor", "Contact-gradient drawing color."))
    , d_debug(initData(&d_debug, true, "debug", "Publish compact contact/compliance diagnostics."))
    , d_publishContactData(initData(&d_publishContactData, true,
          "publishContactData", "Publish point-indexed contact diagnostics."))
    , d_contactStatus(initData(&d_contactStatus, "contactStatus", "0=Active, 1=Pinned, 2=InvalidGeometry."))
    , d_contactGapGradient(initData(&d_contactGapGradient, "contactGapGradient", "Point-indexed H=dg/dx."))
    , d_contactGap(initData(&d_contactGap, "contactGap", "Point-indexed gaps."))
    , d_contactLambda(initData(&d_contactLambda, "contactLambda", "Point-indexed multipliers."))
    , d_contactComplianceScale(initData(&d_contactComplianceScale, "contactComplianceScale", "Point-indexed scalar r."))
    , d_contactScaledLambda(initData(&d_contactScaledLambda, "contactScaledLambda", "Point-indexed r*lambda."))
    , d_contactPhi(initData(&d_contactPhi, "contactPhi", "Point-indexed row residuals."))
    , d_contactDPhiDgap(initData(&d_contactDPhiDgap, "contactDPhiDgap", "Point-indexed dphi/dg."))
    , d_contactDPhiDlambda(initData(&d_contactDPhiDlambda, "contactDPhiDlambda", "Point-indexed dphi/dlambda."))
    , d_activeContactCount(initData(&d_activeContactCount, sofa::Size(0), "activeContactCount", "Active FB rows."))
    , d_pinnedContactCount(initData(&d_pinnedContactCount, sofa::Size(0), "pinnedContactCount", "Pinned rows."))
    , d_invalidContactCount(initData(&d_invalidContactCount, sofa::Size(0), "invalidContactCount", "Invalid rows."))
{
    static_assert(T1::spatial_dimensions >= TranslationalDim,
        "Object1 must expose at least three translational DOFs.");
    static_assert(T2::spatial_dimensions == 1,
        "Object2 must be a scalar Vec1-like multiplier state.");

    d_contactStatus.setReadOnly(true);
    d_contactGapGradient.setReadOnly(true);
    d_contactGap.setReadOnly(true);
    d_contactLambda.setReadOnly(true);
    d_contactComplianceScale.setReadOnly(true);
    d_contactScaledLambda.setReadOnly(true);
    d_contactPhi.setReadOnly(true);
    d_contactDPhiDgap.setReadOnly(true);
    d_contactDPhiDlambda.setReadOnly(true);
    d_activeContactCount.setReadOnly(true);
    d_pinnedContactCount.setReadOnly(true);
    d_invalidContactCount.setReadOnly(true);
}

template<class T1, class T2>
FischerBurmeisterContactForceField<T1, T2>::FischerBurmeisterContactForceField(
    core::behavior::MechanicalState<DataTypes1>* object1,
    core::behavior::MechanicalState<DataTypes2>* object2)
    : FischerBurmeisterContactForceField()
{
    this->mstate1.set(object1);
    this->mstate2.set(object2);
}

template<class T1, class T2>
bool FischerBurmeisterContactForceField<T1, T2>::initializeContactRows()
{
    if (!this->mstate1 || !this->mstate2)
    {
        msg_error() << "object1 and object2 MechanicalState links are required.";
        m_validState = false;
        return false;
    }

    if (this->mstate1->getSize() != this->mstate2->getSize())
    {
        msg_error() << "object1 size (" << this->mstate1->getSize()
                    << ") and lambda size (" << this->mstate2->getSize()
                    << ") must match.";
        m_validState = false;
        return false;
    }

    const unsigned int mode = d_complianceMode.getValue();
    if (mode > static_cast<unsigned int>(Current))
    {
        msg_error() << "complianceMode must be 0=fixed, 1=lagged or 2=current.";
        m_validState = false;
        return false;
    }

    if (mode == static_cast<unsigned int>(Current))
    {
        msg_error() << "complianceMode=2 is reserved. Same-linearization compliance requires a base-residual refresh.";
        m_validState = false;
        return false;
    }

    if (!(d_fixedComplianceScale.getValue() > Real(0))
        || !std::isfinite(static_cast<double>(d_fixedComplianceScale.getValue())))
    {
        msg_error() << "fixedComplianceScale must be finite and positive.";
        m_validState = false;
        return false;
    }

    if (!(d_fbEpsilon.getValue() >= Real(0))
        || !std::isfinite(static_cast<double>(d_fbEpsilon.getValue())))
    {
        msg_error() << "fbEpsilon must be finite and nonnegative.";
        m_validState = false;
        return false;
    }

    m_contacts.clear();
    m_validState = false;
    return true;
}

template<class T1, class T2>
void FischerBurmeisterContactForceField<T1, T2>::init()
{
    Inherit::init();

    if (!l_beamForceField)
    {
        l_beamForceField.set(this->getContext()->template get<BeamForceField>(
            this->getContext()->getTags(), core::objectmodel::BaseContext::SearchDown));
    }

    if (!l_fixedConstraint)
    {
        l_fixedConstraint.set(this->getContext()->template get<FixedConstraint>(
            this->getContext()->getTags(), core::objectmodel::BaseContext::SearchDown));
    }

    const bool valid = initializeContactRows();

    if (valid && usesLaggedCompliance())
    {
        if constexpr (!std::is_same_v<DataTypes1, sofa::defaulttype::Rigid3Types>)
        {
            msg_error() << "complianceMode=1 currently requires object1=Rigid3.";
            this->d_componentState.setValue(core::objectmodel::ComponentState::Invalid);
            return;
        }

        if (!l_beamForceField)
        {
            msg_error() << "complianceMode=1 requires a BeamFEMForceField link.";
            this->d_componentState.setValue(core::objectmodel::ComponentState::Invalid);
            return;
        }

        if (!l_fixedConstraint)
        {
            msg_error() << "complianceMode=1 requires a FixedProjectiveConstraint link.";
            this->d_componentState.setValue(core::objectmodel::ComponentState::Invalid);
            return;
        }
    }

    this->d_componentState.setValue(valid
        ? core::objectmodel::ComponentState::Valid
        : core::objectmodel::ComponentState::Invalid);
}

template<class T1, class T2>
void FischerBurmeisterContactForceField<T1, T2>::reinit()
{
    Inherit::reinit();

    // Keep the cached metric/scales when beam K0 and constrained indices did not change.
    // beginNonlinearSolve() checks both cache signatures before the next promotion.
    const bool valid = initializeContactRows();

    this->d_componentState.setValue(valid
        ? core::objectmodel::ComponentState::Valid
        : core::objectmodel::ComponentState::Invalid);
}

template<class T1, class T2>
bool FischerBurmeisterContactForceField<T1, T2>::usesLaggedCompliance() const
{
    return d_complianceMode.getValue() == static_cast<unsigned int>(Lagged);
}

template<class T1, class T2>
bool FischerBurmeisterContactForceField<T1, T2>::usesCurrentCompliance() const
{
    return d_complianceMode.getValue() == static_cast<unsigned int>(Current);
}

template<class T1, class T2>
void FischerBurmeisterContactForceField<T1, T2>::beginNonlinearSolve()
{
    if (!usesLaggedCompliance())
        return;

    if (!ensureReferenceComplianceCache())
    {
        msg_warning() << "[NCP SCALE] Reference elastic metric unavailable; keeping current/fixed r.";
        return;
    }

    const sofa::Size expectedPoints = this->mstate1 ? this->mstate1->getSize() : 0;

    if (m_hasNextCompliance && m_nextCompliance.size() == expectedPoints)
    {
        m_currentCompliance = std::move(m_nextCompliance);
        m_nextCompliance.clear();
        m_hasNextCompliance = false;
        m_hasCurrentCompliance = true;
        ++m_complianceGeneration;

        if (d_debug.getValue())
            msg_info() << "[NCP SCALE PROMOTE] generation=" << m_complianceGeneration
                       << " points=" << m_currentCompliance.size();
    }

    if (m_hasCurrentCompliance && m_currentCompliance.size() != expectedPoints)
    {
        m_currentCompliance.clear();
        m_hasCurrentCompliance = false;
    }
}

template<class T1, class T2>
void FischerBurmeisterContactForceField<T1, T2>::discardLaggedComplianceSnapshot()
{
    m_nextCompliance.clear();
    m_hasNextCompliance = false;

    if (d_debug.getValue() && usesLaggedCompliance())
        msg_info() << "[NCP SCALE DISCARD]";
}

template<class T1, class T2>
typename FischerBurmeisterContactForceField<T1, T2>::Real
FischerBurmeisterContactForceField<T1, T2>::positiveOrFallback(Real value, Real fallback)
{
    return std::isfinite(static_cast<double>(value)) && value > Real(0) ? value : fallback;
}

template<class T1, class T2>
typename FischerBurmeisterContactForceField<T1, T2>::Vec3
FischerBurmeisterContactForceField<T1, T2>::extractPosition(const Coord1& coordinate)
{
    const auto p = DataTypes1::getCPos(coordinate);
    return Vec3(static_cast<Real>(p[0]), static_cast<Real>(p[1]), static_cast<Real>(p[2]));
}

template<class T1, class T2>
bool FischerBurmeisterContactForceField<T1, T2>::computeGapHessian(const Vec3&, Mat3& hessian) const
{
    hessian.clear();
    return false;
}

template<class T1, class T2>
bool FischerBurmeisterContactForceField<T1, T2>::hasValidKinematics(const Contact& c)
{
    const Real norm2 = c.gapGradient.norm2();
    return std::isfinite(static_cast<double>(c.gap))
        && std::isfinite(static_cast<double>(c.lambda))
        && std::isfinite(static_cast<double>(c.gapGradient[0]))
        && std::isfinite(static_cast<double>(c.gapGradient[1]))
        && std::isfinite(static_cast<double>(c.gapGradient[2]))
        && std::isfinite(static_cast<double>(norm2))
        && norm2 > Real(1e-30);
}

template<class T1, class T2>
void FischerBurmeisterContactForceField<T1, T2>::updateFischerBurmeisterTerms(Contact& c) const
{
    const Real eps = d_fbEpsilon.getValue();
    const Real r = c.complianceScale;
    const Real contactCompliance = 0;

    // Compliant unilateral gap:
    
    //     g_eff = g + c_n * lambda
    
    // Active contact therefore tends toward
    
    //     g + c_n * lambda = 0
    
    // or
    
    //     lambda = -g / c_n.
    const Real effectiveGap = c.gap + contactCompliance * c.lambda;
    const Real s = r * c.lambda;
    const Real n = std::sqrt(effectiveGap * effectiveGap+ s * s+ eps * eps);

    c.scaledLambda = s;

    const Real phi = effectiveGap + s - n;

    // d phi / d g_eff
    const Real a = Real(1) - effectiveGap / n;

    // Contribution from s = r * lambda.
    const Real b = r * (Real(1) - s / n);
    const Real invR = Real(1) / r;

    c.phi = phi * invR;
    c.dPhiDgap = a * invR;
    c.dPhiDlambda = (a * contactCompliance + b) * invR;
}

template<class T1, class T2>
typename FischerBurmeisterContactForceField<T1, T2>::Real
FischerBurmeisterContactForceField<T1, T2>::complianceForPoint(sofa::Index pointIndex, Real fallback) const
{
    if (!usesLaggedCompliance() || !m_hasCurrentCompliance || pointIndex >= m_currentCompliance.size())
        return fallback;

    return positiveOrFallback(m_currentCompliance[pointIndex], fallback);
}

template<class T1, class T2>
void FischerBurmeisterContactForceField<T1, T2>::finalizeContactRow(Contact& c, ContactStatus geometryStatus, Real fixedR) const
{
    if (geometryStatus == ContactStatus::Pinned)
    {
        c.status = ContactStatus::Pinned;
        c.gap = Real(0);
        c.gapGradient.clear();
        c.complianceScale = Real(1);
        c.scaledLambda = c.lambda;
        c.phi = c.lambda;
        c.dPhiDgap = Real(0);
        c.dPhiDlambda = Real(1);
        return;
    }

    if (geometryStatus != ContactStatus::Active || !hasValidKinematics(c))
    {
        c.status = ContactStatus::InvalidGeometry;
        c.gap = Real(0);
        c.gapGradient.clear();
        c.complianceScale = Real(1);
        c.scaledLambda = c.lambda;
        c.phi = c.lambda;
        c.dPhiDgap = Real(0);
        c.dPhiDlambda = Real(1);
        return;
    }

    c.status = ContactStatus::Active;
    c.complianceScale = complianceForPoint(c.pointIndex, fixedR);
    updateFischerBurmeisterTerms(c);
}

template<class T1, class T2>
bool FischerBurmeisterContactForceField<T1, T2>::rebuildCurrentContacts(const VecCoord1& x1, const VecCoord2& x2)
{
    if (x1.size() != x2.size() || x1.empty())
    {
        msg_error() << "Cannot build contact rows: object1 has " << x1.size()
                    << " entries and object2 has " << x2.size() << ".";
        m_contacts.clear();
        m_validState = false;
        return false;
    }

    if (m_contacts.size() != x1.size())
        m_contacts.resize(x1.size());

    const Real fixedR = d_fixedComplianceScale.getValue();
    bool valid = true;

    for (sofa::Index i = 0; i < x1.size(); ++i)
    {
        Contact& c = m_contacts[i];
        c = Contact{};
        c.pointIndex = i;
        c.lambdaIndex = i;
        c.lambda = static_cast<Real>(x2[i][0]);

        const ContactStatus status = computeContactKinematics(extractPosition(x1[i]), c);
        finalizeContactRow(c, status, fixedR);
        valid = valid && c.status != ContactStatus::InvalidGeometry;
    }

    m_validState = valid;
    publishDebugData();
    return valid;
}

template<class T1, class T2>
void FischerBurmeisterContactForceField<T1, T2>::invalidateReferenceComplianceCache()
{
    m_referencePointCompliance.clear();
    m_cachedReferenceMetricVersion = std::numeric_limits<sofa::Size>::max();
    m_cachedConstraintSignature = 0;
    m_referenceComplianceCacheValid = false;
}

template<class T1, class T2>
std::size_t FischerBurmeisterContactForceField<T1, T2>::fixedConstraintSignature() const
{
    const auto* constraint = l_fixedConstraint.get();
    if (!constraint)
        return 0;

    std::size_t hash = static_cast<std::size_t>(1469598103934665603ULL);

    auto mix = [&hash](std::size_t value)
    {
        hash ^= value + static_cast<std::size_t>(0x9e3779b97f4a7c15ULL) + (hash << 6) + (hash >> 2);
    };

    mix(constraint->fixAllDOFs() ? 1u : 0u);

    for (const sofa::Index point : constraint->d_indices.getValue())
        mix(static_cast<std::size_t>(point));

    return hash;
}

template<class T1, class T2>
bool FischerBurmeisterContactForceField<T1, T2>::ensureReferenceComplianceCache()
{
    if (!this->mstate1 || !l_beamForceField || !l_fixedConstraint)
        return false;

    const auto* beam = l_beamForceField.get();
    const sofa::Size version = beam->getReferenceElasticMetricVersion();
    const std::size_t constraintSignature = fixedConstraintSignature();
    const sofa::Size pointCount = this->mstate1->getSize();

    if (m_referenceComplianceCacheValid
        && m_cachedReferenceMetricVersion == version
        && m_cachedConstraintSignature == constraintSignature
        && m_referencePointCompliance.size() == pointCount)
    {
        return true;
    }

    return rebuildReferenceComplianceCache();
}

template<class T1, class T2>
bool FischerBurmeisterContactForceField<T1, T2>::rebuildReferenceComplianceCache()
{
    if constexpr (!std::is_same_v<DataTypes1, sofa::defaulttype::Rigid3Types>)
    {
        return false;
    }
    else
    {
        auto* beam = l_beamForceField.get();
        auto* constraint = l_fixedConstraint.get();

        if (!beam || !constraint || !this->mstate1)
            return false;

        if (constraint->fixAllDOFs())
        {
            msg_error() << "[NCP SCALE] All mechanical points are fixed; reference compliance is undefined.";
            return false;
        }

        const sofa::Size pointCount = this->mstate1->getSize();
        const sofa::Size objectDofs = pointCount * DerivDim1;

        if (pointCount == 0 || DerivDim1 != 6)
            return false;

        sofa::type::vector<bool> constrainedPoint(pointCount, false);

        for (const sofa::Index point : constraint->d_indices.getValue())
        {
            if (point >= pointCount)
            {
                msg_warning() << "[NCP SCALE] Ignoring constrained point=" << point
                              << " pointCount=" << pointCount;
                continue;
            }

            constrainedPoint[point] = true;
        }

        // Full scalar DoF -> reduced free DoF.
        sofa::type::vector<int> freeIndex(objectDofs, -1);
        int freeDofCount = 0;

        for (sofa::Index point = 0; point < pointCount; ++point)
        {
            if (constrainedPoint[point])
                continue;

            for (sofa::Size d = 0; d < DerivDim1; ++d)
                freeIndex[point * DerivDim1 + d] = freeDofCount++;
        }

        if (freeDofCount == 0)
            return false;

        using SparseMatrix = Eigen::SparseMatrix<Real>;
        using Triplet = Eigen::Triplet<Real>;
        using DenseMatrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

        const sofa::Size elementCount = beam->getReferenceElasticMetricElementCount();
        std::vector<Triplet> triplets;
        triplets.reserve(static_cast<std::size_t>(elementCount) * 12u * 12u);

        typename BeamForceField::StiffnessMatrix Ke;

        for (sofa::Size metricElement = 0; metricElement < elementCount; ++metricElement)
        {
            sofa::Index a = 0;
            sofa::Index b = 0;

            if (!beam->getReferenceElasticMetricElement(metricElement, a, b, Ke))
                return false;

            if (a >= pointCount || b >= pointCount)
            {
                msg_error() << "[NCP SCALE] Beam metric element references node outside object1: "
                            << a << ", " << b << " pointCount=" << pointCount;
                return false;
            }

            sofa::Index global[12];

            for (sofa::Size d = 0; d < 6; ++d)
            {
                global[d] = a * 6 + d;
                global[6 + d] = b * 6 + d;
            }

            for (sofa::Size row = 0; row < 12; ++row)
            {
                const int reducedRow = freeIndex[global[row]];
                if (reducedRow < 0)
                    continue;

                for (sofa::Size col = 0; col < 12; ++col)
                {
                    const int reducedCol = freeIndex[global[col]];
                    if (reducedCol < 0)
                        continue;

                    // K0 is symmetric; averaging removes only round-off asymmetry
                    // introduced by the block rotations.
                    const Real value = Real(0.5) * (Ke(row, col) + Ke(col, row));

                    if (value != Real(0))
                        triplets.emplace_back(reducedRow, reducedCol, value);
                }
            }
        }

        SparseMatrix Kfree(freeDofCount, freeDofCount);
        Kfree.setFromTriplets(
            triplets.begin(), triplets.end(),
            [](const Real& a, const Real& b) { return a + b; });
        Kfree.makeCompressed();

        Eigen::SimplicialLDLT<SparseMatrix> factorization;
        factorization.compute(Kfree);

        if (factorization.info() != Eigen::Success)
        {
            msg_error() << "[NCP SCALE] Reference elastic metric factorization failed.";
            return false;
        }

        const auto D = factorization.vectorD();

        Real minD = std::numeric_limits<Real>::infinity();
        Real maxD = Real(0);

        for (Eigen::Index i = 0; i < D.size(); ++i)
        {
            const Real value = D[i];

            if (!std::isfinite(static_cast<double>(value)) || value <= Real(0))
            {
                msg_error() << "[NCP SCALE] Reference elastic metric is not positive definite. "
                            << "D[" << i << "]=" << value;
                return false;
            }

            minD = std::min(minD, value);
            maxD = std::max(maxD, value);
        }

        // Build all translational unit loads once. The resulting 3x3 diagonal
        // inverse block at each point is enough for every future H C_i H^T.
        sofa::type::vector<int> basisColumn(pointCount * TranslationalDim, -1);
        int basisCount = 0;

        for (sofa::Index point = 0; point < pointCount; ++point)
        {
            for (sofa::Size d = 0; d < TranslationalDim; ++d)
            {
                const int reduced = freeIndex[point * DerivDim1 + d];

                if (reduced >= 0)
                    basisColumn[point * TranslationalDim + d] = basisCount++;
            }
        }

        DenseMatrix rhs(freeDofCount, basisCount);
        rhs.setZero();

        for (sofa::Index point = 0; point < pointCount; ++point)
        {
            for (sofa::Size d = 0; d < TranslationalDim; ++d)
            {
                const int reduced = freeIndex[point * DerivDim1 + d];
                const int column = basisColumn[point * TranslationalDim + d];

                if (reduced >= 0 && column >= 0)
                    rhs(reduced, column) = Real(1);
            }
        }

        const DenseMatrix response = factorization.solve(rhs);

        if (factorization.info() != Eigen::Success || !response.allFinite())
        {
            msg_error() << "[NCP SCALE] Reference elastic metric solve failed.";
            return false;
        }

        sofa::type::vector<Mat3> pointCompliance(pointCount);

        for (sofa::Index point = 0; point < pointCount; ++point)
        {
            Mat3 C;
            C.clear();

            for (sofa::Size row = 0; row < TranslationalDim; ++row)
            {
                const int reducedRow = freeIndex[point * DerivDim1 + row];

                if (reducedRow < 0)
                    continue;

                for (sofa::Size col = 0; col < TranslationalDim; ++col)
                {
                    const int column = basisColumn[point * TranslationalDim + col];

                    if (column >= 0)
                        C(row, col) = response(reducedRow, column);
                }
            }

            // Preserve a symmetric positive metric despite solve round-off.
            for (sofa::Size row = 0; row < TranslationalDim; ++row)
            {
                for (sofa::Size col = row + 1; col < TranslationalDim; ++col)
                {
                    const Real average = Real(0.5) * (C(row, col) + C(col, row));
                    C(row, col) = average;
                    C(col, row) = average;
                }
            }

            pointCompliance[point] = C;
        }

        m_referencePointCompliance = std::move(pointCompliance);
        m_cachedReferenceMetricVersion = beam->getReferenceElasticMetricVersion();
        m_cachedConstraintSignature = fixedConstraintSignature();
        m_referenceComplianceCacheValid = true;

        // Values computed from an older metric must never survive a rebuild.
        m_currentCompliance.clear();
        m_nextCompliance.clear();
        m_hasCurrentCompliance = false;
        m_hasNextCompliance = false;

        if (d_debug.getValue())
        {
            msg_info() << "[NCP SCALE CACHE] elements=" << elementCount
                       << " points=" << pointCount
                       << " freeDofs=" << freeDofCount
                       << " translationalRhs=" << basisCount
                       << " LDLT_D[min/max]=" << minD << "/" << maxD
                       << " metricVersion=" << m_cachedReferenceMetricVersion;
        }

        return true;
    }
}

template<class T1, class T2>
bool FischerBurmeisterContactForceField<T1, T2>::computeLaggedComplianceFromCurrentContacts(
    sofa::type::vector<Real>& candidate) const
{
    if (!m_referenceComplianceCacheValid || !this->mstate1)
        return false;

    const sofa::Size pointCount = this->mstate1->getSize();

    if (m_referencePointCompliance.size() != pointCount)
        return false;

    const Real fallback = d_fixedComplianceScale.getValue();

    if (m_hasCurrentCompliance && m_currentCompliance.size() == pointCount)
        candidate = m_currentCompliance;
    else
        candidate.assign(pointCount, fallback);

    Real minR = std::numeric_limits<Real>::infinity();
    Real maxR = Real(0);
    Real sumR = Real(0);
    sofa::Size activeCount = 0;

    for (const Contact& c : m_contacts)
    {
        if (c.status != ContactStatus::Active)
            continue;

        if (c.pointIndex >= m_referencePointCompliance.size())
            return false;

        const Mat3& C = m_referencePointCompliance[c.pointIndex];
        Vec3 Ch;
        Ch.clear();

        for (sofa::Size row = 0; row < TranslationalDim; ++row)
            for (sofa::Size col = 0; col < TranslationalDim; ++col)
                Ch[row] += C(row, col) * c.gapGradient[col];

        const Real r = c.gapGradient * Ch;

        if (!std::isfinite(static_cast<double>(r)) || r <= Real(0))
        {
            msg_info() << "[NCP SCALE] Non-positive reference compliance at point="
                          << c.pointIndex << " r=" << r
                          << ". Keeping the previous/fixed scale for the next solve.";
            continue;
        }

        candidate[c.pointIndex] = r;
        minR = std::min(minR, r);
        maxR = std::max(maxR, r);
        sumR += r;
        ++activeCount;
    }

    if (d_debug.getValue())
    {
        if (activeCount > 0)
        {
            msg_info() << "[NCP SCALE COMPUTE] active=" << activeCount
                       << " r[min/mean/max]=" << minR << "/"
                       << (sumR / static_cast<Real>(activeCount)) << "/" << maxR;
        }
        else
        {
            msg_info() << "[NCP SCALE COMPUTE] active=0";
        }
    }

    return true;
}

template<class T1, class T2>
bool FischerBurmeisterContactForceField<T1, T2>::commitLaggedComplianceSnapshot()
{
    if (!usesLaggedCompliance())
        return false;

    if (!ensureReferenceComplianceCache())
        return false;

    sofa::type::vector<Real> candidate;

    if (!computeLaggedComplianceFromCurrentContacts(candidate))
        return false;

    m_nextCompliance = std::move(candidate);
    m_hasNextCompliance = true;

    if (d_debug.getValue())
        msg_info() << "[NCP SCALE COMMIT] points=" << m_nextCompliance.size()
                   << " targetGeneration=" << (m_complianceGeneration + 1);

    return true;
}

template<class T1, class T2>
typename FischerBurmeisterContactForceField<T1, T2>::ResidualBlockNorms
FischerBurmeisterContactForceField<T1, T2>::currentResidualBlockNorms() const
{
    ResidualBlockNorms result;

    if (this->mstate1)
    {
        const auto force = this->mstate1->read(core::vec_id::read_access::force);
        for (const Deriv1& value : force->getValue())
            for (sofa::Size d = 0; d < DerivDim1; ++d)
                result.mechanicalSquaredNorm += static_cast<SReal>(value[d]) * static_cast<SReal>(value[d]);
    }

    if (this->mstate2)
    {
        const auto force = this->mstate2->read(core::vec_id::read_access::force);
        for (const Deriv2& value : force->getValue())
            for (sofa::Size d = 0; d < DerivDim2; ++d)
                result.complementaritySquaredNorm += static_cast<SReal>(value[d]) * static_cast<SReal>(value[d]);
    }

    return result;
}

template<class T1, class T2>
typename FischerBurmeisterContactForceField<T1, T2>::ContactDiagnostics
FischerBurmeisterContactForceField<T1, T2>::summarizeContacts() const
{
    ContactDiagnostics out;
    Real minGap = std::numeric_limits<Real>::infinity();
    Real minLambda = std::numeric_limits<Real>::infinity();
    Real maxLambda = -std::numeric_limits<Real>::infinity();

    for (const Contact& c : m_contacts)
    {
        switch (c.status)
        {
        case ContactStatus::Active:
            ++out.activeCount;
            out.hasActiveContact = true;
            minGap = std::min(minGap, c.gap);
            break;
        case ContactStatus::Pinned:
            ++out.pinnedCount;
            break;
        case ContactStatus::InvalidGeometry:
            ++out.invalidCount;
            break;
        }

        minLambda = std::min(minLambda, c.lambda);
        maxLambda = std::max(maxLambda, c.lambda);
        out.phiSquaredNorm += c.phi * c.phi;
    }

    if (out.hasActiveContact)
    {
        out.minimumActiveGap = minGap;
        out.maximumPenetration = std::max(Real(0), -minGap);
    }

    if (!m_contacts.empty())
    {
        out.minimumLambda = minLambda;
        out.maximumLambda = maxLambda;
    }

    return out;
}

template<class T1, class T2>
typename FischerBurmeisterContactForceField<T1, T2>::ContactDiagnostics
FischerBurmeisterContactForceField<T1, T2>::currentContactDiagnostics() const
{
    return summarizeContacts();
}

template<class T1, class T2>
void FischerBurmeisterContactForceField<T1, T2>::storeFiniteDifferenceBase()
{
    m_fdBaseContacts.clear();

    if (!m_validState || !this->mstate1)
        return;

    const auto x1Data = this->mstate1->read(core::vec_id::read_access::position);
    const VecCoord1& x1 = x1Data->getValue();
    if (x1.size() != m_contacts.size())
        return;

    m_fdBaseContacts.resize(m_contacts.size());
    for (sofa::Index i = 0; i < m_contacts.size(); ++i)
    {
        m_fdBaseContacts[i].contact = m_contacts[i];
        m_fdBaseContacts[i].position = extractPosition(x1[i]);
    }
}

template<class T1, class T2>
void FischerBurmeisterContactForceField<T1, T2>::logFiniteDifferenceTrial(Real alpha) const
{
    if (!(alpha > Real(0)) || !this->mstate1 || m_fdBaseContacts.size() != m_contacts.size())
        return;

    const auto x1Data = this->mstate1->read(core::vec_id::read_access::position);
    const VecCoord1& x1 = x1Data->getValue();
    if (x1.size() != m_contacts.size())
        return;

    const Real eps = Real(1e-14);

    for (sofa::Index i = 0; i < m_contacts.size(); ++i)
    {
        const auto& snapshot = m_fdBaseContacts[i];
        const Contact& base = snapshot.contact;
        const Contact& trial = m_contacts[i];

        if (base.status != trial.status)
        {
            msg_info() << "[NCP FD CONTACT] point=" << i
                       << " alpha=" << alpha
                       << " status=" << static_cast<unsigned int>(base.status)
                       << "->" << static_cast<unsigned int>(trial.status);
            continue;
        }

        if (base.status != ContactStatus::Active)
            continue;

        const Real H0Norm = std::sqrt(base.gapGradient.norm2());
        const Real H1Norm = std::sqrt(trial.gapGradient.norm2());

        Vec3 deltaH;
        deltaH.clear();

        Real Hdot = Real(0);

        for (sofa::Size d = 0; d < TranslationalDim; ++d)
        {
            deltaH[d] = trial.gapGradient[d] - base.gapGradient[d];
            Hdot += base.gapGradient[d] * trial.gapGradient[d];
        }

        const Real deltaHNorm = std::sqrt(deltaH.norm2());

        const Real gradientRelativeChange =
            deltaHNorm / std::max(H0Norm, eps);

        Real normalCosine = Real(1);

        if (H0Norm > eps && H1Norm > eps)
        {
            normalCosine = Hdot / (H0Norm * H1Norm);
            normalCosine = std::max(Real(-1), std::min(Real(1), normalCosine));
        }

        const Real normalAngleDeg =
            std::acos(normalCosine) * Real(57.29577951308232);

        // Norm of the spatial part a*H of the FB row.
        const Real spatialRowNorm0 =
            std::abs(base.dPhiDgap) * H0Norm;

        const Real spatialRowNorm1 =
            std::abs(trial.dPhiDgap) * H1Norm;

        // Full FB row norm:
        //     [ dPhi/dg * H , dPhi/dlambda ]
        const Real fbRowNorm0 = std::sqrt(
            spatialRowNorm0 * spatialRowNorm0
            + base.dPhiDlambda * base.dPhiDlambda);

        const Real fbRowNorm1 = std::sqrt(
            spatialRowNorm1 * spatialRowNorm1
            + trial.dPhiDlambda * trial.dPhiDlambda);

        msg_info() << "[NCP H CONTINUITY]"
                << " point=" << i
                << " alpha=" << alpha
                << " H0=" << base.gapGradient
                << " H1=" << trial.gapGradient
                << " |H0|=" << H0Norm
                << " |H1|=" << H1Norm
                << " |dH|=" << deltaHNorm
                << " relDH=" << gradientRelativeChange
                << " angleDeg=" << normalAngleDeg
                << " a0=" << base.dPhiDgap
                << " a1=" << trial.dPhiDgap
                << " b0=" << base.dPhiDlambda
                << " b1=" << trial.dPhiDlambda
                << " |aH0|=" << spatialRowNorm0
                << " |aH1|=" << spatialRowNorm1
                << " |FBrow0|=" << fbRowNorm0
                << " |FBrow1|=" << fbRowNorm1;

        const Vec3 trialPosition = extractPosition(x1[i]);
        Vec3 dx, dHfd, dHj, forceBase, forceTrial, dForceFD, dForceJ;
        dx.clear(); dHfd.clear(); dHj.clear(); forceBase.clear(); forceTrial.clear(); dForceFD.clear(); dForceJ.clear();

        Real dgJ = Real(0);
        for (sofa::Size d = 0; d < TranslationalDim; ++d)
        {
            dx[d] = (trialPosition[d] - snapshot.position[d]) / alpha;
            dHfd[d] = (trial.gapGradient[d] - base.gapGradient[d]) / alpha;
            forceBase[d] = base.lambda * base.gapGradient[d];
            forceTrial[d] = trial.lambda * trial.gapGradient[d];
            dForceFD[d] = (forceTrial[d] - forceBase[d]) / alpha;
            dgJ += base.gapGradient[d] * dx[d];
        }

        const Real dgFD = (trial.gap - base.gap) / alpha;
        const Real dLambda = (trial.lambda - base.lambda) / alpha;
        const Real dr = (trial.complianceScale - base.complianceScale) / alpha;

        Mat3 hessian;
        const bool hasHessian = computeGapHessian(snapshot.position, hessian);
        if (hasHessian)
        {
            for (sofa::Size row = 0; row < TranslationalDim; ++row)
                for (sofa::Size col = 0; col < TranslationalDim; ++col)
                    dHj[row] += hessian(row, col) * dx[col];
        }

        for (sofa::Size d = 0; d < TranslationalDim; ++d)
            dForceJ[d] = base.gapGradient[d] * dLambda + base.lambda * dHj[d];

        const Real dPhiFD = (trial.phi - base.phi) / alpha;
        const Real dPhiJ = base.dPhiDgap * dgJ + base.dPhiDlambda * dLambda;
        Real dPhiJWithDr = dPhiJ;
        if (std::abs(base.complianceScale) > eps)
            dPhiJWithDr += (base.dPhiDlambda / base.complianceScale) * base.lambda * dr;

        if (std::abs(dPhiFD -dPhiJ) > 1.0)
        {
            msg_info() << "[NCP FD CONTACT] point=" << i
                    << " alpha=" << alpha
                    << " g0=" << base.gap
                    << " g1=" << trial.gap
                    << " lambda0=" << base.lambda
                    << " lambda1=" << trial.lambda
                    << " dLambda=" << dLambda
                    << " s0=" << base.complianceScale * base.lambda
                    << " s1=" << trial.complianceScale * trial.lambda
                    << " phi0=" << base.phi
                    << " phi1=" << trial.phi
                    << " dPhiDg=" << base.dPhiDgap
                    << " dPhiDl=" << base.dPhiDlambda
                    << " dPhiFD=" << dPhiFD
                    << " dPhiJ=" << dPhiJ
                    << " dr=" << dr;
        }

        msg_info() << "[NCP FD CONTACT] point=" << i
                   << " alpha=" << alpha
                   << " dgFD=" << dgFD
                   << " Hdx=" << dgJ
                   << " dgErr=" << (dgFD - dgJ)
                   << " dHFD=" << dHfd
                   << " HessDx=" << dHj
                   << " dHErr=" << (dHfd - dHj)
                   << " dFfd=" << dForceFD
                   << " dFJ=" << dForceJ
                   << " dFErr=" << (dForceFD - dForceJ)
                   << " dPhiFD=" << dPhiFD
                   << " dPhiJ=" << dPhiJ
                   << " dPhiJ+dr=" << dPhiJWithDr
                   << " dr=" << dr
                   << " Hessian=" << hasHessian;
    }
}

template<class T1, class T2>
void FischerBurmeisterContactForceField<T1, T2>::publishDebugData()
{
    const ContactDiagnostics diagnostics = summarizeContacts();
    d_activeContactCount.setValue(diagnostics.activeCount);
    d_pinnedContactCount.setValue(diagnostics.pinnedCount);
    d_invalidContactCount.setValue(diagnostics.invalidCount);

    if (!d_publishContactData.getValue())
        return;

    const sofa::Size n = m_contacts.size();
    sofa::type::vector<unsigned int> status(n);
    sofa::type::vector<Vec3> gradient(n);
    sofa::type::vector<Real> gap(n), lambda(n), r(n), scaledLambda(n), phi(n), dPhiDgap(n), beta(n);

    for (sofa::Index i = 0; i < n; ++i)
    {
        const Contact& c = m_contacts[i];
        status[i] = static_cast<unsigned int>(c.status);
        gradient[i] = c.gapGradient;
        gap[i] = c.gap;
        lambda[i] = c.lambda;
        r[i] = c.complianceScale;
        scaledLambda[i] = c.scaledLambda;
        phi[i] = c.phi;
        dPhiDgap[i] = c.dPhiDgap;
        beta[i] = c.dPhiDlambda;
    }

    d_contactStatus.setValue(status);
    d_contactGapGradient.setValue(gradient);
    d_contactGap.setValue(gap);
    d_contactLambda.setValue(lambda);
    d_contactComplianceScale.setValue(r);
    d_contactScaledLambda.setValue(scaledLambda);
    d_contactPhi.setValue(phi);
    d_contactDPhiDgap.setValue(dPhiDgap);
    d_contactDPhiDlambda.setValue(beta);
}

template<class T1, class T2>
void FischerBurmeisterContactForceField<T1, T2>::addForce(const sofa::core::MechanicalParams*,DataVecDeriv1& dataF1,DataVecDeriv2& dataF2,const DataVecCoord1& dataX1,const DataVecCoord2& dataX2,const DataVecDeriv1&,const DataVecDeriv2&)
{
    if (!rebuildCurrentContacts(dataX1.getValue(), dataX2.getValue()))
        return;

    VecDeriv1& f1 = *dataF1.beginEdit();
    VecDeriv2& f2 = *dataF2.beginEdit();

    for (const Contact& c : m_contacts)
    {
        if (c.status == ContactStatus::Active)
        {
            // Physical contact residual: R_x^c = H^T lambda.
            for (sofa::Size d = 0; d < TranslationalDim; ++d)
                f1[c.pointIndex][d] += c.lambda * c.gapGradient[d];
        }

        f2[c.lambdaIndex][0] += c.phi;
    }

    dataF1.endEdit();
    dataF2.endEdit();
}

template<class T1, class T2>
void FischerBurmeisterContactForceField<T1, T2>::addDForce(const sofa::core::MechanicalParams* mparams,DataVecDeriv1& dataDF1,DataVecDeriv2& dataDF2,const DataVecDeriv1& dataDX1,const DataVecDeriv2& dataDX2)
{
    if (!m_validState || !this->mstate1)
        return;

    const VecDeriv1& dx1 = dataDX1.getValue();
    const VecDeriv2& dx2 = dataDX2.getValue();
    VecDeriv1& df1 = *dataDF1.beginEdit();
    VecDeriv2& df2 = *dataDF2.beginEdit();
    const Real k = mparams->kFactor();

    // Needed only by the optional J_xx = lambda Hess(g) contribution.
    const auto x1Data = this->mstate1->read(core::vec_id::read_access::position);
    const VecCoord1& x1 = x1Data->getValue();

    for (const Contact& c : m_contacts)
    {
        const Real deltaLambda = dx2[c.lambdaIndex][0];

        Real deltaGap = Real(0);
        for (sofa::Size d = 0; d < TranslationalDim; ++d)
            deltaGap += c.gapGradient[d] * dx1[c.pointIndex][d];

        if (c.status == ContactStatus::Active)
        {
            // d(H^T lambda)/d(lambda) * deltaLambda = H^T deltaLambda.
            for (sofa::Size d = 0; d < TranslationalDim; ++d)
                df1[c.pointIndex][d] += k * c.gapGradient[d] * deltaLambda;

            // Optional d(H^T lambda)/dx * deltaX
            //          = lambda Hess(g) deltaX.
            Mat3 gapHessian;
            if (c.lambda != Real(0)
                && c.pointIndex < x1.size()
                && computeGapHessian(extractPosition(x1[c.pointIndex]), gapHessian))
            {
                for (sofa::Size row = 0; row < TranslationalDim; ++row)
                {
                    Real value = Real(0);
                    for (sofa::Size col = 0; col < TranslationalDim; ++col)
                        value += gapHessian(row, col) * dx1[c.pointIndex][col];

                    df1[c.pointIndex][row] += k * c.lambda * value;
                }
            }
        }

        // d phi = a H deltaX + b deltaLambda.
        df2[c.lambdaIndex][0] += k * (
            c.dPhiDgap * deltaGap
            + c.dPhiDlambda * deltaLambda);
    }

    dataDF1.endEdit();
    dataDF2.endEdit();
}

template<class T1, class T2>
void FischerBurmeisterContactForceField<T1, T2>::addKToMatrix(const sofa::core::MechanicalParams*,const sofa::core::behavior::MultiMatrixAccessor*)
{
}

template<class T1, class T2>
void FischerBurmeisterContactForceField<T1, T2>::buildStiffnessMatrix(core::behavior::StiffnessMatrix* matrix)
{
    if (!matrix || !m_validState || !this->mstate1 || !this->mstate2)
        return;

    // Exact contact Jacobian for
    //
    //     R_x      = H^T lambda
    //     R_lambda = phi(g, r lambda)
    //
    // with r frozen during this linearization:
    //
    //         [ lambda Hess(g)    H^T ]
    //     J = [                       ]
    //         [      a H          b   ]
    //
    // The Hessian hook returns false in the base class, so J_xx is zero until
    // a geometry specialization provides Hess(g).

    auto dRx_dX = matrix->getForceDerivativeIn(this->mstate1.get())
        .withRespectToPositionsIn(this->mstate1.get());
    auto dRx_dLambda = matrix->getForceDerivativeIn(this->mstate1.get())
        .withRespectToPositionsIn(this->mstate2.get());
    auto dPhi_dX = matrix->getForceDerivativeIn(this->mstate2.get())
        .withRespectToPositionsIn(this->mstate1.get());
    auto dPhi_dLambda = matrix->getForceDerivativeIn(this->mstate2.get())
        .withRespectToPositionsIn(this->mstate2.get());

    dRx_dX.checkValidity(this);
    dRx_dLambda.checkValidity(this);
    dPhi_dX.checkValidity(this);
    dPhi_dLambda.checkValidity(this);

    sofa::type::Mat<DerivDim1, DerivDim1, Real> upperLeft;
    sofa::type::Mat<DerivDim1, DerivDim2, Real> upperRight;
    sofa::type::Mat<DerivDim2, DerivDim1, Real> lowerLeft;
    sofa::type::Mat<DerivDim2, DerivDim2, Real> lowerRight;

    const auto x1Data = this->mstate1->read(core::vec_id::read_access::position);
    const VecCoord1& x1 = x1Data->getValue();

    for (const Contact& c : m_contacts)
    {
        upperLeft.clear();
        upperRight.clear();
        lowerLeft.clear();
        lowerRight.clear();

        if (c.status == ContactStatus::Active)
        {
            // J_xlambda = d(H^T lambda)/d(lambda) = H^T.
            for (sofa::Size d = 0; d < TranslationalDim; ++d)
                upperRight(d, 0) = c.gapGradient[d];

            // J_lambdax = dphi/dg * dg/dx = a H.
            for (sofa::Size d = 0; d < TranslationalDim; ++d)
                lowerLeft(0, d) = c.dPhiDgap * c.gapGradient[d];

            // Optional J_xx = lambda Hess(g).
            Mat3 gapHessian;
            if (c.lambda != Real(0)
                && c.pointIndex < x1.size()
                && computeGapHessian(extractPosition(x1[c.pointIndex]), gapHessian))
            {
                for (sofa::Size row = 0; row < TranslationalDim; ++row)
                    for (sofa::Size col = 0; col < TranslationalDim; ++col)
                        upperLeft(row, col) = c.lambda * gapHessian(row, col);
            }
        }

        // Active row: dphi/dlambda = b.
        // Pinned/invalid row: finalizeContactRow() sets b=1 and H=0, yielding the simple equation lambda=0.
        lowerRight(0, 0) = c.dPhiDlambda + 1e-6;

        dRx_dX(
            DerivDim1 * c.pointIndex,
            DerivDim1 * c.pointIndex) += upperLeft;

        dRx_dLambda(
            DerivDim1 * c.pointIndex,
            DerivDim2 * c.lambdaIndex) += upperRight;

        dPhi_dX(
            DerivDim2 * c.lambdaIndex,
            DerivDim1 * c.pointIndex) += lowerLeft;

        dPhi_dLambda(
            DerivDim2 * c.lambdaIndex,
            DerivDim2 * c.lambdaIndex) += lowerRight;
    }
}

template<class T1, class T2>
void FischerBurmeisterContactForceField<T1, T2>::buildDampingMatrix(core::behavior::DampingMatrix*)
{
}

template<class T1, class T2>
SReal FischerBurmeisterContactForceField<T1, T2>::getPotentialEnergy(const sofa::core::MechanicalParams*,const DataVecCoord1&,const DataVecCoord2&) const
{
    return SReal(0);
}

template<class T1, class T2>
void FischerBurmeisterContactForceField<T1, T2>::draw(const core::visual::VisualParams* vparams)
{
    drawActiveNormals(vparams);
}

template<class T1, class T2>
void FischerBurmeisterContactForceField<T1, T2>::drawActiveNormals(const core::visual::VisualParams* vparams) const
{
    if (!vparams || !vparams->drawTool() || !d_showContactGradients.getValue() || !this->mstate1)
        return;

    const auto x1Data = this->mstate1->read(core::vec_id::read_access::position);
    const VecCoord1& x1 = x1Data->getValue();
    std::vector<sofa::type::Vec3> lines;
    lines.reserve(2 * m_contacts.size());
    const Real scale = d_drawGradientScale.getValue();

    for (const Contact& c : m_contacts)
    {
        const Real norm2 = c.gapGradient.norm2();
        if (c.status != ContactStatus::Active
            || c.pointIndex >= x1.size()
            || c.gap > Real(1e-2)
            || norm2 <= Real(1e-30))
        {
            continue;
        }

        const Vec3 p0 = extractPosition(x1[c.pointIndex]);
        const Vec3 p1 = p0 + c.gapGradient * (scale / std::sqrt(norm2));
        lines.emplace_back(static_cast<float>(p0[0]), static_cast<float>(p0[1]), static_cast<float>(p0[2]));
        lines.emplace_back(static_cast<float>(p1[0]), static_cast<float>(p1[1]), static_cast<float>(p1[2]));
    }

    if (!lines.empty())
        vparams->drawTool()->drawLines(lines, 2.0f, d_contactColor.getValue());
}

} // namespace sofa::ncp