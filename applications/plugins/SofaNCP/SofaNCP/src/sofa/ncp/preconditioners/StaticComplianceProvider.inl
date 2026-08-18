#pragma once

#include <sofa/ncp/preconditioners/StaticComplianceProvider.h>

#include <sofa/core/behavior/LinearSolver.h>
#include <sofa/helper/ScopedAdvancedTimer.h>
#include <sofa/helper/logging/Messaging.h>

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace sofa::ncp
{

template<class TDataTypes>
StaticComplianceProvider<TDataTypes>::StaticComplianceProvider()
    : Inherit(nullptr)
    , l_mechanicalState(initLink("mechanicalState", "Mechanical state whose point compliance is extracted."))
    , l_linearSystem(initLink("linearSystem", "Optional explicit linear-system source."))
    , d_debug(this->initData(&d_debug, false, "debug", "Log successful compliance captures."))
{
}

template<class TDataTypes>
void StaticComplianceProvider<TDataTypes>::init()
{
    Inherit::init();

    if (!l_mechanicalState)
        l_mechanicalState.set(this->getContext()->template get<MechanicalState>());

    if (!l_mechanicalState)
    {
        msg_error() << "StaticComplianceProvider requires a mechanicalState link.";
        this->d_componentState.setValue(core::objectmodel::ComponentState::Invalid);
        return;
    }

    this->d_componentState.setValue(core::objectmodel::ComponentState::Valid);
}

template<class TDataTypes>
void StaticComplianceProvider<TDataTypes>::reinit()
{
    msg_warning() << "[COMPLIANCE] reinit -> INVALIDATE"
                  << " generation=" << m_generation
                  << " blocks=" << m_pointCompliance.size();

    Inherit::reinit();
    invalidate();
}

template<class TDataTypes>
core::behavior::BaseMatrixLinearSystem*
StaticComplianceProvider<TDataTypes>::resolveLinearSystem() const
{
    if (auto* system = l_linearSystem.get())
        return system;

    if (auto* solver = this->l_linearSolver.get())
        return solver->getLinearSystem();

    return nullptr;
}

template<class TDataTypes>
void StaticComplianceProvider<TDataTypes>::invalidate()
{
    msg_warning() << "[COMPLIANCE] invalidate"
                  << " generation=" << m_generation
                  << " blocks=" << m_pointCompliance.size();

    m_valid = false;
    m_pointCompliance.clear();
}

template<class TDataTypes>
bool StaticComplianceProvider<TDataTypes>::captureFromCurrentLinearSystem()
{
    SCOPED_TIMER("StaticComplianceCapture");

    auto* system = resolveLinearSystem();
    const auto* matrix = system ? system->getSystemBaseMatrix() : nullptr;

    if (!matrix)
    {
        msg_warning() << "Static compliance capture requires an assembled explicit linear system; previous data retained.";
        return false;
    }

    return extractPointCompliance(*matrix);
}

template<class TDataTypes>
bool StaticComplianceProvider<TDataTypes>::extractPointCompliance(const sofa::linearalgebra::BaseMatrix& matrix)
{
    const auto* state = l_mechanicalState.get();
    const sofa::Size pointCount = state ? state->getSize() : 0;
    const sofa::Size objectDofs = pointCount * DerivDim;
    using MatrixIndex = sofa::linearalgebra::BaseMatrix::Index;

    if (objectDofs == 0 || objectDofs > static_cast<sofa::Size>(std::numeric_limits<MatrixIndex>::max()))
    {
        msg_warning() << "Static compliance capture has an invalid mechanical block size; previous data retained.";
        return false;
    }

    const MatrixIndex requiredDofs = static_cast<MatrixIndex>(objectDofs);
    if (matrix.rowSize() < requiredDofs || matrix.colSize() < requiredDofs)
    {
        msg_warning() << "Static compliance matrix is " << matrix.rowSize() << "x"
                      << matrix.colSize() << ", but the leading mechanical block requires "
                      << objectDofs << "x" << objectDofs << "; previous data retained.";
        return false;
    }

    using DenseMatrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
    using DenseVector = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

    DenseMatrix stiffness(
        static_cast<Eigen::Index>(objectDofs),
        static_cast<Eigen::Index>(objectDofs));

    Real maximumMagnitude = Real(0);
    Real maximumAsymmetry = Real(0);

    for (sofa::Size row = 0; row < objectDofs; ++row)
    {
        for (sofa::Size column = 0; column < objectDofs; ++column)
        {
            const Real value = static_cast<Real>(matrix.element(row, column));
            stiffness(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) = value;
            maximumMagnitude = std::max(maximumMagnitude, std::abs(value));
        }
    }

    for (sofa::Size row = 0; row < objectDofs; ++row)
    {
        for (sofa::Size column = row + 1; column < objectDofs; ++column)
        {
            maximumAsymmetry = std::max(
                maximumAsymmetry,
                std::abs(
                    stiffness(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column))
                    - stiffness(static_cast<Eigen::Index>(column), static_cast<Eigen::Index>(row))));
        }
    }

    const Real tolerance = Real(1e-8) * std::max(maximumMagnitude, Real(1));
    if (maximumAsymmetry > tolerance)
    {
        msg_warning() << "Static compliance leading mechanical block is not symmetric: max asymmetry="
                      << maximumAsymmetry << ", tolerance=" << tolerance
                      << "; previous data retained.";
        return false;
    }

    Eigen::LDLT<DenseMatrix> factorization(stiffness);
    if (factorization.info() != Eigen::Success)
    {
        msg_warning() << "Static compliance LDLT factorization failed; previous data retained.";
        return false;
    }

    sofa::type::vector<Mat3> candidate(pointCount, Mat3());
    DenseVector rhs(static_cast<Eigen::Index>(objectDofs));

    for (sofa::Index point = 0; point < pointCount; ++point)
    {
        for (sofa::Size column = 0; column < TranslationalDim; ++column)
        {
            rhs.setZero();
            rhs(static_cast<Eigen::Index>(localDof(point, column))) = Real(1);
            const DenseVector displacement = factorization.solve(rhs);

            if (factorization.info() != Eigen::Success || !displacement.allFinite())
            {
                msg_warning() << "Static compliance solve failed at point=" << point
                              << ", component=" << column
                              << "; previous data retained.";
                return false;
            }

            for (sofa::Size row = 0; row < TranslationalDim; ++row)
            {
                candidate[point](row, column) = displacement(
                    static_cast<Eigen::Index>(localDof(point, row)));
            }
        }
    }

    m_pointCompliance = std::move(candidate);
    m_valid = true;
    ++m_generation;

    if (d_debug.getValue())
        msg_info() << "Captured " << pointCount << " point-compliance blocks, generation="
                   << m_generation << ".";

    return true;
}

template<class TDataTypes>
bool StaticComplianceProvider<TDataTypes>::hasCompliance(
    sofa::Size expectedPointCount) const
{
    const bool available =
        m_valid && m_pointCompliance.size() == expectedPointCount;

    if (!available)
    {
        msg_warning()
            << "[COMPLIANCE UNAVAILABLE]"
            << " valid=" << m_valid
            << " generation=" << m_generation
            << " capturedPoints=" << m_pointCompliance.size()
            << " expectedPoints=" << expectedPointCount;
    }

    return available;
}

template<class TDataTypes>
typename StaticComplianceProvider<TDataTypes>::Real
StaticComplianceProvider<TDataTypes>::scalarCompliance(sofa::Index pointIndex,const Vec3& gapGradient) const
{
    if (!m_valid || pointIndex >= m_pointCompliance.size())
        return Real(0);

    const Mat3& compliance = m_pointCompliance[pointIndex];
    Real value = Real(0);

    for (sofa::Size row = 0; row < TranslationalDim; ++row)
    {
        for (sofa::Size column = 0; column < TranslationalDim; ++column)
        {
            value += gapGradient[row] * compliance(row, column) * gapGradient[column];
        }
    }

    return value;
}

} // namespace sofa::ncp
