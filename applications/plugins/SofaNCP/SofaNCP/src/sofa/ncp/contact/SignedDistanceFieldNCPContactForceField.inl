/****************************************************************************
* SDF Fischer-Burmeister NCP contact forcefield implementation.
****************************************************************************/
#pragma once

#include <sofa/ncp/contact/SignedDistanceFieldNCPContactForceField.h>
#include <sofa/ncp/contact/FischerBurmeisterContactForceField.inl>
#include <sofa/core/visual/VisualParams.h>

#include <cmath>
#include <fstream>
#include <vector>

namespace sofa::ncp
{

template<class TDataTypes1, class TDataTypes2>
SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::SignedDistanceFieldNCPContactForceField()
    : Base()
    , d_sdfFilename(this->initData(&d_sdfFilename, std::string(), "sdfFilename", "Path to float32 raw phi file. Layout: phi[ix,iy,iz], z fastest."))
    , d_sdfGradientFilename(this->initData(&d_sdfGradientFilename, std::string(), "sdfGradientFilename", "Optional float32 first derivatives. Layout: [fx,fy,fz] per node."))
    , d_sdfHermiteMixedFilename(this->initData(&d_sdfHermiteMixedFilename, std::string(), "sdfHermiteMixedFilename", "Optional float32 mixed derivatives. Layout: [fxy,fxz,fyz,fxyz] per node."))
    , d_sdfHermitePacketFilename(this->initData(&d_sdfHermitePacketFilename, std::string(), "sdfHermitePacketFilename", "Optional float32 full Hermite packet. Layout: [fx,fy,fz,fxy,fxz,fyz,fxyz] per node. Overrides derivative files."))
    , d_origin(this->initData(&d_origin, Vec3(Real(0), Real(0), Real(0)), "sdfOrigin", "World position of grid index 0 0 0."))
    , d_spacing(this->initData(&d_spacing, Vec3(Real(1), Real(1), Real(1)), "sdfSpacing", "Grid spacing dx dy dz."))
    , d_dimensions(this->initData(&d_dimensions, UInt3(0, 0, 0), "sdfDimensions", "Grid dimensions nx ny nz."))
    , d_interpolationMode(this->initData(&d_interpolationMode, static_cast<unsigned int>(Cubic64), "interpolationMode", "0=Auto, 1=Linear8, 2=Cubic64, 3=HermiteFirstDerivatives, 4=HermiteFull."))
    , d_normalizeGradient(this->initData(&d_normalizeGradient, false, "normalizeGradient", "Normalize grad(phi) before returning gapGradient. This changes the Jacobian unless phi is an exact SDF."))
    , d_debugQueries(this->initData(&d_debugQueries, true, "debugQueries", "Log invalid, out-of-grid, boundary-fallback, and interpolation queries."))
    , d_showGridBox(this->initData(&d_showGridBox, false, "showGridBox", "Draw SDF grid bounding box."))
    , d_color(this->initData(&d_color, sofa::type::RGBAColor(0.1f, 0.8f, 1.0f, 1.0f), "color", "SDF grid box draw color."))
{
}

template<class TDataTypes1, class TDataTypes2>
SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::SignedDistanceFieldNCPContactForceField(
    core::behavior::MechanicalState<DataTypes1>* object1,
    core::behavior::MechanicalState<DataTypes2>* object2)
    : SignedDistanceFieldNCPContactForceField()
{
    this->mstate1.set(object1);
    this->mstate2.set(object2);
}

template<class TDataTypes1, class TDataTypes2>
void SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::init()
{
    Base::init();
    loadGridData();
}

template<class TDataTypes1, class TDataTypes2>
void SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::reinit()
{
    Base::reinit();
    loadGridData();
}

template<class TDataTypes1, class TDataTypes2>
const char* SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::interpolationModeName(unsigned int mode)
{
    switch (mode)
    {
    case Auto: return "Auto";
    case Linear8: return "Linear8";
    case Cubic64: return "Cubic64";
    case HermiteFirstDerivatives: return "HermiteFirstDerivatives";
    case HermiteFull: return "HermiteFull";
    default: return "Unknown";
    }
}

template<class TDataTypes1, class TDataTypes2>
void SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::clearGridData()
{
    m_phi.clear();
    m_gradientData.clear();
    m_mixedDerivativeData.clear();

    m_nx = 0;
    m_ny = 0;
    m_nz = 0;

    m_origin = Vec3(Real(0), Real(0), Real(0));
    m_spacing = Vec3(Real(1), Real(1), Real(1));

    m_loaded = false;
    m_hasHermiteFirstDerivatives = false;
    m_hasHermiteMixedDerivatives = false;

    m_warnedMissingFirstDerivatives = false;
    m_warnedMissingFullHermite = false;
    m_warnedCubicBoundaryFallback = false;
}

template<class TDataTypes1, class TDataTypes2>
bool SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::loadGridData()
{
    clearGridData();

    const std::string phiFilename = d_sdfFilename.getValue();
    const std::string gradFilename = d_sdfGradientFilename.getValue();
    const std::string mixedFilename = d_sdfHermiteMixedFilename.getValue();
    const std::string packetFilename = d_sdfHermitePacketFilename.getValue();
    const UInt3 dims = d_dimensions.getValue();
    const Vec3 spacing = d_spacing.getValue();
    const Vec3 origin = d_origin.getValue();

    if (phiFilename.empty())
    {
        msg_error() << "Missing sdfFilename.";
        return false;
    }

    if (dims[0] < 2 || dims[1] < 2 || dims[2] < 2)
    {
        msg_error() << "Invalid SDF dimensions: " << dims << ". Expected nx ny nz >= 2.";
        return false;
    }

    if (spacing[0] <= Real(0) || spacing[1] <= Real(0) || spacing[2] <= Real(0))
    {
        msg_error() << "Invalid SDF spacing: " << spacing << ". Expected strictly positive spacing.";
        return false;
    }

    m_nx = dims[0];
    m_ny = dims[1];
    m_nz = dims[2];
    m_origin = origin;
    m_spacing = spacing;

    const sofa::Size count =
        static_cast<sofa::Size>(m_nx) *
        static_cast<sofa::Size>(m_ny) *
        static_cast<sofa::Size>(m_nz);

    auto readFloatFile = [](const std::string& filename, sofa::Size expectedCount, std::vector<float>& dst) -> bool
    {
        dst.assign(static_cast<std::size_t>(expectedCount), 0.0f);
        std::ifstream file(filename.c_str(), std::ios::binary);
        if (!file)
            return false;

        file.read(reinterpret_cast<char*>(dst.data()), static_cast<std::streamsize>(expectedCount * sizeof(float)));
        return static_cast<bool>(file);
    };

    if (!readFloatFile(phiFilename, count, m_phi))
    {
        msg_error() << "Failed reading SDF phi raw file: " << phiFilename << ". Expected " << count << " float32 values.";
        clearGridData();
        return false;
    }

    if (!packetFilename.empty())
    {
        std::vector<float> packetTmp;
        if (!readFloatFile(packetFilename, count * 7, packetTmp))
        {
            msg_error() << "Failed reading full Hermite packet raw file: " << packetFilename
                        << ". Expected " << count * 7 << " float32 values.";
            clearGridData();
            return false;
        }

        m_gradientData.resize(static_cast<std::size_t>(count * 3));
        m_mixedDerivativeData.resize(static_cast<std::size_t>(count * 4));
        for (sofa::Size i = 0; i < count; ++i)
        {
            const std::size_t packetOffset = static_cast<std::size_t>(7 * i);
            const std::size_t gradientOffset = static_cast<std::size_t>(3 * i);
            const std::size_t mixedOffset = static_cast<std::size_t>(4 * i);
            for (std::size_t component = 0; component < 3; ++component) m_gradientData[gradientOffset + component] = packetTmp[packetOffset + component];
            for (std::size_t component = 0; component < 4; ++component) m_mixedDerivativeData[mixedOffset + component] = packetTmp[packetOffset + 3 + component];
        }

        m_hasHermiteFirstDerivatives = true;
        m_hasHermiteMixedDerivatives = true;
    }
    else
    {
        if (!gradFilename.empty())
        {
            std::vector<float> gradTmp;
            if (!readFloatFile(gradFilename, count * 3, gradTmp))
            {
                msg_error() << "Failed reading SDF gradient raw file: " << gradFilename
                            << ". Expected " << count * 3 << " float32 values.";
                clearGridData();
                return false;
            }

            m_gradientData = std::move(gradTmp);
            m_hasHermiteFirstDerivatives = true;
        }

        if (!mixedFilename.empty())
        {
            std::vector<float> mixedTmp;
            if (!readFloatFile(mixedFilename, count * 4, mixedTmp))
            {
                msg_error() << "Failed reading Hermite mixed derivative raw file: " << mixedFilename
                            << ". Expected " << count * 4 << " float32 values.";
                clearGridData();
                return false;
            }

            m_mixedDerivativeData = std::move(mixedTmp);
            m_hasHermiteMixedDerivatives = true;
        }
    }

    if (m_hasHermiteMixedDerivatives && !m_hasHermiteFirstDerivatives)
    {
        msg_error() << "Hermite mixed derivatives were loaded without first derivatives. "
                    << "Provide sdfGradientFilename or use sdfHermitePacketFilename.";
        clearGridData();
        return false;
    }

    m_loaded = true;

    msg_info() << "Loaded SDF phi='" << phiFilename
               << "' dims=" << m_nx << " " << m_ny << " " << m_nz
               << " origin=" << m_origin
               << " spacing=" << m_spacing
               << " interpolationMode=" << interpolationModeName(d_interpolationMode.getValue())
               << " hasFirstDerivatives=" << m_hasHermiteFirstDerivatives
               << " hasMixedDerivatives=" << m_hasHermiteMixedDerivatives
               << " normalizeGradient=" << d_normalizeGradient.getValue();

    return true;
}

template<class TDataTypes1, class TDataTypes2>
sofa::Size SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::flattenedIndex(
    unsigned int i,
    unsigned int j,
    unsigned int k) const
{
    return (static_cast<sofa::Size>(i) * static_cast<sofa::Size>(m_ny)
          + static_cast<sofa::Size>(j)) * static_cast<sofa::Size>(m_nz)
          + static_cast<sofa::Size>(k);
}

template<class TDataTypes1, class TDataTypes2>
typename SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::Real
SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::phiAt(unsigned int i,unsigned int j,unsigned int k) const
{
    return static_cast<Real>(m_phi[flattenedIndex(i, j, k)]);
}

template<class TDataTypes1, class TDataTypes2>
typename SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::Vec3
SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::gradAt(unsigned int i,unsigned int j,unsigned int k) const
{
    const std::size_t offset = static_cast<std::size_t>(3 * flattenedIndex(i, j, k));
    return Vec3(static_cast<Real>(m_gradientData[offset]), static_cast<Real>(m_gradientData[offset + 1]), static_cast<Real>(m_gradientData[offset + 2]));
}

template<class TDataTypes1, class TDataTypes2>
typename SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::HermiteMixedDerivatives
SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::mixedAt(unsigned int i,unsigned int j,unsigned int k) const
{
    const std::size_t offset = static_cast<std::size_t>(4 * flattenedIndex(i, j, k));
    HermiteMixedDerivatives derivatives;
    derivatives.fxy = static_cast<Real>(m_mixedDerivativeData[offset]);
    derivatives.fxz = static_cast<Real>(m_mixedDerivativeData[offset + 1]);
    derivatives.fyz = static_cast<Real>(m_mixedDerivativeData[offset + 2]);
    derivatives.fxyz = static_cast<Real>(m_mixedDerivativeData[offset + 3]);
    return derivatives;
}

template<class TDataTypes1, class TDataTypes2>
bool SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::validateAndNormalizeGradient(const Vec3& p,Real phi,Vec3& gradPhi) const
{
    if (!std::isfinite(phi) ||
        !std::isfinite(gradPhi[0]) ||
        !std::isfinite(gradPhi[1]) ||
        !std::isfinite(gradPhi[2]))
    {
        if (d_debugQueries.getValue()) msg_warning() << "[SDF invalid] non-finite interpolation p=" << p
                      << " phi=" << phi << " grad=" << gradPhi;
        return false;
    }

    if (!d_normalizeGradient.getValue())
        return true;

    const Real n = gradPhi.norm();
    if (!std::isfinite(n) || n < Real(1e-8))
    {
        if (d_debugQueries.getValue()) msg_info() << "[SDF invalid] near-zero gradient p=" << p
                   << " phi=" << phi << " grad=" << gradPhi << " norm=" << n;
        return false;
    }

    gradPhi /= n;
    return true;
}

template<class TDataTypes1, class TDataTypes2>
bool SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::evaluate(const Vec3& p,Real& phi,Vec3& gradPhi) const
{
    if (!m_loaded || m_phi.empty())
    {
        if (d_debugQueries.getValue()) msg_warning() << "[SDF invalid] not loaded or empty.";
        return false;
    }

    const unsigned int mode = d_interpolationMode.getValue();

    if (mode == Auto)
    {
        if (m_hasHermiteFirstDerivatives && m_hasHermiteMixedDerivatives)
            return evaluateHermiteFull(p, phi, gradPhi);
        return evaluateCubic64(p, phi, gradPhi);
    }

    if (mode == Linear8)
        return evaluateLinear8(p, phi, gradPhi);

    if (mode == Cubic64)
        return evaluateCubic64(p, phi, gradPhi);

    if (mode == HermiteFirstDerivatives)
    {
        if (!m_hasHermiteFirstDerivatives)
        {
            if (!m_warnedMissingFirstDerivatives)
            {
                msg_warning() << "interpolationMode=HermiteFirstDerivatives but no sdfGradientFilename/sdfHermitePacketFilename was loaded. Falling back to Cubic64.";
                m_warnedMissingFirstDerivatives = true;
            }
            return evaluateCubic64(p, phi, gradPhi);
        }
        return evaluateHermiteFirst(p, phi, gradPhi);
    }

    if (mode == HermiteFull)
    {
        if (!m_hasHermiteFirstDerivatives || !m_hasHermiteMixedDerivatives)
        {
            if (!m_warnedMissingFullHermite)
            {
                msg_warning() << "interpolationMode=HermiteFull requires first and mixed derivatives. Falling back to Cubic64.";
                m_warnedMissingFullHermite = true;
            }
            return evaluateCubic64(p, phi, gradPhi);
        }
        return evaluateHermiteFull(p, phi, gradPhi);
    }

    msg_warning() << "Unknown interpolationMode=" << mode << ". Falling back to Cubic64.";
    return evaluateCubic64(p, phi, gradPhi);
}

template<class TDataTypes1, class TDataTypes2>
bool SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::evaluateLinear8(const Vec3& p,Real& phi,Vec3& gradPhi) const
{
    const Real gx = (p[0] - m_origin[0]) / m_spacing[0];
    const Real gy = (p[1] - m_origin[1]) / m_spacing[1];
    const Real gz = (p[2] - m_origin[2]) / m_spacing[2];

    if (gx < Real(0) || gx >= Real(m_nx - 1) ||
        gy < Real(0) || gy >= Real(m_ny - 1) ||
        gz < Real(0) || gz >= Real(m_nz - 1))
    {
        if (d_debugQueries.getValue()) msg_info() << "[SDF invalid] out of grid p=" << p
                   << " grid=(" << gx << ", " << gy << ", " << gz << ")"
                   << " dims=(" << m_nx << ", " << m_ny << ", " << m_nz << ")";
        return false;
    }

    const unsigned int i = static_cast<unsigned int>(std::floor(gx));
    const unsigned int j = static_cast<unsigned int>(std::floor(gy));
    const unsigned int k = static_cast<unsigned int>(std::floor(gz));

    const Real u = gx - Real(i);
    const Real v = gy - Real(j);
    const Real w = gz - Real(k);

    const Real omu = Real(1) - u;
    const Real omv = Real(1) - v;
    const Real omw = Real(1) - w;

    const Real c000 = phiAt(i,     j,     k);
    const Real c100 = phiAt(i + 1, j,     k);
    const Real c010 = phiAt(i,     j + 1, k);
    const Real c110 = phiAt(i + 1, j + 1, k);
    const Real c001 = phiAt(i,     j,     k + 1);
    const Real c101 = phiAt(i + 1, j,     k + 1);
    const Real c011 = phiAt(i,     j + 1, k + 1);
    const Real c111 = phiAt(i + 1, j + 1, k + 1);

    phi =
        c000 * omu * omv * omw +
        c100 * u   * omv * omw +
        c010 * omu * v   * omw +
        c110 * u   * v   * omw +
        c001 * omu * omv * w   +
        c101 * u   * omv * w   +
        c011 * omu * v   * w   +
        c111 * u   * v   * w;

    gradPhi = Vec3(
        ((c100 - c000) * omv * omw +
         (c110 - c010) * v   * omw +
         (c101 - c001) * omv * w   +
         (c111 - c011) * v   * w) / m_spacing[0],

        ((c010 - c000) * omu * omw +
         (c110 - c100) * u   * omw +
         (c011 - c001) * omu * w   +
         (c111 - c101) * u   * w) / m_spacing[1],

        ((c001 - c000) * omu * omv +
         (c101 - c100) * u   * omv +
         (c011 - c010) * omu * v   +
         (c111 - c110) * u   * v) / m_spacing[2]);

    return validateAndNormalizeGradient(p, phi, gradPhi);
}

template<class TDataTypes1, class TDataTypes2>
bool SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::evaluateCubic64(const Vec3& p,Real& phi,Vec3& gradPhi) const
{
    const Real gx = (p[0] - m_origin[0]) / m_spacing[0];
    const Real gy = (p[1] - m_origin[1]) / m_spacing[1];
    const Real gz = (p[2] - m_origin[2]) / m_spacing[2];

    if (gx < Real(0) || gx >= Real(m_nx - 1) ||
        gy < Real(0) || gy >= Real(m_ny - 1) ||
        gz < Real(0) || gz >= Real(m_nz - 1))
    {
        if (d_debugQueries.getValue()) msg_info() << "[SDF invalid] out of grid p=" << p
                   << " grid=(" << gx << ", " << gy << ", " << gz << ")"
                   << " dims=(" << m_nx << ", " << m_ny << ", " << m_nz << ")";
        return false;
    }

    if (m_nx < 4 || m_ny < 4 || m_nz < 4)
        return evaluateLinear8(p, phi, gradPhi);

    const unsigned int i = static_cast<unsigned int>(std::floor(gx));
    const unsigned int j = static_cast<unsigned int>(std::floor(gy));
    const unsigned int k = static_cast<unsigned int>(std::floor(gz));

    if (i < 1 || i + 2 >= m_nx ||
        j < 1 || j + 2 >= m_ny ||
        k < 1 || k + 2 >= m_nz)
    {
        if (!m_warnedCubicBoundaryFallback)
        {
            if (d_debugQueries.getValue()) msg_info() << "Cubic64 query reached one-cell grid boundary. Falling back to Linear8 near boundary.";
            m_warnedCubicBoundaryFallback = true;
        }
        return evaluateLinear8(p, phi, gradPhi);
    }

    const Real u = gx - Real(i);
    const Real v = gy - Real(j);
    const Real w = gz - Real(k);

    auto catmullRomBasis = [](Real t, Real W[4], Real dW[4])
    {
        const Real t2 = t * t;
        const Real t3 = t2 * t;

        W[0] = Real(-0.5) * t + t2 - Real(0.5) * t3;
        W[1] = Real(1) - Real(2.5) * t2 + Real(1.5) * t3;
        W[2] = Real(0.5) * t + Real(2) * t2 - Real(1.5) * t3;
        W[3] = Real(-0.5) * t2 + Real(0.5) * t3;

        dW[0] = Real(-0.5) + Real(2) * t - Real(1.5) * t2;
        dW[1] = Real(-5) * t + Real(4.5) * t2;
        dW[2] = Real(0.5) + Real(4) * t - Real(4.5) * t2;
        dW[3] = -t + Real(1.5) * t2;
    };

    Real wx[4], dwx[4];
    Real wy[4], dwy[4];
    Real wz[4], dwz[4];
    catmullRomBasis(u, wx, dwx);
    catmullRomBasis(v, wy, dwy);
    catmullRomBasis(w, wz, dwz);

    phi = Real(0);
    Real ddu = Real(0);
    Real ddv = Real(0);
    Real ddw = Real(0);

    for (unsigned int a = 0; a < 4; ++a)
    {
        const unsigned int ii = i + a - 1;
        for (unsigned int b = 0; b < 4; ++b)
        {
            const unsigned int jj = j + b - 1;
            for (unsigned int c = 0; c < 4; ++c)
            {
                const unsigned int kk = k + c - 1;
                const Real f = phiAt(ii, jj, kk);

                phi += f * wx[a]  * wy[b]  * wz[c];
                ddu += f * dwx[a] * wy[b]  * wz[c];
                ddv += f * wx[a]  * dwy[b] * wz[c];
                ddw += f * wx[a]  * wy[b]  * dwz[c];
            }
        }
    }

    gradPhi = Vec3(ddu / m_spacing[0], ddv / m_spacing[1], ddw / m_spacing[2]);
    return validateAndNormalizeGradient(p, phi, gradPhi);
}

template<class TDataTypes1, class TDataTypes2>
bool SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::evaluateHermiteFirst(const Vec3& p,Real& phi,Vec3& gradPhi) const
{
    return evaluateHermiteCell(p, false, phi, gradPhi);
}

template<class TDataTypes1, class TDataTypes2>
bool SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::evaluateHermiteFull(const Vec3& p,Real& phi,Vec3& gradPhi) const
{
    return evaluateHermiteCell(p, true, phi, gradPhi);
}

template<class TDataTypes1, class TDataTypes2>
bool SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::evaluateHermiteCell(const Vec3& p,bool useMixedDerivatives,Real& phi,Vec3& gradPhi) const
{
    const Real gx = (p[0] - m_origin[0]) / m_spacing[0];
    const Real gy = (p[1] - m_origin[1]) / m_spacing[1];
    const Real gz = (p[2] - m_origin[2]) / m_spacing[2];

    if (gx < Real(0) || gx >= Real(m_nx - 1) ||
        gy < Real(0) || gy >= Real(m_ny - 1) ||
        gz < Real(0) || gz >= Real(m_nz - 1))
    {
        if (d_debugQueries.getValue()) msg_info() << "[SDF invalid] out of grid p=" << p
                   << " grid=(" << gx << ", " << gy << ", " << gz << ")"
                   << " dims=(" << m_nx << ", " << m_ny << ", " << m_nz << ")";
        return false;
    }

    const unsigned int i0 = static_cast<unsigned int>(std::floor(gx));
    const unsigned int j0 = static_cast<unsigned int>(std::floor(gy));
    const unsigned int k0 = static_cast<unsigned int>(std::floor(gz));

    const Real u = gx - Real(i0);
    const Real v = gy - Real(j0);
    const Real w = gz - Real(k0);

    auto hermiteBasis = [](Real t, Real V[2], Real D[2], Real dV[2], Real dD[2])
    {
        const Real t2 = t * t;
        const Real t3 = t2 * t;

        V[0]  = Real(2) * t3 - Real(3) * t2 + Real(1);
        V[1]  = Real(-2) * t3 + Real(3) * t2;
        D[0]  = t3 - Real(2) * t2 + t;
        D[1]  = t3 - t2;

        dV[0] = Real(6) * t2 - Real(6) * t;
        dV[1] = Real(-6) * t2 + Real(6) * t;
        dD[0] = Real(3) * t2 - Real(4) * t + Real(1);
        dD[1] = Real(3) * t2 - Real(2) * t;
    };

    Real Vx[2], Dx[2], dVx[2], dDx[2];
    Real Vy[2], Dy[2], dVy[2], dDy[2];
    Real Vz[2], Dz[2], dVz[2], dDz[2];
    hermiteBasis(u, Vx, Dx, dVx, dDx);
    hermiteBasis(v, Vy, Dy, dVy, dDy);
    hermiteBasis(w, Vz, Dz, dVz, dDz);

    const Real hx = m_spacing[0];
    const Real hy = m_spacing[1];
    const Real hz = m_spacing[2];
    const Real hxhy = hx * hy;
    const Real hxhz = hx * hz;
    const Real hyhz = hy * hz;
    const Real hxhyhz = hx * hy * hz;

    phi = Real(0);
    Real ddu = Real(0);
    Real ddv = Real(0);
    Real ddw = Real(0);

    for (unsigned int a = 0; a < 2; ++a)
    {
        const unsigned int ii = i0 + a;
        for (unsigned int b = 0; b < 2; ++b)
        {
            const unsigned int jj = j0 + b;
            for (unsigned int c = 0; c < 2; ++c)
            {
                const unsigned int kk = k0 + c;

                const Real F = phiAt(ii, jj, kk);
                const Vec3 G = gradAt(ii, jj, kk);

                const Real F000 = F;
                const Real Fx   = G[0];
                const Real Fy   = G[1];
                const Real Fz   = G[2];

                Real Fxy = Real(0);
                Real Fxz = Real(0);
                Real Fyz = Real(0);
                Real Fxyz = Real(0);
                if (useMixedDerivatives)
                {
                    const HermiteMixedDerivatives M = mixedAt(ii, jj, kk);
                    Fxy = M.fxy;
                    Fxz = M.fxz;
                    Fyz = M.fyz;
                    Fxyz = M.fxyz;
                }

                phi +=
                    F000        * Vx[a] * Vy[b] * Vz[c] +
                    Fx * hx     * Dx[a] * Vy[b] * Vz[c] +
                    Fy * hy     * Vx[a] * Dy[b] * Vz[c] +
                    Fz * hz     * Vx[a] * Vy[b] * Dz[c] +
                    Fxy * hxhy  * Dx[a] * Dy[b] * Vz[c] +
                    Fxz * hxhz  * Dx[a] * Vy[b] * Dz[c] +
                    Fyz * hyhz  * Vx[a] * Dy[b] * Dz[c] +
                    Fxyz * hxhyhz * Dx[a] * Dy[b] * Dz[c];

                ddu +=
                    F000        * dVx[a] * Vy[b] * Vz[c] +
                    Fx * hx     * dDx[a] * Vy[b] * Vz[c] +
                    Fy * hy     * dVx[a] * Dy[b] * Vz[c] +
                    Fz * hz     * dVx[a] * Vy[b] * Dz[c] +
                    Fxy * hxhy  * dDx[a] * Dy[b] * Vz[c] +
                    Fxz * hxhz  * dDx[a] * Vy[b] * Dz[c] +
                    Fyz * hyhz  * dVx[a] * Dy[b] * Dz[c] +
                    Fxyz * hxhyhz * dDx[a] * Dy[b] * Dz[c];

                ddv +=
                    F000        * Vx[a] * dVy[b] * Vz[c] +
                    Fx * hx     * Dx[a] * dVy[b] * Vz[c] +
                    Fy * hy     * Vx[a] * dDy[b] * Vz[c] +
                    Fz * hz     * Vx[a] * dVy[b] * Dz[c] +
                    Fxy * hxhy  * Dx[a] * dDy[b] * Vz[c] +
                    Fxz * hxhz  * Dx[a] * dVy[b] * Dz[c] +
                    Fyz * hyhz  * Vx[a] * dDy[b] * Dz[c] +
                    Fxyz * hxhyhz * Dx[a] * dDy[b] * Dz[c];

                ddw +=
                    F000        * Vx[a] * Vy[b] * dVz[c] +
                    Fx * hx     * Dx[a] * Vy[b] * dVz[c] +
                    Fy * hy     * Vx[a] * Dy[b] * dVz[c] +
                    Fz * hz     * Vx[a] * Vy[b] * dDz[c] +
                    Fxy * hxhy  * Dx[a] * Dy[b] * dVz[c] +
                    Fxz * hxhz  * Dx[a] * Vy[b] * dDz[c] +
                    Fyz * hyhz  * Vx[a] * Dy[b] * dDz[c] +
                    Fxyz * hxhyhz * Dx[a] * Dy[b] * dDz[c];
            }
        }
    }

    gradPhi = Vec3(ddu / hx, ddv / hy, ddw / hz);
    return validateAndNormalizeGradient(p, phi, gradPhi);
}

template<class TDataTypes1, class TDataTypes2>
bool SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::computeGapHessian(const Vec3& p,Mat3& hessian) const
{
    hessian.clear();

    if (!m_loaded || m_phi.empty())
        return false;

    // When normalizeGradient=true, computeContactKinematics returns
    //     H = -grad(phi) / ||grad(phi)||,
    // while gap itself remains g=-phi. H is therefore no longer dg/dx in
    // general, so -Hess(phi) would not be the derivative of the force
    // direction used by the base class. Do not inject an inconsistent J_xx.
    if (d_normalizeGradient.getValue())
        return false;

    Mat3 hessianPhi;
    hessianPhi.clear();

    const unsigned int mode = d_interpolationMode.getValue();
    bool valid = false;

    if (mode == Auto)
    {
        valid = m_hasHermiteFirstDerivatives && m_hasHermiteMixedDerivatives
            ? evaluateHermiteCellHessian(p, true, hessianPhi)
            : evaluateCubic64Hessian(p, hessianPhi);
    }
    else if (mode == Linear8)
    {
        valid = evaluateLinear8Hessian(p, hessianPhi);
    }
    else if (mode == Cubic64)
    {
        valid = evaluateCubic64Hessian(p, hessianPhi);
    }
    else if (mode == HermiteFirstDerivatives)
    {
        valid = m_hasHermiteFirstDerivatives
            ? evaluateHermiteCellHessian(p, false, hessianPhi)
            : evaluateCubic64Hessian(p, hessianPhi);
    }
    else if (mode == HermiteFull)
    {
        valid = m_hasHermiteFirstDerivatives && m_hasHermiteMixedDerivatives
            ? evaluateHermiteCellHessian(p, true, hessianPhi)
            : evaluateCubic64Hessian(p, hessianPhi);
    }
    else
    {
        valid = evaluateCubic64Hessian(p, hessianPhi);
    }

    if (!valid)
        return false;

    // NCP convention: g = -phi, hence Hess(g) = -Hess(phi).
    for (sofa::Size row = 0; row < 3; ++row)
    {
        for (sofa::Size col = 0; col < 3; ++col)
        {
            const Real value = -hessianPhi(row, col);
            if (!std::isfinite(static_cast<double>(value)))
            {
                hessian.clear();
                return false;
            }
            hessian(row, col) = value;
        }
    }

    return true;
}

template<class TDataTypes1, class TDataTypes2>
bool SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::evaluateLinear8Hessian(
    const Vec3& p,
    Mat3& hessianPhi) const
{
    hessianPhi.clear();

    const Real gx = (p[0] - m_origin[0]) / m_spacing[0];
    const Real gy = (p[1] - m_origin[1]) / m_spacing[1];
    const Real gz = (p[2] - m_origin[2]) / m_spacing[2];

    if (gx < Real(0) || gx >= Real(m_nx - 1) ||
        gy < Real(0) || gy >= Real(m_ny - 1) ||
        gz < Real(0) || gz >= Real(m_nz - 1))
        return false;

    const unsigned int i = static_cast<unsigned int>(std::floor(gx));
    const unsigned int j = static_cast<unsigned int>(std::floor(gy));
    const unsigned int k = static_cast<unsigned int>(std::floor(gz));

    const Real u = gx - Real(i);
    const Real v = gy - Real(j);
    const Real w = gz - Real(k);

    const Real wx[2] = { Real(1) - u, u };
    const Real wy[2] = { Real(1) - v, v };
    const Real wz[2] = { Real(1) - w, w };
    const Real dw[2] = { Real(-1), Real(1) };

    Real duv = Real(0);
    Real duw = Real(0);
    Real dvw = Real(0);

    for (unsigned int a = 0; a < 2; ++a)
    {
        for (unsigned int b = 0; b < 2; ++b)
        {
            for (unsigned int c = 0; c < 2; ++c)
            {
                const Real f = phiAt(i + a, j + b, k + c);
                duv += f * dw[a] * dw[b] * wz[c];
                duw += f * dw[a] * wy[b] * dw[c];
                dvw += f * wx[a] * dw[b] * dw[c];
            }
        }
    }

    const Real hxy = duv / (m_spacing[0] * m_spacing[1]);
    const Real hxz = duw / (m_spacing[0] * m_spacing[2]);
    const Real hyz = dvw / (m_spacing[1] * m_spacing[2]);

    // A trilinear polynomial is affine in each coordinate separately, so the
    // three pure second derivatives are zero inside a cell.
    hessianPhi(0, 0) = Real(0);
    hessianPhi(1, 1) = Real(0);
    hessianPhi(2, 2) = Real(0);
    hessianPhi(0, 1) = hessianPhi(1, 0) = hxy;
    hessianPhi(0, 2) = hessianPhi(2, 0) = hxz;
    hessianPhi(1, 2) = hessianPhi(2, 1) = hyz;

    return true;
}

template<class TDataTypes1, class TDataTypes2>
bool SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::evaluateCubic64Hessian(
    const Vec3& p,
    Mat3& hessianPhi) const
{
    hessianPhi.clear();

    const Real gx = (p[0] - m_origin[0]) / m_spacing[0];
    const Real gy = (p[1] - m_origin[1]) / m_spacing[1];
    const Real gz = (p[2] - m_origin[2]) / m_spacing[2];

    if (gx < Real(0) || gx >= Real(m_nx - 1) ||
        gy < Real(0) || gy >= Real(m_ny - 1) ||
        gz < Real(0) || gz >= Real(m_nz - 1))
        return false;

    if (m_nx < 4 || m_ny < 4 || m_nz < 4)
        return evaluateLinear8Hessian(p, hessianPhi);

    const unsigned int i = static_cast<unsigned int>(std::floor(gx));
    const unsigned int j = static_cast<unsigned int>(std::floor(gy));
    const unsigned int k = static_cast<unsigned int>(std::floor(gz));

    if (i < 1 || i + 2 >= m_nx ||
        j < 1 || j + 2 >= m_ny ||
        k < 1 || k + 2 >= m_nz)
        return evaluateLinear8Hessian(p, hessianPhi);

    const Real u = gx - Real(i);
    const Real v = gy - Real(j);
    const Real w = gz - Real(k);

    auto catmullRomBasis = [](Real t, Real W[4], Real dW[4], Real ddW[4])
    {
        const Real t2 = t * t;
        const Real t3 = t2 * t;

        W[0] = Real(-0.5) * t + t2 - Real(0.5) * t3;
        W[1] = Real(1) - Real(2.5) * t2 + Real(1.5) * t3;
        W[2] = Real(0.5) * t + Real(2) * t2 - Real(1.5) * t3;
        W[3] = Real(-0.5) * t2 + Real(0.5) * t3;

        dW[0] = Real(-0.5) + Real(2) * t - Real(1.5) * t2;
        dW[1] = Real(-5) * t + Real(4.5) * t2;
        dW[2] = Real(0.5) + Real(4) * t - Real(4.5) * t2;
        dW[3] = -t + Real(1.5) * t2;

        ddW[0] = Real(2) - Real(3) * t;
        ddW[1] = Real(-5) + Real(9) * t;
        ddW[2] = Real(4) - Real(9) * t;
        ddW[3] = Real(-1) + Real(3) * t;
    };

    Real wx[4], dwx[4], ddwx[4];
    Real wy[4], dwy[4], ddwy[4];
    Real wz[4], dwz[4], ddwz[4];
    catmullRomBasis(u, wx, dwx, ddwx);
    catmullRomBasis(v, wy, dwy, ddwy);
    catmullRomBasis(w, wz, dwz, ddwz);

    Real duu = Real(0);
    Real dvv = Real(0);
    Real dww = Real(0);
    Real duv = Real(0);
    Real duw = Real(0);
    Real dvw = Real(0);

    for (unsigned int a = 0; a < 4; ++a)
    {
        const unsigned int ii = i + a - 1;
        for (unsigned int b = 0; b < 4; ++b)
        {
            const unsigned int jj = j + b - 1;
            for (unsigned int c = 0; c < 4; ++c)
            {
                const unsigned int kk = k + c - 1;
                const Real f = phiAt(ii, jj, kk);

                duu += f * ddwx[a] * wy[b]   * wz[c];
                dvv += f * wx[a]   * ddwy[b] * wz[c];
                dww += f * wx[a]   * wy[b]   * ddwz[c];
                duv += f * dwx[a]  * dwy[b]  * wz[c];
                duw += f * dwx[a]  * wy[b]   * dwz[c];
                dvw += f * wx[a]   * dwy[b]  * dwz[c];
            }
        }
    }

    const Real hx = m_spacing[0];
    const Real hy = m_spacing[1];
    const Real hz = m_spacing[2];

    hessianPhi(0, 0) = duu / (hx * hx);
    hessianPhi(1, 1) = dvv / (hy * hy);
    hessianPhi(2, 2) = dww / (hz * hz);
    hessianPhi(0, 1) = hessianPhi(1, 0) = duv / (hx * hy);
    hessianPhi(0, 2) = hessianPhi(2, 0) = duw / (hx * hz);
    hessianPhi(1, 2) = hessianPhi(2, 1) = dvw / (hy * hz);

    return true;
}

template<class TDataTypes1, class TDataTypes2>
bool SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::evaluateHermiteCellHessian(
    const Vec3& p,
    bool useMixedDerivatives,
    Mat3& hessianPhi) const
{
    hessianPhi.clear();

    if (!m_hasHermiteFirstDerivatives)
        return false;
    if (useMixedDerivatives && !m_hasHermiteMixedDerivatives)
        return false;

    const Real gx = (p[0] - m_origin[0]) / m_spacing[0];
    const Real gy = (p[1] - m_origin[1]) / m_spacing[1];
    const Real gz = (p[2] - m_origin[2]) / m_spacing[2];

    if (gx < Real(0) || gx >= Real(m_nx - 1) ||
        gy < Real(0) || gy >= Real(m_ny - 1) ||
        gz < Real(0) || gz >= Real(m_nz - 1))
        return false;

    const unsigned int i0 = static_cast<unsigned int>(std::floor(gx));
    const unsigned int j0 = static_cast<unsigned int>(std::floor(gy));
    const unsigned int k0 = static_cast<unsigned int>(std::floor(gz));

    const Real u = gx - Real(i0);
    const Real v = gy - Real(j0);
    const Real w = gz - Real(k0);

    auto hermiteBasis = [](Real t, Real V[2], Real D[2], Real dV[2], Real dD[2], Real ddV[2], Real ddD[2])
    {
        const Real t2 = t * t;
        const Real t3 = t2 * t;

        V[0] = Real(2) * t3 - Real(3) * t2 + Real(1);
        V[1] = Real(-2) * t3 + Real(3) * t2;
        D[0] = t3 - Real(2) * t2 + t;
        D[1] = t3 - t2;

        dV[0] = Real(6) * t2 - Real(6) * t;
        dV[1] = Real(-6) * t2 + Real(6) * t;
        dD[0] = Real(3) * t2 - Real(4) * t + Real(1);
        dD[1] = Real(3) * t2 - Real(2) * t;

        ddV[0] = Real(12) * t - Real(6);
        ddV[1] = Real(-12) * t + Real(6);
        ddD[0] = Real(6) * t - Real(4);
        ddD[1] = Real(6) * t - Real(2);
    };

    Real Vx[2], Dx[2], dVx[2], dDx[2], ddVx[2], ddDx[2];
    Real Vy[2], Dy[2], dVy[2], dDy[2], ddVy[2], ddDy[2];
    Real Vz[2], Dz[2], dVz[2], dDz[2], ddVz[2], ddDz[2];
    hermiteBasis(u, Vx, Dx, dVx, dDx, ddVx, ddDx);
    hermiteBasis(v, Vy, Dy, dVy, dDy, ddVy, ddDy);
    hermiteBasis(w, Vz, Dz, dVz, dDz, ddVz, ddDz);

    const Real hx = m_spacing[0];
    const Real hy = m_spacing[1];
    const Real hz = m_spacing[2];

    Real duu = Real(0);
    Real dvv = Real(0);
    Real dww = Real(0);
    Real duv = Real(0);
    Real duw = Real(0);
    Real dvw = Real(0);

    for (unsigned int a = 0; a < 2; ++a)
    {
        for (unsigned int b = 0; b < 2; ++b)
        {
            for (unsigned int c = 0; c < 2; ++c)
            {
                const unsigned int ii = i0 + a;
                const unsigned int jj = j0 + b;
                const unsigned int kk = k0 + c;

                const Real F = phiAt(ii, jj, kk);
                const Vec3 G = gradAt(ii, jj, kk);

                Real derivative[8] = {
                    F,
                    G[0],
                    G[1],
                    Real(0),
                    G[2],
                    Real(0),
                    Real(0),
                    Real(0)
                };

                if (useMixedDerivatives)
                {
                    const HermiteMixedDerivatives M = mixedAt(ii, jj, kk);
                    derivative[3] = M.fxy;
                    derivative[5] = M.fxz;
                    derivative[6] = M.fyz;
                    derivative[7] = M.fxyz;
                }

                for (unsigned int mask = 0; mask < 8; ++mask)
                {
                    const bool xDerivative = (mask & 1u) != 0;
                    const bool yDerivative = (mask & 2u) != 0;
                    const bool zDerivative = (mask & 4u) != 0;

                    const Real bx = xDerivative ? Dx[a] : Vx[a];
                    const Real by = yDerivative ? Dy[b] : Vy[b];
                    const Real bz = zDerivative ? Dz[c] : Vz[c];
                    const Real dbx = xDerivative ? dDx[a] : dVx[a];
                    const Real dby = yDerivative ? dDy[b] : dVy[b];
                    const Real dbz = zDerivative ? dDz[c] : dVz[c];
                    const Real ddbx = xDerivative ? ddDx[a] : ddVx[a];
                    const Real ddby = yDerivative ? ddDy[b] : ddVy[b];
                    const Real ddbz = zDerivative ? ddDz[c] : ddVz[c];

                    Real coefficient = derivative[mask];
                    if (xDerivative) coefficient *= hx;
                    if (yDerivative) coefficient *= hy;
                    if (zDerivative) coefficient *= hz;

                    duu += coefficient * ddbx * by * bz;
                    dvv += coefficient * bx * ddby * bz;
                    dww += coefficient * bx * by * ddbz;
                    duv += coefficient * dbx * dby * bz;
                    duw += coefficient * dbx * by * dbz;
                    dvw += coefficient * bx * dby * dbz;
                }
            }
        }
    }

    hessianPhi(0, 0) = duu / (hx * hx);
    hessianPhi(1, 1) = dvv / (hy * hy);
    hessianPhi(2, 2) = dww / (hz * hz);
    hessianPhi(0, 1) = hessianPhi(1, 0) = duv / (hx * hy);
    hessianPhi(0, 2) = hessianPhi(2, 0) = duw / (hx * hz);
    hessianPhi(1, 2) = hessianPhi(2, 1) = dvw / (hy * hz);

    return true;
}

template<class TDataTypes1, class TDataTypes2>
typename SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::ContactStatus
SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::computeContactKinematics(const Vec3& p,Contact& contact) const
{
    Real phi = Real(0);
    Vec3 gradPhi(Real(0), Real(0), Real(0));

    if (!evaluate(p, phi, gradPhi))
    {
        contact.gap = Real(0);
        contact.gapGradient = Vec3(Real(0), Real(0), Real(0));
        return ContactStatus::Pinned;
        // return ContactStatus::InvalidGeometry;
    }

    contact.gap = -phi;
    contact.gapGradient = gradPhi * Real(-1);

    return ContactStatus::Active;
}

template<class TDataTypes1, class TDataTypes2>
void SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::draw(const core::visual::VisualParams* vparams)
{
    Base::draw(vparams);

    if (d_showGridBox.getValue())
        drawGridBox(vparams);
}

template<class TDataTypes1, class TDataTypes2>
void SignedDistanceFieldNCPContactForceField<TDataTypes1, TDataTypes2>::drawGridBox(const core::visual::VisualParams* vparams) const
{
    if (!vparams || !vparams->drawTool() || !m_loaded)
        return;

    auto toDrawVec = [](const Vec3& v) -> sofa::type::Vec3
    {
        return sofa::type::Vec3(
            static_cast<float>(v[0]),
            static_cast<float>(v[1]),
            static_cast<float>(v[2]));
    };

    const Vec3 minP = m_origin;
    const Vec3 maxP = Vec3(
        m_origin[0] + m_spacing[0] * Real(m_nx - 1),
        m_origin[1] + m_spacing[1] * Real(m_ny - 1),
        m_origin[2] + m_spacing[2] * Real(m_nz - 1));

    const Vec3 p000(minP[0], minP[1], minP[2]);
    const Vec3 p100(maxP[0], minP[1], minP[2]);
    const Vec3 p010(minP[0], maxP[1], minP[2]);
    const Vec3 p110(maxP[0], maxP[1], minP[2]);
    const Vec3 p001(minP[0], minP[1], maxP[2]);
    const Vec3 p101(maxP[0], minP[1], maxP[2]);
    const Vec3 p011(minP[0], maxP[1], maxP[2]);
    const Vec3 p111(maxP[0], maxP[1], maxP[2]);

    std::vector<sofa::type::Vec3> lines;
    lines.reserve(24);

    auto addEdge = [&](const Vec3& a, const Vec3& b)
    {
        lines.push_back(toDrawVec(a));
        lines.push_back(toDrawVec(b));
    };

    addEdge(p000, p100);
    addEdge(p100, p110);
    addEdge(p110, p010);
    addEdge(p010, p000);

    addEdge(p001, p101);
    addEdge(p101, p111);
    addEdge(p111, p011);
    addEdge(p011, p001);

    addEdge(p000, p001);
    addEdge(p100, p101);
    addEdge(p110, p111);
    addEdge(p010, p011);

    vparams->drawTool()->drawLines(lines, 1.0f, d_color.getValue());
}

} // namespace sofa::ncp
