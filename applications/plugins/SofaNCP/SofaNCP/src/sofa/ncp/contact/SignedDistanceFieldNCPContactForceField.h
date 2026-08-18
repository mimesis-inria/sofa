/****************************************************************************
* SDF Fischer-Burmeister NCP contact forcefield.
*
* Geometry specialization only. The generic FB/KKT algebra stays in
* FischerBurmeisterContactForceField.
*
* Sign convention of the raw scalar field:
*   phi < 0  inside lumen
*   phi = 0  wall
*   phi > 0  outside lumen / violation
*
* NCP convention given to the base class:
*   gap         = -phi
*   gapGradient = -grad(phi)
*
* Interpolation modes:
*   0 Auto
*       Use full tricubic Hermite if full derivatives are loaded; otherwise
*       use phi-only Catmull-Rom cubic interpolation.
*
*   1 Linear8
*       Trilinear phi and exact derivative of that trilinear interpolant.
*       Mostly useful as a debugging baseline.
*
*   2 Cubic64
*       Phi-only 4x4x4 Catmull-Rom cubic interpolation. The returned gradient
*       is the derivative of the same cubic interpolant. This is the preferred
*       simple mode for the next SDF contact pass.
*
*   3 HermiteFirstDerivatives
*       Tensor-product cubic Hermite using only f, fx, fy, fz at the 8 cell
*       corners. Mixed derivatives are assumed zero. This is kept as a debug
*       / compatibility mode, but it is not full tricubic Hermite.
*
*   4 HermiteFull
*       Full tensor-product tricubic Hermite using f, fx, fy, fz, fxy, fxz,
*       fyz, fxyz at the 8 cell corners.
*
* Raw files:
*   sdfFilename                 : float32 phi[ix,iy,iz], C-order, z fastest.
*   sdfGradientFilename          : optional float32 [fx,fy,fz] per grid node.
*   sdfHermiteMixedFilename      : optional float32 [fxy,fxz,fyz,fxyz] per node.
*   sdfHermitePacketFilename     : optional float32 [fx,fy,fz,fxy,fxz,fyz,fxyz]
*                                  per node. If provided, it overrides the two
*                                  derivative files above.
****************************************************************************/
#pragma once

#include <sofa/ncp/config.h>
#include <sofa/ncp/contact/FischerBurmeisterContactForceField.h>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/type/RGBAColor.h>
#include <sofa/type/Vec.h>
#include <sofa/type/vector.h>

#include <string>
#include <vector>

namespace sofa::core { class ObjectFactory; }

namespace sofa::ncp
{

template<class TDataTypes1, class TDataTypes2>
class SignedDistanceFieldNCPContactForceField
    : public FischerBurmeisterContactForceField<TDataTypes1, TDataTypes2>
{
public:
    SOFA_CLASS(
        SOFA_TEMPLATE2(SignedDistanceFieldNCPContactForceField, TDataTypes1, TDataTypes2),
        SOFA_TEMPLATE2(FischerBurmeisterContactForceField, TDataTypes1, TDataTypes2)
    );

    using Base = FischerBurmeisterContactForceField<TDataTypes1, TDataTypes2>;
    using Inherit = Base;

    using DataTypes1 = TDataTypes1;
    using DataTypes2 = TDataTypes2;
    using Real = typename Base::Real;
    using Vec3 = typename Base::Vec3;
    using Mat3 = typename Base::Mat3;
    using Contact = typename Base::Contact;
    using ContactStatus = typename Base::ContactStatus;

    using UInt3 = sofa::type::Vec<3, unsigned int>;

    enum InterpolationMode : unsigned int
    {
        Auto = 0,
        Linear8 = 1,
        Cubic64 = 2,
        HermiteFirstDerivatives = 3,
        HermiteFull = 4
    };

    struct HermiteMixedDerivatives
    {
        Real fxy = Real(0);
        Real fxz = Real(0);
        Real fyz = Real(0);
        Real fxyz = Real(0);
    };

    Data<std::string> d_sdfFilename;
    Data<std::string> d_sdfGradientFilename;
    Data<std::string> d_sdfHermiteMixedFilename;
    Data<std::string> d_sdfHermitePacketFilename;

    Data<Vec3> d_origin;
    Data<Vec3> d_spacing;
    Data<UInt3> d_dimensions;

    /// 0=Auto, 1=Linear8, 2=Cubic64, 3=HermiteFirstDerivatives, 4=HermiteFull.
    Data<unsigned int> d_interpolationMode;

    /// Normalize grad(phi) before returning gapGradient. Keep false for an exact derivative of the interpolated gap.
    Data<bool> d_normalizeGradient;

    Data<bool> d_debugQueries;
    Data<bool> d_showGridBox;
    Data<sofa::type::RGBAColor> d_color;

protected:
    SignedDistanceFieldNCPContactForceField();

    SignedDistanceFieldNCPContactForceField(
        core::behavior::MechanicalState<DataTypes1>* object1,
        core::behavior::MechanicalState<DataTypes2>* object2);

public:
    void init() override;
    void reinit() override;
    void draw(const core::visual::VisualParams* vparams) override;

protected:
    ContactStatus computeContactKinematics(const Vec3& p, Contact& contact) const override;
    bool computeGapHessian(const Vec3& p, Mat3& hessian) const override;

    bool loadGridData();
    void clearGridData();

    bool evaluate(const Vec3& p, Real& phi, Vec3& gradPhi) const;
    bool evaluateLinear8(const Vec3& p, Real& phi, Vec3& gradPhi) const;
    bool evaluateCubic64(const Vec3& p, Real& phi, Vec3& gradPhi) const;
    bool evaluateHermiteFirst(const Vec3& p, Real& phi, Vec3& gradPhi) const;
    bool evaluateHermiteFull(const Vec3& p, Real& phi, Vec3& gradPhi) const;
    bool evaluateHermiteCell(const Vec3& p, bool useMixedDerivatives, Real& phi, Vec3& gradPhi) const;

    bool evaluateLinear8Hessian(const Vec3& p, Mat3& hessianPhi) const;
    bool evaluateCubic64Hessian(const Vec3& p, Mat3& hessianPhi) const;
    bool evaluateHermiteCellHessian(const Vec3& p, bool useMixedDerivatives, Mat3& hessianPhi) const;

    bool validateAndNormalizeGradient(const Vec3& p, Real phi, Vec3& gradPhi) const;

    sofa::Size flattenedIndex(unsigned int i, unsigned int j, unsigned int k) const;

    Real phiAt(unsigned int i, unsigned int j, unsigned int k) const;
    Vec3 gradAt(unsigned int i, unsigned int j, unsigned int k) const;
    HermiteMixedDerivatives mixedAt(unsigned int i, unsigned int j, unsigned int k) const;

    static const char* interpolationModeName(unsigned int mode);

    void drawGridBox(const core::visual::VisualParams* vparams) const;

    std::vector<float> m_phi;
    std::vector<float> m_gradientData;
    std::vector<float> m_mixedDerivativeData;

    unsigned int m_nx = 0;
    unsigned int m_ny = 0;
    unsigned int m_nz = 0;

    Vec3 m_origin = Vec3(Real(0), Real(0), Real(0));
    Vec3 m_spacing = Vec3(Real(1), Real(1), Real(1));
    bool m_loaded = false;
    bool m_hasHermiteFirstDerivatives = false;
    bool m_hasHermiteMixedDerivatives = false;

    mutable bool m_warnedMissingFirstDerivatives = false;
    mutable bool m_warnedMissingFullHermite = false;
    mutable bool m_warnedCubicBoundaryFallback = false;
};

#if !defined(SOFANCP_SDF_CONTACT_FORCE_FIELD_CPP)
extern template class SOFANCP_API SignedDistanceFieldNCPContactForceField<defaulttype::Vec3Types, defaulttype::Vec1Types>;
extern template class SOFANCP_API SignedDistanceFieldNCPContactForceField<defaulttype::Rigid3Types, defaulttype::Vec1Types>;
#endif

SOFANCP_API void registerSignedDistanceFieldNCPContactForceField(sofa::core::ObjectFactory* factory);

} // namespace sofa::ncp
