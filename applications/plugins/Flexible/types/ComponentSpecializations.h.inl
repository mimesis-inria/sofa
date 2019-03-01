

#include <Flexible/config.h>


#include <SofaBaseMechanics/MechanicalObject.h>


#include <SofaBaseMechanics/AddMToMatrixFunctor.h>
#include <SofaBaseMechanics/UniformMass.h>

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/Mass.h>
#include <sofa/core/behavior/ProjectiveConstraintSet.inl>
#include <sofa/core/behavior/ConstraintCorrection.inl>
#include <sofa/core/Mapping.h>


#ifdef SOFA_HAVE_IMAGE
#include "../mass/ImageDensityMass.h"
#endif


#include "ComponentSpecializationsDefines.h"





namespace sofa
{


// ==========================================================================
// Mechanical Object

namespace component
{
namespace container
{

template <> SOFA_Flexible_API
void MechanicalObject<defaulttype::TYPEABSTRACTNAME3dTypes>::draw(const core::visual::VisualParams* vparams);


#if !defined(FLEXIBLE_COMPILING_CPP)
extern template class SOFA_Flexible_API MechanicalObjectInternalData<defaulttype::TYPEABSTRACTNAME3dTypes>;
extern template class SOFA_Flexible_API MechanicalObject<defaulttype::TYPEABSTRACTNAME3dTypes>;
#endif


} // namespace container







// ==========================================================================
// Uniform Mass


namespace mass
{

template<int N, typename Real>
class AddMToMatrixFunctor< typename defaulttype::StdTYPEABSTRACTNAMETypes<N,Real>::Deriv, defaulttype::DeformableFrameMass<N,defaulttype::StdTYPEABSTRACTNAMETypes<N,Real>::deriv_total_size,Real> >
{
public:
    void operator()(defaulttype::BaseMatrix * mat, const defaulttype::DeformableFrameMass<N,defaulttype::StdTYPEABSTRACTNAMETypes<N,Real>::deriv_total_size,Real>& mass, int pos, SReal fact)
    {
        typedef defaulttype::DeformableFrameMass<N,defaulttype::StdTYPEABSTRACTNAMETypes<N,Real>::deriv_total_size,Real> TYPEABSTRACTNAMEMass;
        for( unsigned i=0; i<TYPEABSTRACTNAMEMass::VSize; ++i )
            for( unsigned j=0; j<TYPEABSTRACTNAMEMass::VSize; ++j )
            {
                mat->add(pos+i, pos+j, mass[i][j]*fact);
                //            cerr<<"AddMToMatrixFunctor< defaulttype::Vec<N,Real>, defaulttype::Mat<N,N,Real> >::operator(), add "<< mass[i][j]*fact << " in " << pos+i <<","<< pos+j <<endl;
            }
    }
};


template <> SOFA_Flexible_API
void UniformMass<defaulttype::TYPEABSTRACTNAME3dTypes, defaulttype::TYPEABSTRACTNAME3dMass>::constructor_message();
template <> SOFA_Flexible_API
void UniformMass<defaulttype::TYPEABSTRACTNAME3dTypes, defaulttype::TYPEABSTRACTNAME3dMass>::draw( const core::visual::VisualParams* vparams );
template <> SOFA_Flexible_API
SReal UniformMass<defaulttype::TYPEABSTRACTNAME3dTypes, defaulttype::TYPEABSTRACTNAME3dMass>::getPotentialEnergy( const core::MechanicalParams*, const DataVecCoord& vx ) const;





#if !defined(FLEXIBLE_COMPILING_CPP)
#ifdef SOFA_HAVE_IMAGE
extern template class SOFA_Flexible_API ImageDensityMass<defaulttype::TYPEABSTRACTNAME3dTypes,core::behavior::ShapeFunction3d,defaulttype::TYPEABSTRACTNAME3dMass>;

extern template class SOFA_Flexible_API UniformMass<defaulttype::TYPEABSTRACTNAME3dTypes,defaulttype::TYPEABSTRACTNAME3dMass>;
#endif

#endif



} // namespace mass



} // namespace component



namespace core
{

namespace behavior
{

#if !defined(FLEXIBLE_COMPILING_CPP)
extern template class SOFA_Flexible_API ForceField<defaulttype::TYPEABSTRACTNAME3dTypes>;
extern template class SOFA_Flexible_API Mass<defaulttype::TYPEABSTRACTNAME3dTypes>;
extern template class SOFA_Flexible_API ConstraintCorrection<defaulttype::TYPEABSTRACTNAME3dTypes>;
extern template class SOFA_Flexible_API ProjectiveConstraintSet<defaulttype::TYPEABSTRACTNAME3dTypes>;
#endif


} // namespace behavior



#if !defined(FLEXIBLE_COMPILING_CPP)
extern template class SOFA_Flexible_API Mapping<defaulttype::TYPEABSTRACTNAME3dTypes,defaulttype::Rigid3Types>;

#endif


} // namespace core




}// namespace sofa


#include "ComponentSpecializationsUndef.h"

