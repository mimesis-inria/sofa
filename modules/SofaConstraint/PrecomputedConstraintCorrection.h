/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2018 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_CORE_COLLISION_CONTACTCORRECTION_H
#define SOFA_CORE_COLLISION_CONTACTCORRECTION_H
#include "config.h"

#include <sofa/core/behavior/ConstraintCorrection.h>
#include <sofa/core/objectmodel/DataFileName.h>

#include <SofaBaseLinearSolver/FullMatrix.h>

#include <sofa/defaulttype/Mat.h>
#include <sofa/defaulttype/Vec.h>

namespace sofa
{

namespace component
{

namespace constraintset
{

/**
 *  \brief Component computing constraint forces within a simulated body using the compliance method.
 */
template<class TDataTypes>
class PrecomputedConstraintCorrection : public sofa::core::behavior::ConstraintCorrection< TDataTypes >
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(PrecomputedConstraintCorrection,TDataTypes), SOFA_TEMPLATE(core::behavior::ConstraintCorrection, TDataTypes));

    typedef TDataTypes DataTypes;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::MatrixDeriv::RowConstIterator MatrixDerivRowConstIterator;
    typedef typename DataTypes::MatrixDeriv::ColConstIterator MatrixDerivColConstIterator;
    typedef typename DataTypes::MatrixDeriv::RowIterator MatrixDerivRowIterator;
    typedef typename DataTypes::MatrixDeriv::ColIterator MatrixDerivColIterator;

    typedef sofa::core::behavior::ConstraintCorrection< TDataTypes > Inherit;

    typedef typename Coord::value_type Real;
    typedef sofa::defaulttype::MatNoInit<3, 3, Real> Transformation;

    Data<bool> m_rotations;
    Data<bool> m_restRotations;

    Data<bool> recompute; ///< if true, always recompute the compliance
	Data<double> debugViewFrameScale; ///< Scale on computed node's frame
	sofa::core::objectmodel::DataFileName f_fileCompliance; ///< Precomputed compliance matrix data file
	Data<std::string> fileDir; ///< If not empty, the compliance will be saved in this repertory
    
protected:
    PrecomputedConstraintCorrection(sofa::core::behavior::MechanicalState<DataTypes> *mm = NULL);

    virtual ~PrecomputedConstraintCorrection();
public:
    virtual void bwdInit() override;

    virtual void addComplianceInConstraintSpace(const sofa::core::ConstraintParams *cparams, sofa::defaulttype::BaseMatrix* W) override;

    virtual void getComplianceMatrix(defaulttype::BaseMatrix* m) const override;

    virtual void computeMotionCorrection(const core::ConstraintParams*, core::MultiVecDerivId dx, core::MultiVecDerivId f) override;

    virtual void applyMotionCorrection(const sofa::core::ConstraintParams *cparams, sofa::Data< VecCoord > &x, sofa::Data< VecDeriv > &v, sofa::Data< VecDeriv > &dx , const sofa::Data< VecDeriv > & correction) override;

    virtual void applyPositionCorrection(const sofa::core::ConstraintParams *cparams, sofa::Data< VecCoord > &x, sofa::Data< VecDeriv > &dx, const sofa::Data< VecDeriv > & correction) override;

    virtual void applyVelocityCorrection(const sofa::core::ConstraintParams *cparams, sofa::Data< VecDeriv > &v, sofa::Data< VecDeriv > &dv, const sofa::Data< VecDeriv > & correction) override;

    /// @name Deprecated API
    /// @{

    virtual void applyContactForce(const defaulttype::BaseVector *f) override;

    virtual void resetContactForce() override;

    /// @}

    virtual void rotateConstraints(bool back);

    virtual void rotateResponse();

    virtual void draw(const core::visual::VisualParams* vparams) override;

    /// @name Unbuilt constraint system during resolution
    /// @{

    virtual void resetForUnbuiltResolution(double * f, std::list<unsigned int>& /*renumbering*/) override;

    virtual bool hasConstraintNumber(int index) override;  // virtual ???

    virtual void addConstraintDisplacement(double *d, int begin,int end) override;

    virtual void setConstraintDForce(double *df, int begin, int end, bool update) override;

    virtual void getBlockDiagonalCompliance(defaulttype::BaseMatrix* W, int begin, int end) override;

    /// @}

public:

    struct InverseStorage
    {
        Real* data;
        int nbref;
        InverseStorage() : data(NULL), nbref(0) {}
    };

    std::string invName;
    InverseStorage* invM;
    Real* appCompliance;
    unsigned int dimensionAppCompliance;

    static std::map<std::string, InverseStorage>& getInverseMap()
    {
        static std::map<std::string, InverseStorage> registry;
        return registry;
    }

    static InverseStorage* getInverse(std::string name);

    static void releaseInverse(std::string name, InverseStorage* inv);

    unsigned int nbRows, nbCols, dof_on_node, nbNodes;
    helper::vector<int> _indexNodeSparseCompliance;
    helper::vector<Deriv> _sparseCompliance;
    Real Fbuf[6], DXbuf;

    // new :  for non building the constraint system during solving process //
    //VecDeriv constraint_disp, constraint_force;
    helper::vector<int> id_to_localIndex;	// table that gives the local index of a constraint given its id
    helper::vector<unsigned int> localIndex_to_id; //inverse table that gives the id of a constraint given its local index
    std::list<unsigned int> active_local_force; // table of local index of the non-null forces;
    linearsolver::FullMatrix< Real > localW;
    double* constraint_force;

    // NEW METHOD FOR UNBUILT
    // new :  for non building the constraint system during solving process //
    VecDeriv constraint_D, constraint_F;
    std::list<int> constraint_dofs;		// list of indices of each point which is involve with constraint

public:
    Real* getInverse()
    {
        if(invM->data)
            return invM->data;
        else
            serr<<"Inverse is not computed yet"<<sendl;
        return NULL;
    }

protected:
    /**
     * @brief Load compliance matrix from memory or external file according to fileName.
     *
     * @return Loading success.
     */
    bool loadCompliance(std::string fileName);

    /**
     * @brief Save compliance matrix into a file.
     */
    void saveCompliance(const std::string& fileName);

    /**
     * @brief Builds the compliance file name using the SOFA component internal data.
     */
    std::string buildFileName();

    /**
     * @brief Compute dx correction from motion space force vector.
     */
    void computeDx(Data<VecDeriv>& dx, const Data< VecDeriv > &f, const std::list< int > &activeDofs);


    /// the list of 
    std::list< int > m_activeDofs;
};


/////////////////////////////////////////////////////////////////////////////////


#ifndef SOFA_FLOAT

template<>
void PrecomputedConstraintCorrection<defaulttype::Rigid3dTypes>::rotateConstraints(bool back);

template<>
void PrecomputedConstraintCorrection<defaulttype::Vec1dTypes>::rotateConstraints(bool back);

template<>
void PrecomputedConstraintCorrection<defaulttype::Rigid3dTypes>::rotateResponse();

template<>
void PrecomputedConstraintCorrection<defaulttype::Vec1dTypes>::rotateResponse();

template<>
void PrecomputedConstraintCorrection<defaulttype::Rigid3dTypes>::draw(const core::visual::VisualParams* vparams);

template<>
void PrecomputedConstraintCorrection<defaulttype::Vec1dTypes>::draw(const core::visual::VisualParams* vparams);

#endif

#ifndef SOFA_DOUBLE

template<>
void PrecomputedConstraintCorrection<defaulttype::Rigid3fTypes>::rotateConstraints(bool back);

template<>
void PrecomputedConstraintCorrection<defaulttype::Vec1fTypes>::rotateConstraints(bool back);

template<>
void PrecomputedConstraintCorrection<defaulttype::Rigid3fTypes>::rotateResponse();

template<>
void PrecomputedConstraintCorrection<defaulttype::Vec1fTypes>::rotateResponse();

template<>
void PrecomputedConstraintCorrection<defaulttype::Rigid3fTypes>::draw(const core::visual::VisualParams* vparams);

template<>
void PrecomputedConstraintCorrection<sofa::defaulttype::Vec1fTypes>::draw(const sofa::core::visual::VisualParams* vparams);

#endif

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_CONSTRAINTSET_PRECOMPUTEDCONSTRAINTCORRECTION_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_CONSTRAINT_API PrecomputedConstraintCorrection<defaulttype::Vec3dTypes>;
extern template class SOFA_CONSTRAINT_API PrecomputedConstraintCorrection<defaulttype::Vec1dTypes>;
extern template class SOFA_CONSTRAINT_API PrecomputedConstraintCorrection<defaulttype::Rigid3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_CONSTRAINT_API PrecomputedConstraintCorrection<defaulttype::Vec3fTypes>;
extern template class SOFA_CONSTRAINT_API PrecomputedConstraintCorrection<defaulttype::Vec1fTypes>;
extern template class SOFA_CONSTRAINT_API PrecomputedConstraintCorrection<defaulttype::Rigid3fTypes>;
#endif
#endif


} // namespace collision

} // namespace component

} // namespace sofa

#endif
