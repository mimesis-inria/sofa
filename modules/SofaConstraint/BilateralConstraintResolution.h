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
#ifndef SOFA_COMPONENT_CONSTRAINTSET_BILATERALCONSTRAINTRESOLUTION_H
#define SOFA_COMPONENT_CONSTRAINTSET_BILATERALCONSTRAINTRESOLUTION_H
#include "config.h"


#include <sofa/defaulttype/Mat.h>
#include <sofa/core/behavior/BaseConstraint.h>

#include <Eigen/Core>
#include <Eigen/Cholesky>

namespace sofa
{

namespace component
{

namespace constraintset
{

namespace bilateralconstraintresolution
{

using sofa::core::behavior::ConstraintResolution ;


class BilateralConstraintResolution : public ConstraintResolution
{
public:
    BilateralConstraintResolution(double* initF=NULL) 
        : ConstraintResolution(1)
        , _f(initF) {}
    virtual void resolution(int line, double** w, double* d, double* force, double *dfree)
    {
        SOFA_UNUSED(dfree);
        force[line] -= d[line] / w[line][line];
    }

    virtual void init(int line, double** /*w*/, double* force)
    {
        if(_f) { force[line] = *_f; }
    }

    virtual void initForce(int line, double* force)
    {
        if(_f) { force[line] = *_f; }
    }

    void store(int line, double* force, bool /*convergence*/)
    {
        if(_f) *_f = force[line];
    }

protected:
    double* _f;
};

class BilateralConstraintResolution3Dof : public ConstraintResolution
{
public:

    BilateralConstraintResolution3Dof(sofa::defaulttype::Vec3d* vec = NULL)
        : ConstraintResolution(3)
        , _f(vec)
    {
    }
    virtual void init(int line, double** w, double *force)
    {
        sofa::defaulttype::Mat<3,3,double> temp;
        temp[0][0] = w[line][line];
        temp[0][1] = w[line][line+1];
        temp[0][2] = w[line][line+2];
        temp[1][0] = w[line+1][line];
        temp[1][1] = w[line+1][line+1];
        temp[1][2] = w[line+1][line+2];
        temp[2][0] = w[line+2][line];
        temp[2][1] = w[line+2][line+1];
        temp[2][2] = w[line+2][line+2];

        invertMatrix(invW, temp);

        if(_f)
        {
            for(int i=0; i<3; i++)
                force[line+i] = (*_f)[i];
        }
    }

    virtual void initForce(int line, double* force)
    {
        if(_f)
        {
            for(int i=0; i<3; i++)
                force[line+i] = (*_f)[i];
        }
    }

    virtual void resolution(int line, double** /*w*/, double* d, double* force, double * dFree)
    {
        SOFA_UNUSED(dFree);
        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
                force[line+i] -= d[line+j] * invW[i][j];
        }
    }

    void store(int line, double* force, bool /*convergence*/)
    {
        if(_f)
        {
            for(int i=0; i<3; i++)
                (*_f)[i] = force[line+i];
        }
    }

protected:
    sofa::defaulttype::Mat<3,3,double> invW;
    sofa::defaulttype::Vec3d* _f;
};

class BilateralConstraintResolutionNDof : public ConstraintResolution
{
public:
    BilateralConstraintResolutionNDof(unsigned blockSize ) 
    : ConstraintResolution(blockSize)
    , wBlock(Eigen::MatrixXd(blockSize, blockSize))
    {
    }
    virtual void init(int line, double** w, double * /*force*/)
    {
        for (auto i = 0; i < wBlock.rows(); ++i)   
        {
            wBlock.row(i) = Eigen::VectorXd::Map(&w[line + i][line], wBlock.cols());
        }
        wBlockInv.compute(wBlock);
    }

    virtual void resolution(int line, double** /*w*/, double* displacement, double* force, double* /*dFree*/)
    {
        Eigen::Map< Eigen::VectorXd > f(&force[line], wBlock.cols());
        Eigen::Map< Eigen::VectorXd > d(&displacement[line], wBlock.cols());
        Eigen::VectorXd f_local = wBlockInv.solve(d);
        f -= f_local;
    }

protected:
    Eigen::MatrixXd  wBlock;
    Eigen::LDLT< Eigen::MatrixXd > wBlockInv;
};

} // namespace bilateralconstraintresolution

/// Import the following into the constraintset namespace to preserve
/// compatibility with the existing sofa source code.
using bilateralconstraintresolution::BilateralConstraintResolution ;
using bilateralconstraintresolution::BilateralConstraintResolution3Dof ;
using bilateralconstraintresolution::BilateralConstraintResolutionNDof ;


} // namespace constraintset

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONSTRAINTSET_BILATERALCONSTRAINTRESOLUTION_H
