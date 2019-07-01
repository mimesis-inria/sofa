/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2019 INRIA, USTL, UJF, CNRS, MGH                    *
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
#ifndef SOFA_CORE_BEHAVIOR_SingleMatrixAccessor_H
#define SOFA_CORE_BEHAVIOR_SingleMatrixAccessor_H
#include "config.h"

#include <sofa/core/behavior/MultiMatrixAccessor.h>
#include <sofa/defaulttype/BaseMatrix.h>

namespace sofa
{

namespace component
{

namespace linearsolver
{

/** Special case to access a single square matrix.
*/
class SOFA_BASE_LINEAR_SOLVER_API SingleMatrixAccessor : public core::behavior::MultiMatrixAccessor
{
public:
    typedef defaulttype::BaseMatrix BaseMatrix;

    SingleMatrixAccessor( BaseMatrix* m=nullptr ) { setMatrix(m); }
    ~SingleMatrixAccessor() override;

    void setMatrix( BaseMatrix* m );
    BaseMatrix* getMatrix() { return matrix; }
    const BaseMatrix* getMatrix() const { return matrix; }


    int getGlobalDimension() const override { return matrix->rowSize(); }
    int getGlobalOffset(const core::behavior::BaseMechanicalState*) const override { return 0; }
    MatrixRef getMatrix(const core::behavior::BaseMechanicalState*) const override;


    InteractionMatrixRef getMatrix(const core::behavior::BaseMechanicalState* mstate1, const core::behavior::BaseMechanicalState* mstate2) const override;

protected:
    BaseMatrix* matrix;   ///< The single matrix
    MatrixRef matRef; ///< The accessor to the single matrix

};

} // namespace behavior

} // namespace core

} // namespace sofa

#endif
