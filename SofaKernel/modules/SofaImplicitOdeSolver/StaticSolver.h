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
// Author: François Faure, INRIA-UJF, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
#ifndef SOFA_COMPONENT_ODESOLVER_STATICSOLVER_H
#define SOFA_COMPONENT_ODESOLVER_STATICSOLVER_H
#include "config.h"

#include <sofa/core/behavior/OdeSolver.h>


namespace sofa
{

namespace component
{

namespace odesolver
{

/** Finds the static equilibrium of a system. Can diverge when there are an infinity of solutions. */
class SOFA_IMPLICIT_ODE_SOLVER_API StaticSolver : public sofa::core::behavior::OdeSolver
{

public:
    SOFA_CLASS(StaticSolver, sofa::core::behavior::OdeSolver);
protected:
    StaticSolver();
public:
    void solve (const core::ExecParams* params, SReal dt,
                sofa::core::MultiVecCoordId xResult, sofa::core::MultiVecDerivId vResult) override;

    /// Given an input derivative order (0 for position, 1 for velocity, 2 for acceleration),
    /// how much will it affect the output derivative of the given order.
    double getIntegrationFactor(int inputDerivative, int outputDerivative) const override;

    /// Given a solution of the linear system,
    /// how much will it affect the output derivative of the given order.
    double getSolutionIntegrationFactor(int outputDerivative) const override;

    Data<SReal> massCoef; ///< factor associated with the mass matrix in the equation system
    Data<SReal> dampingCoef; ///< factor associated with the mass matrix in the equation system
    Data<SReal> stiffnessCoef; ///< factor associated with the mass matrix in the equation system
    Data<bool> applyIncrementFactor; ///< multiply the solution by dt. Default: false
    Data<bool> d_threadSafeVisitor;
};

} // namespace odesolver

} // namespace component

} // namespace sofa

#endif
