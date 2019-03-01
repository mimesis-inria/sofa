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
#include <SofaImplicitOdeSolver/StaticSolver.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/MechanicalVisitor.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/VectorOperations.h>
#include <sofa/core/ObjectFactory.h>
#include <math.h>
#include <iostream>





namespace sofa
{

namespace component
{

namespace odesolver
{
using core::VecId;
using namespace sofa::defaulttype;
using namespace core::behavior;

StaticSolver::StaticSolver()
    : massCoef( initData(&massCoef,(SReal)0.0,"massCoef","factor associated with the mass matrix in the equation system") )
    , dampingCoef( initData(&dampingCoef,(SReal)0.0,"dampingCoef","factor associated with the mass matrix in the equation system") )
    , stiffnessCoef( initData(&stiffnessCoef,(SReal)1.0,"stiffnessCoef","factor associated with the mass matrix in the equation system") )
    , applyIncrementFactor( initData(&applyIncrementFactor,false,"applyIncrementFactor","multiply the solution by dt before adding it to the current state") )
    , d_threadSafeVisitor(initData(&d_threadSafeVisitor, false, "threadSafeVisitor", "If true, do not use realloc and free visitors in fwdInteractionForceField."))
{
}

void StaticSolver::solve(const core::ExecParams* params, SReal dt, sofa::core::MultiVecCoordId xResult, sofa::core::MultiVecDerivId /*vResult*/)
{
    sofa::simulation::common::VectorOperations vop( params, this->getContext() );
    sofa::simulation::common::MechanicalOperations mop( params, this->getContext() );
    MultiVecCoord pos(&vop, core::VecCoordId::position() );
    MultiVecDeriv force(&vop, core::VecDerivId::force() );
    MultiVecCoord pos2(&vop, xResult /*core::VecCoordId::position()*/ );
    //MultiVecDeriv vel2(&vop, vResult /*core::VecDerivId::velocity()*/ );

    //MultiVecDeriv b(&vop);
    MultiVecDeriv x(&vop);

    // dx is no longer allocated by default (but it will be deleted automatically by the mechanical objects)
    MultiVecDeriv dx(&vop, core::VecDerivId::dx()); dx.realloc(&vop, !d_threadSafeVisitor.getValue(), true);
    mop->setImplicit(true); // this solver is implicit
    mop.addSeparateGravity(dt);	// v += dt*g . Used if mass wants to add G to v separately from the other forces.

    // compute the right-hand term of the equation system
    mop.computeForce(force);             // b = f0
    mop.projectResponse(force);         // b is projected to the constrained space
    //    b.teq(-1);

    dmsg_info() << "StaticSolver, f0 = "<< force ;
    core::behavior::MultiMatrix<simulation::common::MechanicalOperations> matrix(&mop);
    //matrix = MechanicalMatrix::K;
    matrix = MechanicalMatrix(massCoef.getValue(),dampingCoef.getValue(),stiffnessCoef.getValue());

    dmsg_info() <<" matrix = "<< (MechanicalMatrix::K) << " = " << matrix ;

    matrix.solve(x,force);
    // x is the opposite solution of the system

    // apply the solution
    /*    serr<<"StaticSolver::solve, nb iter = "<<nb_iter<<sendl;
     serr<<"StaticSolver::solve, solution = "<<x<<sendl;*/

    dmsg_info() <<" opposite solution = "<< x ;

    if(applyIncrementFactor.getValue()==true )
        pos2.eq( pos, x, -dt );
    else
        pos2.eq( pos, x, -1 );

    mop.solveConstraint(pos2, core::ConstraintParams::POS);

}

/// Given an input derivative order (0 for position, 1 for velocity, 2 for acceleration),
/// how much will it affect the output derivative of the given order.
double StaticSolver::getIntegrationFactor(int inputDerivative, int outputDerivative) const
{
    double matrix[3][3] =
    {
        { 1, 0, 0},
        { 0, 1, 0},
        { 0, 0, 0}
    };
    if (inputDerivative >= 3 || outputDerivative >= 3)
        return 0;
    else
        return matrix[outputDerivative][inputDerivative];
}

/// Given a solution of the linear system,
/// how much will it affect the output derivative of the given order.
double StaticSolver::getSolutionIntegrationFactor(int outputDerivative) const
{
    double vect[3] = { 1, 0, 0};
    if (outputDerivative >= 3)
        return 0;
    else
        return vect[outputDerivative];
}



int StaticSolverClass = core::RegisterObject("A solver which seeks the static equilibrium of the scene it monitors")
        .add< StaticSolver >();

} // namespace odesolver

} // namespace component

} // namespace sofa

