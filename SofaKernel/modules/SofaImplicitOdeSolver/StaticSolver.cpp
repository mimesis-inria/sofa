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

#include <SofaImplicitOdeSolver/StaticSolver.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/simulation/MechanicalOperations.h>
#include <sofa/simulation/VectorOperations.h>
#include <sofa/simulation/PropagateEventVisitor.h>

namespace sofa
{

namespace component
{

namespace odesolver
{

using sofa::core::VecId;
using namespace sofa::defaulttype;
using namespace sofa::core::behavior;

StaticSolver::StaticSolver()
    : d_newton_iterations(initData(&d_newton_iterations,
            (unsigned) 1,
            "newton_iterations",
            "Number of newton iterations between each load increments (normally, one load increment per simulation time-step."))
    , d_correction_tolerance_threshold(initData(&d_correction_tolerance_threshold,
            (double) 1e-5,
            "correction_tolerance_threshold",
            "Convergence criterion: The newton iterations will stop when the norm of correction |du| reach this threshold."))
    , d_residual_tolerance_threshold( initData(&d_residual_tolerance_threshold,
            (double) 1e-5,
            "residual_tolerance_threshold",
            "Convergence criterion: The newton iterations will stop when the norm of the residual |f - K(u)| reach this threshold. "
            "Use a negative value to disable this criterion."))
    , d_should_diverge_when_residual_is_growing( initData(&d_should_diverge_when_residual_is_growing,
            false,
            "should_diverge_when_residual_is_growing",
            "Divergence criterion: The newton iterations will stop when the residual is greater than the one from the previous iteration."))
{}

void StaticSolver::parse(sofa::core::objectmodel::BaseObjectDescription* arg)
{
    /// Now handling backward compatibility with old scenes.
    /// point is deprecated since '19.06'
    /// massCoef, dampingCoef, stiffnessCoef, threadSafeVisitor
    const char* val=arg->getAttribute("massCoef",nullptr) ;
    if(val)
    {
        msg_deprecated() << "The attribute 'massCoef' is deprecated since Sofa 19.06'" << msgendl
                         << "This data was previously used for stabilization purposes but it was preventing" << msgendl
                         << "from computing a strictly-static system (Use the forum for any question)";

    }
    val=arg->getAttribute("dampingCoef",nullptr) ;
    if(val)
    {
        msg_deprecated() << "The attribute 'dampingCoef' is deprecated since Sofa 19.06'" << msgendl
                         << "This data was previously used for stabilization purposes but it was preventing" << msgendl
                         << "from computing a strictly-static system (Use the forum for any question)";

    }
    val=arg->getAttribute("stiffnessCoef",nullptr) ;
    if(val)
    {
        msg_deprecated() << "The attribute 'stiffnessCoef' is deprecated since Sofa 19.06'" << msgendl
                         << "This data was previously used for stabilization purposes but it was preventing" << msgendl
                         << "from computing a strictly-static system (Use the forum for any question)";

    }
    val=arg->getAttribute("applyIncrementFactor",nullptr) ;
    if(val)
    {
        msg_deprecated() << "The attribute 'applyIncrementFactor' is deprecated since Sofa 19.06'" << msgendl
                         << "The incremental loading is now supposed to be done within the desired ForceField." << msgendl
                         << "(Use the forum for any question)";

    }
    sofa::core::behavior::OdeSolver::parse(arg) ;
}

void StaticSolver::solve(const sofa::core::ExecParams* params, double dt, sofa::core::MultiVecCoordId xResult, sofa::core::MultiVecDerivId /*vResult*/) {

    SOFA_UNUSED(dt);

    sofa::simulation::common::VectorOperations vop( params, this->getContext() );
    sofa::simulation::common::MechanicalOperations mop( params, this->getContext() );

    MultiVecCoord x_start(&vop, sofa::core::VecCoordId::position() );
    MultiVecCoord x(&vop, xResult );
    MultiVecDeriv force( &vop, sofa::core::VecDerivId::force() );
    dx.realloc( &vop, true );

    // MO vector dx is not allocated by default, it will seg fault if the CG is used (dx is taken by default) with an IdentityMapping
    MultiVecDeriv tempdx(&vop, sofa::core::VecDerivId::dx() ); tempdx.realloc( &vop, true, true );

    // Set implicit param to true to trigger nonlinear stiffness matrix recomputation
    mop->setImplicit(true);

    msg_info() << "======= Starting static ODE solver in time step " << this->getTime();
    msg_info() << "(doing a maximum of " << d_newton_iterations.getValue() << " newton iterations)";

    unsigned n_it=0;
    double dx_norm = -1.0, f_norm;

    sofa::helper::AdvancedTimer::stepBegin("StaticSolver::Solve");

    // compute addForce, in mapped: addForce + applyJT (vec)
    // Initial computation
    force.clear();
    mop.computeForce(force);
    mop.projectResponse(force);
    f_norm = sqrt(force.dot(force));

    if (d_residual_tolerance_threshold.getValue() > 0 && f_norm <= d_residual_tolerance_threshold.getValue()) {
        msg_info() << "The ODE has already reached an equilibrium state";
    } else {

        while (n_it < d_newton_iterations.getValue()) {
            std::string stepname = "step_" + std::to_string(n_it);
            sofa::helper::AdvancedTimer::stepBegin(stepname.c_str());


            // Assemble matrix, CG: does nothing
            // LDL non-mapped: addKToMatrix added to system matrix
            // LDL mapped: addKToMatrix not added to the system matrix, needs mapped FF (TODO)
            sofa::core::behavior::MultiMatrix<sofa::simulation::common::MechanicalOperations> matrix(&mop);
            matrix = MechanicalMatrix::K * -1.0;

            // for CG: calls iteratively addDForce, mapped:  [applyJ, addDForce, applyJt(vec)]+
            // for LDL: solves the system, everything's already assembled
            matrix.solve(dx, force);

            x.eq(x_start, dx, 1);
            mop.solveConstraint(x, sofa::core::ConstraintParams::POS);

            // Propagate positions to mapped nodes: taken from AnimateVisitor::processNodeTopDown executed by the animation loop
            // calls apply, applyJ
            sofa::core::MechanicalParams mp;
            sofa::simulation::MechanicalPropagateOnlyPositionAndVelocityVisitor(&mp).execute(
                this->getContext()); // propagate the changes to mappings below

            // Compute addForce, in mapped: addForce + applyJT (vec)
            force.clear();
            mop.computeForce(force);
            mop.projectResponse(force);
            double f_cur_norm = sqrt(force.dot(force));

            dx_norm = sqrt(dx.dot(dx));


            msg_info() << "Newton iteration #" << n_it << ": |f - K(x0 + dx)| = " << f_cur_norm << " |dx| = " << dx_norm;
            sofa::helper::AdvancedTimer::valSet("residual", f_cur_norm);
            sofa::helper::AdvancedTimer::valSet("correction", dx_norm);
            sofa::helper::AdvancedTimer::stepEnd(stepname.c_str());

            if (d_should_diverge_when_residual_is_growing.getValue() && f_cur_norm > f_norm && n_it>1) {
                msg_info() << "[DIVERGED] residual's norm increased";
                break;
            }

            if (dx_norm <= this->d_correction_tolerance_threshold.getValue()) {
                msg_info() << "[CONVERGED] The correction's norm |dx| is smaller than the threshold of "
                           << d_correction_tolerance_threshold;
                break;
            }

            if (d_residual_tolerance_threshold.getValue() > 0 && f_cur_norm <= d_residual_tolerance_threshold.getValue()) {
                msg_info() << "[CONVERGED] The residual's norm |f - K(x0 + dx)| is smaller than the threshold of "
                           << d_residual_tolerance_threshold;
                break;
            }

            f_norm = f_cur_norm;
            n_it++;

        } // End while (n_it < d_newton_iterations.getValue())
    }

    if (n_it >= d_newton_iterations.getValue()) {
        n_it--;
        msg_info() << "[DIVERGED] The number of Newton iterations reached the threshold of " << d_newton_iterations << " iterations";
    }

    sofa::helper::AdvancedTimer::valSet("nb_iterations", n_it+1);
    sofa::helper::AdvancedTimer::valSet("residual", f_norm);
    sofa::helper::AdvancedTimer::valSet("correction", dx_norm);
    sofa::helper::AdvancedTimer::stepEnd("StaticSolver::Solve");
}


int StaticSolverClass = sofa::core::RegisterObject("Static ODE Solver")
    .add< StaticSolver >()
;

} // namespace odesolver

} // namespace component

} // namespace sofa

