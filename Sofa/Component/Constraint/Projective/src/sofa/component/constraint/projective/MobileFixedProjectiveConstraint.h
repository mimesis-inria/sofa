/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or at      *
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

#pragma once

#include <sofa/component/constraint/projective/config.h>

#include <sofa/core/behavior/ProjectiveConstraintSet.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/topology/TopologySubsetIndices.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>

#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/linearalgebra/BaseMatrix.h>
#include <sofa/linearalgebra/BaseVector.h>

#include <sofa/type/Quat.h>
#include <sofa/type/Vec.h>
#include <sofa/type/vector.h>
#include <sofa/type/BoundingBox.h>

namespace sofa::component::constraint::projective
{

class SOFA_COMPONENT_CONSTRAINT_PROJECTIVE_API MobileFixedProjectiveConstraint
    : public sofa::core::behavior::ProjectiveConstraintSet<sofa::defaulttype::Rigid3Types>
{
public:
    SOFA_CLASS(MobileFixedProjectiveConstraint,
               SOFA_TEMPLATE(sofa::core::behavior::ProjectiveConstraintSet, sofa::defaulttype::Rigid3Types));

    using DataTypes       = sofa::defaulttype::Rigid3Types;
    using Real            = DataTypes::Real;
    using Coord           = DataTypes::Coord;
    using Deriv           = DataTypes::Deriv;
    using VecCoord        = DataTypes::VecCoord;
    using VecDeriv        = DataTypes::VecDeriv;
    using MatrixDeriv     = DataTypes::MatrixDeriv;
    using Index           = sofa::Index;

    using DataVecCoord    = Data<VecCoord>;
    using DataVecDeriv    = Data<VecDeriv>;
    using DataMatrixDeriv = Data<MatrixDeriv>;

    using SetIndexArray   = sofa::type::vector<Index>;
    using SetIndex        = sofa::core::topology::TopologySubsetIndices;

    using Vec3            = sofa::type::Vec<3, Real>;
    using Quat            = sofa::type::Quat<Real>;
    using BoolVector      = sofa::type::vector<bool>;

public:
    MobileFixedProjectiveConstraint();
    ~MobileFixedProjectiveConstraint() override = default;

    void init() override;
    void reinit() override;

    void projectResponse(const sofa::core::MechanicalParams* mparams, DataVecDeriv& resData) override;
    void projectVelocity(const sofa::core::MechanicalParams* mparams, DataVecDeriv& vData) override;
    void projectPosition(const sofa::core::MechanicalParams* mparams, DataVecCoord& xData) override;
    void projectJacobianMatrix(const sofa::core::MechanicalParams* mparams, DataMatrixDeriv& cData) override;

    void applyConstraint(const sofa::core::MechanicalParams* mparams,
                         const sofa::core::behavior::MultiMatrixAccessor* matrix) override;

    void applyConstraint(const sofa::core::MechanicalParams* mparams,
                         sofa::linearalgebra::BaseVector* vect,
                         const sofa::core::behavior::MultiMatrixAccessor* matrix) override;

    void applyConstraint(sofa::core::behavior::ZeroDirichletCondition* matrix) override;
    void projectMatrix(sofa::linearalgebra::BaseMatrix* M, unsigned offset) override;

    void computeBBox(const sofa::core::ExecParams* params, bool onlyVisible) override;
    void draw(const sofa::core::visual::VisualParams* vparams) override;

    void clearConstraints();
    void addConstraint(Index index);
    void removeConstraint(Index index);

protected:
    void checkIndices();
    void checkMask();
    void captureInitialCoordinates();

    bool isDofConstrained(unsigned int localDof) const;
    bool hasRHSForce() const;

    unsigned int getIterationFromTime() const;
    Real getElapsedTime() const;

    Coord computeTargetCoord(const Coord& initialCoord) const;
    void setConstraintOnMatrixRowCol(sofa::linearalgebra::BaseMatrix* matrix, unsigned int globalDofIndex) const;
    void computeBBoxFromIndices(const SetIndexArray& indices);

protected:
    SetIndex d_indices;

    Data<bool>  d_fixAll;
    Data<bool>  d_projectVelocity;
    Data<bool>  d_relativeMovement;
    Data<bool>  d_showObject;
    Data<SReal> d_drawSize;

    // Dynamic prescribed motion
    Data<Vec3> d_linearVelocity;     // translation velocity [vx vy vz]
    Data<Vec3> d_angularVelocity;    // angular velocity [wx wy wz] in rad/s
    Data<Real> d_startTime;

    // Optional RHS-force mode.
    // When useRHSForce=false, the component behaves as a usual projective constraint.
    // When useRHSForce=true, the selected rigid DoFs are not projected out of the matrix;
    // rhsForce is added to their response/RHS vector instead.
    Data<bool>  d_useRHSForce;
    Data<Deriv> d_rhsForce;          // generalized force [Fx Fy Fz Mx My Mz]

    // Optional per-rigid-DoF mask [tx ty tz rx ry rz]
    Data<BoolVector> d_constrainedDofMask;

    SingleLink<MobileFixedProjectiveConstraint,
               sofa::core::topology::BaseMeshTopology,
               BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_topology;

    VecCoord m_initialCoordinates;
};

} // namespace sofa::component::constraint::projective