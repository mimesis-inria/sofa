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
#define SOFA_COMPONENT_PROJECTIVECONSTRAINTSET_MOBILEFIXEDPROJECTIVECONSTRAINT_CPP

#include <sofa/component/constraint/projective/MobileFixedProjectiveConstraint.h>

#include <sofa/core/behavior/MultiMatrixAccessor.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/helper/accessor.h>
#include <sofa/type/vector_algorithm.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace sofa::component::constraint::projective
{

using sofa::core::objectmodel::ComponentState;
using sofa::defaulttype::Rigid3Types;

MobileFixedProjectiveConstraint::MobileFixedProjectiveConstraint()
    : sofa::core::behavior::ProjectiveConstraintSet<Rigid3Types>(nullptr)
    , d_indices(initData(&d_indices, "indices", "Indices of constrained rigid DoFs"))
    , d_fixAll(initData(&d_fixAll, false, "fixAll", "If true, constrain all rigid DoFs"))
    , d_projectVelocity(initData(&d_projectVelocity, true, "projectVelocity", "If true, zero constrained velocities"))
    , d_relativeMovement(initData(&d_relativeMovement, true, "relativeMovement", "If true, motion is applied relative to initial pose"))
    , d_showObject(initData(&d_showObject, true, "showObject", "Show constrained rigid points"))
    , d_drawSize(initData(&d_drawSize, (SReal)0.0, "drawSize", "Drawing size"))
    , d_linearVelocity(initData(&d_linearVelocity, Vec3(0, 0, 0), "linearVelocity", "Prescribed translational velocity [vx vy vz]"))
    , d_angularVelocity(initData(&d_angularVelocity, Vec3(0, 0, 0), "angularVelocity", "Prescribed angular velocity [wx wy wz] in rad/s"))
    , d_startTime(initData(&d_startTime, (Real)0, "startTime", "Motion start time"))
    , d_useRHSForce(initData(&d_useRHSForce, false, "useRHSForce", "If true, selected rigid DoFs are force-driven instead of projectively fixed"))
    , d_rhsForce(initData(&d_rhsForce, Deriv(), "rhsForce", "Generalized RHS force [Fx Fy Fz Mx My Mz] added to selected rigid DoFs when useRHSForce is true"))
    , d_constrainedDofMask(initData(&d_constrainedDofMask, BoolVector{true, true, true, true, true, true},
                                    "constrainedDofMask", "Mask [tx ty tz rx ry rz]"))
    , l_topology(initLink("topology", "Link to the topology container"))
{
    this->addUpdateCallback("updateIndices", {&d_indices}, [this](const core::DataTracker&)
    {
        checkIndices();
        return ComponentState::Valid;
    }, {});

    this->addUpdateCallback("updateMask", {&d_constrainedDofMask}, [this](const core::DataTracker&)
    {
        checkMask();
        return ComponentState::Valid;
    }, {});
}

void MobileFixedProjectiveConstraint::init()
{
    this->d_componentState.setValue(ComponentState::Invalid);
    this->core::behavior::ProjectiveConstraintSet<Rigid3Types>::init();

    if (!this->mstate.get())
    {
        msg_warning() << "Missing MechanicalState.";
        return;
    }

    if (l_topology.empty())
    {
        l_topology.set(this->getContext()->getMeshTopologyLink());
    }

    if (auto* topology = l_topology.get())
    {
        d_indices.createTopologyHandler(topology);
    }

    checkIndices();
    checkMask();
    captureInitialCoordinates();

    this->d_componentState.setValue(ComponentState::Valid);
}

void MobileFixedProjectiveConstraint::reinit()
{
    checkIndices();
    checkMask();
    captureInitialCoordinates();
}

void MobileFixedProjectiveConstraint::clearConstraints()
{
    auto* indices = d_indices.beginEdit();
    indices->clear();
    d_indices.endEdit();
}

void MobileFixedProjectiveConstraint::addConstraint(Index index)
{
    auto* indices = d_indices.beginEdit();
    indices->push_back(index);
    d_indices.endEdit();
}

void MobileFixedProjectiveConstraint::removeConstraint(Index index)
{
    auto* indices = d_indices.beginEdit();
    sofa::type::removeValue(*indices, index);
    d_indices.endEdit();
}

void MobileFixedProjectiveConstraint::checkIndices()
{
    if (!this->mstate.get())
        return;

    const Index maxIndex = this->mstate->getSize();
    const auto& indices = d_indices.getValue();

    SetIndexArray invalid;
    for (Index idx : indices)
    {
        if (idx >= maxIndex)
            invalid.push_back(idx);
    }

    for (Index idx : invalid)
    {
        msg_warning() << "Removing invalid index " << idx << ", valid range is [0," << maxIndex << ").";
        removeConstraint(idx);
    }
}

void MobileFixedProjectiveConstraint::checkMask()
{
    auto mask = d_constrainedDofMask.getValue();
    if (mask.size() != 6)
    {
        msg_warning() << "constrainedDofMask must have 6 entries [tx ty tz rx ry rz]. Resetting to all true.";
        d_constrainedDofMask.setValue(BoolVector{true, true, true, true, true, true});
    }
}

void MobileFixedProjectiveConstraint::captureInitialCoordinates()
{
    if (!this->mstate.get())
        return;

    m_initialCoordinates = this->mstate->read(core::vec_id::read_access::position)->getValue();
}

bool MobileFixedProjectiveConstraint::isDofConstrained(unsigned int localDof) const
{
    const auto& mask = d_constrainedDofMask.getValue();
    return (localDof < mask.size()) ? mask[localDof] : true;
}

bool MobileFixedProjectiveConstraint::hasRHSForce() const
{
    return d_useRHSForce.getValue();
}

unsigned int MobileFixedProjectiveConstraint::getIterationFromTime() const
{
    const SReal dt = this->getContext()->getDt();
    if (dt <= std::numeric_limits<SReal>::epsilon())
        return 0u;

    const SReal time = this->getContext()->getTime();
    return static_cast<unsigned int>(std::floor(time / dt + (SReal)1e-9));
}

MobileFixedProjectiveConstraint::Real MobileFixedProjectiveConstraint::getElapsedTime() const
{
    const Real time = static_cast<Real>(this->getContext()->getTime());
    const Real t0   = d_startTime.getValue();
    return std::max((Real)0, time - t0);
}

MobileFixedProjectiveConstraint::Coord
MobileFixedProjectiveConstraint::computeTargetCoord(const Coord& initialCoord) const
{
    Coord target = initialCoord;

    const Real t = getElapsedTime();
    const Vec3 translation = d_linearVelocity.getValue() * t;
    const Vec3 angles = d_angularVelocity.getValue() * t;

    Quat qDelta = Quat::createQuaterFromEuler(angles);
    qDelta.normalize();

    if (d_relativeMovement.getValue())
    {
        target.getCenter() = initialCoord.getCenter() + translation;
        target.getOrientation() = initialCoord.getOrientation() * qDelta;
        target.getOrientation().normalize();
    }
    else
    {
        target.getCenter() = translation;
        target.getOrientation() = qDelta;
        target.getOrientation().normalize();
    }

    return target;
}

void MobileFixedProjectiveConstraint::projectResponse(const core::MechanicalParams* mparams, DataVecDeriv& resData)
{
    SOFA_UNUSED(mparams);
    helper::WriteAccessor<DataVecDeriv> response(resData);

    msg_info() << "Called ProjectResponse.";

    auto clearConstrained = [this](Deriv& d)
    {
        for (unsigned int i = 0; i < Deriv::size(); ++i)
        {
            if (isDofConstrained(i))
                d[i] = (Real)0;
        }
    };


    if (hasRHSForce())
    {
        const Deriv& rhsForce = d_rhsForce.getValue();

        auto applyRHSForce = [this, &rhsForce, &clearConstrained](Deriv& d)
        {
            // First recover the normal projective behavior.
            clearConstrained(d);

            // Then prescribe a nonzero RHS value on the constrained entries.
            // If rhsForce = 0, this is exactly the original constrained behavior.
            for (unsigned int i = 0; i < Deriv::size(); ++i)
            {
                if (isDofConstrained(i))
                    d[i] += rhsForce[i];
            }
        };

        if (d_fixAll.getValue())
        {
            for (auto& r : *response)
                applyRHSForce(r);
        }
        else
        {
            for (Index idx : d_indices.getValue())
            {
                if (idx < response.size())
                    applyRHSForce((*response)[idx]);
                else
                    msg_warning() << "Ignoring rhsForce for invalid index " << idx << ".";
            }
        }

        return;
    }

    if (d_fixAll.getValue())
    {
        for (auto& r : *response)
            clearConstrained(r);
    }
    else
    {
        for (Index idx : d_indices.getValue())
        {
            if (idx < response.size())
                clearConstrained((*response)[idx]);
            else
                msg_warning() << "Ignoring projectResponse for invalid index " << idx << ".";
        }
    }
}

void MobileFixedProjectiveConstraint::projectVelocity(const core::MechanicalParams* mparams, DataVecDeriv& vData)
{
    SOFA_UNUSED(mparams);

    msg_info() << "Called projectVelocity.";

    if (hasRHSForce())
        return;

    if (!d_projectVelocity.getValue())
        return;

    helper::WriteAccessor<DataVecDeriv> velocity(vData);

    auto clearConstrained = [this](Deriv& d)
    {
        for (unsigned int i = 0; i < Deriv::size(); ++i)
        {
            if (isDofConstrained(i))
                d[i] = (Real)0;
        }
    };

    if (d_fixAll.getValue())
    {
        for (auto& v : *velocity)
            clearConstrained(v);
    }
    else
    {
        for (Index idx : d_indices.getValue())
        {
            if (idx < velocity.size())
                clearConstrained((*velocity)[idx]);
            else
                msg_warning() << "Ignoring projectVelocity for invalid index " << idx << ".";
        }
    }
}

void MobileFixedProjectiveConstraint::projectPosition(const core::MechanicalParams* mparams, DataVecCoord& xData)
{
    SOFA_UNUSED(mparams);

    msg_info() << "Called projectPosition.";

    if (hasRHSForce())
        return;

    helper::WriteAccessor<DataVecCoord> position(xData);

    auto projectOne = [this](Coord& current, const Coord& initial)
    {
        const Coord target = computeTargetCoord(initial);

        const bool lockT = isDofConstrained(0) && isDofConstrained(1) && isDofConstrained(2);
        const bool lockR = isDofConstrained(3) && isDofConstrained(4) && isDofConstrained(5);

        if (lockT)
            current.getCenter() = target.getCenter();

        if (lockR)
            current.getOrientation() = target.getOrientation();
    };

    if (d_fixAll.getValue())
    {
        const Size n = std::min(position.size(), m_initialCoordinates.size());
        for (Size i = 0; i < n; ++i)
            projectOne((*position)[i], m_initialCoordinates[i]);
    }
    else
    {
        for (Index idx : d_indices.getValue())
        {
            if (idx < position.size() && idx < m_initialCoordinates.size())
                projectOne((*position)[idx], m_initialCoordinates[idx]);
            else
                msg_warning() << "Ignoring projectPosition for invalid index " << idx << ".";
        }
    }
}

void MobileFixedProjectiveConstraint::projectJacobianMatrix(const core::MechanicalParams* mparams, DataMatrixDeriv& cData)
{
    SOFA_UNUSED(mparams);
    helper::WriteAccessor<DataMatrixDeriv> jacobian(cData);

    msg_info() << "Called projectJacobianMatrix.";

    // if (hasRHSForce())
    //     return;

    if (d_fixAll.getValue())
    {
        jacobian->clear();
        return;
    }

    for (Index idx : d_indices.getValue())
        jacobian->clearColBlock(idx);
}

void MobileFixedProjectiveConstraint::setConstraintOnMatrixRowCol(sofa::linearalgebra::BaseMatrix* matrix,
                                                                  unsigned int globalDofIndex) const
{
    matrix->clearRowCol(globalDofIndex);
    matrix->set(globalDofIndex, globalDofIndex, 1.0);
}

void MobileFixedProjectiveConstraint::applyConstraint(const core::MechanicalParams* mparams,
                                                      const core::behavior::MultiMatrixAccessor* matrix)
{
    SOFA_UNUSED(mparams);

    msg_info() << "Called applyConstraint wo linear algebra.";

    // if (hasRHSForce())
    //     return;

    if (const auto r = matrix->getMatrix(this->mstate.get()))
    {
        static constexpr unsigned int blockSize = Deriv::size();

        auto constrainRigid = [this, &r](unsigned int rigidIndex)
        {
            for (unsigned int c = 0; c < blockSize; ++c)
            {
                if (isDofConstrained(c))
                    setConstraintOnMatrixRowCol(r.matrix, r.offset + blockSize * rigidIndex + c);
            }
        };

        if (d_fixAll.getValue())
        {
            const unsigned int n = static_cast<unsigned int>(this->mstate->getSize());
            for (unsigned int i = 0; i < n; ++i)
                constrainRigid(i);
        }
        else
        {
            for (Index idx : d_indices.getValue())
                constrainRigid(static_cast<unsigned int>(idx));
        }
    }
}

void MobileFixedProjectiveConstraint::applyConstraint(const core::MechanicalParams* mparams,
                                                      sofa::linearalgebra::BaseVector* vect,
                                                      const core::behavior::MultiMatrixAccessor* matrix)
{
    msg_info() << "Called applyConstraint.";

    SOFA_UNUSED(mparams);

    // if (hasRHSForce())
    //     return;

    const int globalOffset = matrix->getGlobalOffset(this->mstate.get());
    if (globalOffset < 0)
        return;

    static constexpr unsigned int blockSize = Deriv::size();
    const unsigned int offset = static_cast<unsigned int>(globalOffset);

    auto clearRigid = [this, vect, offset](unsigned int rigidIndex)
    {
        for (unsigned int c = 0; c < blockSize; ++c)
        {
            if (isDofConstrained(c))
                vect->clear(offset + blockSize * rigidIndex + c);
        }
    };

    if (d_fixAll.getValue())
    {
        const unsigned int n = static_cast<unsigned int>(this->mstate->getSize());
        for (unsigned int i = 0; i < n; ++i)
            clearRigid(i);
    }
    else
    {
        for (Index idx : d_indices.getValue())
            clearRigid(static_cast<unsigned int>(idx));
    }
}

void MobileFixedProjectiveConstraint::applyConstraint(sofa::core::behavior::ZeroDirichletCondition* matrix)
{
    static constexpr unsigned int blockSize = Deriv::size();

    msg_info() << "Called applyConstraint zero drilecht.";

    // if (hasRHSForce())
    //     return;

    auto discardRigid = [this, matrix](unsigned int rigidIndex)
    {
        for (unsigned int c = 0; c < blockSize; ++c)
        {
            if (isDofConstrained(c))
                matrix->discardRowCol(blockSize * rigidIndex + c, blockSize * rigidIndex + c);
        }
    };

    if (d_fixAll.getValue())
    {
        const unsigned int n = static_cast<unsigned int>(this->mstate->getSize());
        for (unsigned int i = 0; i < n; ++i)
            discardRigid(i);
    }
    else
    {
        for (Index idx : d_indices.getValue())
            discardRigid(static_cast<unsigned int>(idx));
    }
}

void MobileFixedProjectiveConstraint::projectMatrix(sofa::linearalgebra::BaseMatrix* M, unsigned offset)
{
    msg_info() << "Called projectMatrix.";

    // if (hasRHSForce())
    //     return;

    static constexpr unsigned int blockSize = Deriv::size();

    auto clearRigid = [this, M, offset](unsigned int rigidIndex)
    {
        for (unsigned int c = 0; c < blockSize; ++c)
        {
            if (isDofConstrained(c))
                M->clearRowsCols(offset + rigidIndex * blockSize + c,
                                 offset + rigidIndex * blockSize + c + 1);
        }
    };

    if (d_fixAll.getValue())
    {
        const unsigned int n = static_cast<unsigned int>(this->mstate->getSize());
        for (unsigned int i = 0; i < n; ++i)
            clearRigid(i);
    }
    else
    {
        for (Index idx : d_indices.getValue())
            clearRigid(static_cast<unsigned int>(idx));
    }
}

void MobileFixedProjectiveConstraint::computeBBoxFromIndices(const SetIndexArray& indices)
{
    const auto drawSize = static_cast<Real>(d_drawSize.getValue());
    const VecCoord& x = this->mstate->read(core::vec_id::read_access::position)->getValue();

    sofa::type::BoundingBox bbox;
    for (Index idx : indices)
    {
        const auto p = DataTypes::getCPos(x[idx]);
        for (unsigned int i = 0; i < 3; ++i)
        {
            bbox.maxBBox()[i] = std::max(static_cast<SReal>(p[i] + drawSize), bbox.maxBBox()[i]);
            bbox.minBBox()[i] = std::min(static_cast<SReal>(p[i] - drawSize), bbox.minBBox()[i]);
        }
    }
    this->f_bbox.setValue(bbox);
}

void MobileFixedProjectiveConstraint::computeBBox(const core::ExecParams* params, bool onlyVisible)
{
    SOFA_UNUSED(params);

    if (onlyVisible && !d_showObject.getValue())
        return;
    if (this->d_componentState.getValue() == ComponentState::Invalid)
        return;

    if (d_fixAll.getValue())
    {
        this->f_bbox.setValue(this->mstate->computeBBox());
    }
    else if (!d_indices.getValue().empty())
    {
        computeBBoxFromIndices(d_indices.getValue());
    }
}

void MobileFixedProjectiveConstraint::draw(const core::visual::VisualParams* vparams)
{
    if (this->d_componentState.getValue() != ComponentState::Valid) return;
    if (!vparams->displayFlags().getShowBehaviorModels()) return;
    if (!d_showObject.getValue()) return;
    if (!this->isActive()) return;

    const auto stateLifeCycle = vparams->drawTool()->makeStateLifeCycle();
    SOFA_UNUSED(stateLifeCycle);

    const VecCoord& x = this->mstate->read(core::vec_id::read_access::position)->getValue();

    std::vector<sofa::type::Vec3> points;
    if (d_fixAll.getValue())
    {
        for (unsigned int i = 0; i < x.size(); ++i)
            points.push_back(sofa::type::toVec3(DataTypes::getCPos(x[i])));
    }
    else
    {
        for (Index idx : d_indices.getValue())
            points.push_back(sofa::type::toVec3(DataTypes::getCPos(x[idx])));
    }

    if (d_drawSize.getValue() == 0)
    {
        vparams->drawTool()->drawPoints(points, 10.0f, sofa::type::RGBAColor(1, 0.5f, 0.5f, 1));
    }
    else
    {
        vparams->drawTool()->setLightingEnabled(true);
        vparams->drawTool()->drawSpheres(points, (float)d_drawSize.getValue(),
                                         sofa::type::RGBAColor(1.0f, 0.35f, 0.35f, 1.0f));
    }
}

void registerMobileFixedProjectiveConstraint(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(core::ObjectRegistrationData(
    "Dynamic rigid prescribed-motion projective constraint without keyTimes/finished state.")
    .add<MobileFixedProjectiveConstraint>());
}

} // namespace sofa::component::constraint::projective