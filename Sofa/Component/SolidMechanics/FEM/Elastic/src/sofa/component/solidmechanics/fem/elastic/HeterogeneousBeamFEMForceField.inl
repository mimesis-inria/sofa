/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
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
#pragma once
#include <sofa/component/solidmechanics/fem/elastic/HeterogeneousBeamFEMForceField.h>
#include <sofa/component/solidmechanics/fem/elastic/BaseLinearElasticityFEMForceField.inl>
#include <sofa/core/topology/TopologyData.inl>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/helper/rmath.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/behavior/MultiMatrixAccessor.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/behavior/BaseLocalForceFieldMatrix.h>

namespace sofa::component::solidmechanics::fem::elastic::_heterogeneousbeamfemforcefield_
{

template<class DataTypes>
HeterogeneousBeamFEMForceField<DataTypes>::HeterogeneousBeamFEMForceField()
    : Inherit()
    , d_sectionStart(this->initData(&d_sectionStart, "sectionStart",
        "Inclusive start edge index for each section."))
    , d_sectionEnd(this->initData(&d_sectionEnd, "sectionEnd",
        "Exclusive end edge index for each section."))
    , d_sectionYoung(this->initData(&d_sectionYoung, "sectionYoung",
        "Young modulus E for each section."))
    , d_sectionPoisson(this->initData(&d_sectionPoisson, "sectionPoisson",
        "Poisson ratio nu for each section."))
    , d_sectionRadius(this->initData(&d_sectionRadius, "sectionRadius",
        "Outer radius for each section."))
    , d_sectionInnerRadius(this->initData(&d_sectionInnerRadius, "sectionInnerRadius",
        "Inner radius for each section (0 for solid)."))
    , d_edgeToSection(this->initData(&d_edgeToSection, "Edge to Section Map",
        "Map between edges and sections (-1 if no section found)"))
    , l_sectionYoungLinks(this->initLink(
      "sectionYoungLinks",
      "Optional list of linked parameter components, one per section, used instead of sectionYoung."))
{
    this->d_radius.setRequired(false);
}

template<class DataTypes>
void HeterogeneousBeamFEMForceField<DataTypes>::init()
{
    Inherit::init();
}

template<class DataTypes>
bool HeterogeneousBeamFEMForceField<DataTypes>::validateSections(Index nbEdges) const
{
    const auto& start = d_sectionStart.getValue();
    const auto& end   = d_sectionEnd.getValue();

    if (start.size() != end.size())
        return false;

    const size_t nSec = start.size();
    if (nSec == 0)
        return true;

    auto okSize = [&](const auto& v) { return v.size() == nSec; };

    const bool useLinkedYoung = !l_sectionYoungLinks.empty();

    if (!useLinkedYoung)
    {
        if (!okSize(d_sectionYoung.getValue())) return false;
    }
    else
    {
        if (l_sectionYoungLinks.size() != nSec) return false;
    }
    if (!okSize(d_sectionPoisson.getValue()))     return false;
    if (!okSize(d_sectionRadius.getValue()))      return false;
    if (!okSize(d_sectionInnerRadius.getValue())) return false;

    for (size_t s = 0; s < nSec; ++s)
    {
        if (start[s] > end[s]) return false;
        if (end[s] > nbEdges)  return false;
    }
    return true;
}

template<class DataTypes>
void HeterogeneousBeamFEMForceField<DataTypes>::buildEdgeToSectionMap()
{
    const Index nbEdges = static_cast<Index>(this->m_indexedElements->size());
    m_edgeToSection.assign(nbEdges, -1);

    const auto& start = d_sectionStart.getValue();
    const auto& end   = d_sectionEnd.getValue();

    for (size_t s = 0; s < start.size(); ++s)
        for (Index e = start[s]; e < end[s]; ++e)
            m_edgeToSection[e] = static_cast<int>(s);
}

template<class DataTypes>
int HeterogeneousBeamFEMForceField<DataTypes>::sectionForEdge(Index edgeIndex) const
{
    if (edgeIndex >= m_edgeToSection.size())
        return -1;
    return m_edgeToSection[edgeIndex];
}

template<class DataTypes>
void HeterogeneousBeamFEMForceField<DataTypes>::reinit()
{
    if (!this->m_indexedElements)
    {
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }

    const Index nbEdges = static_cast<Index>(this->m_indexedElements->size());

    if (!validateSections(nbEdges))
    {
        msg_error() << "Invalid section definition: check sectionStart/sectionEnd and parameter vector sizes.";
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }
    
    buildEdgeToSectionMap();
    d_edgeToSection.setValue(m_edgeToSection);
    m_resolvedSectionYoung.clear();

    if (!l_sectionYoungLinks.empty())
    {
        const size_t nSec = d_sectionStart.getValue().size();
        m_resolvedSectionYoung.resize(nSec);

        for (size_t s = 0; s < nSec; ++s)
        {
            auto* obj = l_sectionYoungLinks[s];
            if (!obj)
            {
                msg_error() << "sectionYoungLinks[" << s << "] is null.";
                this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
                return;
            }

            // Replace this block with the actual type/data accessor of your parameter component.
            // Example strategy: read its "value" Data and extract one scalar.
            auto* valueData = obj->findData("value");
            if (!valueData)
            {
                msg_error() << "Linked object '" << obj->getName()
                            << "' has no Data named 'value'.";
                this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
                return;
            }

            std::stringstream ss(valueData->getValueString());
            SReal E;
            ss >> E;
            if (ss.fail())
            {
                msg_error() << "Could not parse scalar Young modulus from linked object '"
                            << obj->getName() << "'.";
                this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
                return;
            }

            m_resolvedSectionYoung[s] = E;
        }
    }
    else
    {
        m_resolvedSectionYoung = d_sectionYoung.getValue();
    }

    Inherit::reinit();
}

template<class DataTypes>
void HeterogeneousBeamFEMForceField<DataTypes>::reinitBeam(Index i)
{
    // If no sections or edge uncovered, fall back to base behavior
    if (d_sectionStart.getValue().empty())
    {
        Inherit::reinitBeam(i);
        return;
    }

    const int s = sectionForEdge(i);
    if (s < 0)
    {
        Inherit::reinitBeam(i);
        return;
    }

    // Edge nodes
    const auto& [a, b] = (*this->m_indexedElements)[i].array();

    // Rest positions -> length (same style as BeamFEMForceField)
    const VecCoord& x0 = this->mstate->read(core::vec_id::read_access::restPosition)->getValue();
    const SReal length = (x0[a].getCenter() - x0[b].getCenter()).norm();

    const auto& R  = d_sectionRadius.getValue();
    const auto& Ri = d_sectionInnerRadius.getValue();
    const auto& Ey = m_resolvedSectionYoung.empty() ? d_sectionYoung.getValue()
                                                : m_resolvedSectionYoung;
    const auto& Poissony = d_sectionPoisson.getValue();

    // assumes validateSections() ensured these sizes:
    SReal E      = Ey[s];
    SReal nu     = Poissony[s] ;

    SReal radius = R[s];
    SReal inner  = Ri[s];
    // Same pipeline as base class
    this->setBeam(i, E, length, nu, radius, inner);
    this->computeStiffness(i, a, b);
    this->initLarge(i, a, b);
}

inline sofa::type::Quat<SReal> qDiff(sofa::type::Quat<SReal> a, const sofa::type::Quat<SReal>& b)
{
    if (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3] < 0)
    {
        a[0] = -a[0]; a[1] = -a[1]; a[2] = -a[2]; a[3] = -a[3];
    }
    const sofa::type::Quat<SReal> q = b.inverse() * a;
    return q;
}


template<class DataTypes>
void HeterogeneousBeamFEMForceField<DataTypes>::LocalForceLarge( const VecCoord & x, int i, Index a, Index b, Deriv& fa, Deriv& fb)
{
    const VecCoord& x0 = this->mstate->read(core::vec_id::read_access::restPosition)->getValue();

    this->beamQuat(i)= x[a].getOrientation();
    this->beamQuat(i).normalize();

    type::Vec<3,Real> u, P1P2, P1P2_0;

    // local displacement
    Displacement depl;

    // translations //
    P1P2_0 = x0[b].getCenter() - x0[a].getCenter();
    P1P2_0 = x0[a].getOrientation().inverseRotate(P1P2_0);
    P1P2 = x[b].getCenter() - x[a].getCenter();
    P1P2 = x[a].getOrientation().inverseRotate(P1P2);
    u = P1P2 - P1P2_0;

    depl[0] = 0.0; 	depl[1] = 0.0; 	depl[2] = 0.0;
    depl[6] = u[0]; depl[7] = u[1]; depl[8] = u[2];

    // rotations //
    type::Quat<SReal> dQ0, dQ;

    dQ0 = qDiff(x0[b].getOrientation(), x0[a].getOrientation());
    dQ =  qDiff(x[b].getOrientation(), x[a].getOrientation());

    dQ0.normalize();
    dQ.normalize();

    type::Quat<SReal> tmpQ = qDiff(dQ,dQ0);
    tmpQ.normalize();

    // TODO(e.coevoet) remove before v20.12
    // Use of quatToRotationVector instead of toEulerVector: u = tmpQ.quatToRotationVector();
    // this is done to keep the old behavior (before the
    // correction of the toEulerVector  function). If the
    // purpose was to obtain the Eulerian vector and not the
    // rotation vector please use the following line instead
    // u = tmpQ.toEulerVector();
    u = tmpQ.quatToRotationVector();

    depl[3] = 0.0; 	depl[4] = 0.0; 	depl[5] = 0.0;
    depl[9] = u[0]; depl[10]= u[1]; depl[11]= u[2];

    // this computation can be optimised: (we know that half of "depl" is null)
    Displacement force = this->d_beamsData.getValue()[i]._k_loc * depl;


    // Apply lambda transpose (we use the rotation value of point a for the beam)
    const Vec3 fa1 = x[a].getOrientation().rotate(type::Vec3d(force[0],force[1],force[2]));
    const Vec3 fa2 = x[a].getOrientation().rotate(type::Vec3d(force[3],force[4],force[5]));

    const Vec3 fb1 = x[a].getOrientation().rotate(type::Vec3d(force[6],force[7],force[8]));
    const Vec3 fb2 = x[a].getOrientation().rotate(type::Vec3d(force[9],force[10],force[11]));

    fa = Deriv(-fa1, -fa2);
    fb = Deriv(-fb1, -fb2);

}

} // namespace sofa::component::solidmechanics::fem::elastic::_heterogeneousbeamfemforcefield_