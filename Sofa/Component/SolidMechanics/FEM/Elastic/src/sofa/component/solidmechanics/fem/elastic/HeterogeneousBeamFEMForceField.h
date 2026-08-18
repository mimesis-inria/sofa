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
#include <sofa/component/solidmechanics/fem/elastic/config.h>

#include <sofa/component/solidmechanics/fem/elastic/BeamFEMForceField.h>
#include <sofa/type/vector.h>

namespace sofa::component::solidmechanics::fem::elastic
{

namespace _heterogeneousbeamfemforcefield_
{

using type::Vec3;

template<class DataTypes>
class HeterogeneousBeamFEMForceField
    : public _beamfemforcefield_::BeamFEMForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(HeterogeneousBeamFEMForceField, DataTypes),
               SOFA_TEMPLATE(_beamfemforcefield_::BeamFEMForceField, DataTypes));

    using Inherit = _beamfemforcefield_::BeamFEMForceField<DataTypes>;
    using typename Inherit::Real;
    using typename Inherit::Index;
    using typename Inherit::VecCoord;
    using typename Inherit::Deriv;
    using typename Inherit::Displacement;

    HeterogeneousBeamFEMForceField();
    ~HeterogeneousBeamFEMForceField() override = default;

    void init() override;
    void reinit() override;

    // BeamFEMForceField has reinitBeam(Index). We override it to apply per-section params.
    void reinitBeam(Index i) override;

    int getSectionForEdge(Index edgeIndex) const { return sectionForEdge(edgeIndex); };

    void LocalForceLarge(const VecCoord & x, int i, Index a, Index b, Deriv& fa, Deriv& fb);

protected:
    Data< sofa::type::vector<Index> > d_sectionStart;       ///< inclusive start edge index per section
    Data< sofa::type::vector<Index> > d_sectionEnd;         ///< exclusive end edge index per section

    Data< sofa::type::vector<SReal> > d_sectionYoung;       ///< E per section
    Data< sofa::type::vector<SReal> > d_sectionPoisson;     ///< nu per section
    Data< sofa::type::vector<SReal> > d_sectionRadius;      ///< radius per section
    Data< sofa::type::vector<SReal> > d_sectionInnerRadius; ///< inner radius per section

    sofa::core::objectmodel::MultiLink<HeterogeneousBeamFEMForceField<DataTypes>,sofa::core::objectmodel::BaseObject,sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK> l_sectionYoungLinks;

    sofa::type::vector<SReal> m_resolvedSectionYoung;
    sofa::type::vector<int> m_edgeToSection; ///< cache: edgeIndex -> sectionId, -1 if none
    Data< sofa::type::vector<int> > d_edgeToSection;
    void buildEdgeToSectionMap();
    int  sectionForEdge(Index edgeIndex) const;
    bool validateSections(Index nbEdges) const;
};

} // namespace _heterogeneousbeamfemforcefield_

// Export in the same way BeamFEMForceField does:
using _heterogeneousbeamfemforcefield_::HeterogeneousBeamFEMForceField;

} // namespace sofa::component::solidmechanics::fem::elastic