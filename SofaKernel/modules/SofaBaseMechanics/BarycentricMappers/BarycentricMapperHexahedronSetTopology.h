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
#ifndef SOFA_COMPONENT_MAPPING_BARYCENTRICMAPPERHEXAHEDRONSETTOPOLOGY_H
#define SOFA_COMPONENT_MAPPING_BARYCENTRICMAPPERHEXAHEDRONSETTOPOLOGY_H

#include <SofaBaseMechanics/BarycentricMappers/BarycentricMapperTopologyContainer.h>
#include <SofaBaseTopology/HexahedronSetTopologyContainer.h>
#include <SofaBaseTopology/HexahedronSetGeometryAlgorithms.h>

namespace sofa
{

namespace component
{

namespace mapping
{


using sofa::defaulttype::Mat3x3d;
using sofa::defaulttype::Vector3;
using sofa::defaulttype::Vec3dTypes;
using sofa::defaulttype::Vec3fTypes;
using sofa::defaulttype::ExtVec3Types;
typedef typename sofa::core::topology::BaseMeshTopology::Hexahedron Hexahedron;

/// Class allowing barycentric mapping computation on a HexahedronSetTopology
template<class In, class Out>
class BarycentricMapperHexahedronSetTopology : public BarycentricMapperTopologyContainer<In,Out,typename BarycentricMapper<In, Out>::MappingData3D,Hexahedron>
{
    typedef typename BarycentricMapper<In, Out>::MappingData3D MappingData;

public:
    SOFA_CLASS(SOFA_TEMPLATE2(BarycentricMapperHexahedronSetTopology,In,Out),
               SOFA_TEMPLATE4(BarycentricMapperTopologyContainer,In,Out,MappingData,Hexahedron));

    typedef typename Inherit1::Real Real;

    ~BarycentricMapperHexahedronSetTopology() override ;
    virtual helper::vector<Hexahedron> getElements() override;
    virtual helper::vector<SReal> getBaryCoef(const Real* f) override;
    helper::vector<SReal> getBaryCoef(const Real fx, const Real fy, const Real fz);
    void computeBase(Mat3x3d& base, const typename In::VecCoord& in, const Hexahedron& element) override;
    void computeCenter(Vector3& center, const typename In::VecCoord& in, const Hexahedron& element) override;
    void computeDistance(double& d, const Vector3& v) override;
    void addPointInElement(const int elementIndex, const SReal* baryCoords) override;

    int addPointInCube(const int index, const SReal* baryCoords) override;
    int setPointInCube(const int pointIndex, const int cubeIndex, const SReal* baryCoords) override;
    void applyOnePoint( const unsigned int& hexaId, typename Out::VecCoord& out, const typename In::VecCoord& in) override;
    void handleTopologyChange(core::topology::Topology* t) override;

protected:
    BarycentricMapperHexahedronSetTopology();
    BarycentricMapperHexahedronSetTopology(topology::HexahedronSetTopologyContainer* fromTopology,
                                           topology::PointSetTopologyContainer* toTopology);

    void setTopology(topology::HexahedronSetTopologyContainer* topology);

    topology::HexahedronSetTopologyContainer*		m_fromContainer {nullptr};
    topology::HexahedronSetGeometryAlgorithms<In>*	m_fromGeomAlgo  {nullptr};
    std::set<int> m_invalidIndex;

    using Inherit1::d_map;
    using Inherit1::m_matrixJ;
    using Inherit1::m_updateJ;
    using Inherit1::m_fromTopology;
};

#if !defined(SOFA_COMPONENT_MAPPING_BARYCENTRICMAPPERHEXAHEDRONSETTOPOLOGY_CPP)
extern template class SOFA_BASE_MECHANICS_API BarycentricMapperHexahedronSetTopology< Vec3dTypes, Vec3dTypes >;
extern template class SOFA_BASE_MECHANICS_API BarycentricMapperHexahedronSetTopology< Vec3dTypes, ExtVec3Types >;


#endif

}}}

#endif
