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
#define SOFA_COMPONENT_MAPPING_DISTANCETOPLANEMAPPING_CPP
#include <sofa/component/mapping/linear/DistanceToPlaneMapping.inl>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/component/mapping/linear/config.h>

namespace sofa::component::mapping::linear
{

using namespace sofa::defaulttype;

void registerDistanceToPlaneMapping(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(core::ObjectRegistrationData("Mapping that computes the distance to a plane")
        .add< DistanceToPlaneMapping< Vec3Types > >()
        .add< DistanceToPlaneMapping< Vec2Types > >()
        .add< DistanceToPlaneMapping< Vec6Types > >()
        .add< DistanceToPlaneMapping< Rigid3Types > >()
        .add< DistanceToPlaneMapping< Rigid2Types > >());
}

template class SOFA_COMPONENT_MAPPING_LINEAR_API DistanceToPlaneMapping< Vec3Types > ;
template class SOFA_COMPONENT_MAPPING_LINEAR_API DistanceToPlaneMapping< Vec2Types > ;
template class SOFA_COMPONENT_MAPPING_LINEAR_API DistanceToPlaneMapping< Vec6Types > ;
template class SOFA_COMPONENT_MAPPING_LINEAR_API DistanceToPlaneMapping< Rigid3Types >;
template class SOFA_COMPONENT_MAPPING_LINEAR_API DistanceToPlaneMapping< Rigid2Types >;

} // namespace sofa::component::mapping::linear
