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
#define SOFA_COMPONENT_MAPPING_TUBULARMAPPING_CPP

#include <SofaMiscMapping/TubularMapping.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/behavior/MechanicalState.h>

namespace sofa
{

namespace component
{

namespace mapping
{

using namespace sofa::defaulttype;
using namespace core;
using namespace core::behavior;

// Register in the Factory
int TubularMappingClass = core::RegisterObject("Create a Tube around rigid points")
        .add< TubularMapping< Rigid3Types, Vec3dTypes > > ()
        .add< TubularMapping< Rigid3Types, ExtVec3Types > >()



        ;

template class SOFA_MISC_MAPPING_API TubularMapping< Rigid3Types, Vec3dTypes >;
template class SOFA_MISC_MAPPING_API TubularMapping< Rigid3Types, ExtVec3Types >;





} // namespace mapping

} // namespace component

} // namespace sofa

