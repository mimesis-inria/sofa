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
#define SOFA_COMPONENT_ENGINE_SUBSETTOPOLOGY_CPP
#include <sofa/component/engine/select/SubsetTopology.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>

namespace sofa::component::engine::select
{

using namespace sofa::defaulttype;

void registerSubsetTopology(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(core::ObjectRegistrationData("Engine used to create subset topology given box, sphere, plan, ...")
        .add< SubsetTopology<Vec3Types> >()
        .add< SubsetTopology<Rigid3Types> >());
}

template class SOFA_COMPONENT_ENGINE_SELECT_API SubsetTopology<Vec3Types>;
template class SOFA_COMPONENT_ENGINE_SELECT_API SubsetTopology<Rigid3Types>;

} //namespace sofa::component::engine::select
