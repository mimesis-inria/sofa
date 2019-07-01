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
#define SOFA_COMPONENT_ENGINE_MERGEMESHES_CPP
#include <SofaGeneralEngine/MergeMeshes.inl>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace engine
{

int MergeMeshesClass = core::RegisterObject("Merge several meshes")
        .add< MergeMeshes<defaulttype::Vec3Types> >(true) // default template
        .add< MergeMeshes<defaulttype::Vec1Types> >()
        .add< MergeMeshes<defaulttype::Vec2Types> >()
        .add< MergeMeshes<defaulttype::Rigid2Types> >()
        .add< MergeMeshes<defaulttype::Rigid3Types> >()
 
        ;

template class SOFA_GENERAL_ENGINE_API MergeMeshes<defaulttype::Vec1Types>;
template class SOFA_GENERAL_ENGINE_API MergeMeshes<defaulttype::Vec2Types>;
template class SOFA_GENERAL_ENGINE_API MergeMeshes<defaulttype::Vec3Types>;
template class SOFA_GENERAL_ENGINE_API MergeMeshes<defaulttype::Rigid2Types>;
template class SOFA_GENERAL_ENGINE_API MergeMeshes<defaulttype::Rigid3Types>;
 


} // namespace engine

} // namespace component

} // namespace sofa

