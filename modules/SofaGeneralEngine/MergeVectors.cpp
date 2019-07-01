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
#define SOFA_COMPONENT_ENGINE_MERGEVECTORS_CPP
#include <SofaGeneralEngine/MergeVectors.inl>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace engine
{

int MergeVectorsClass = core::RegisterObject("Apply a merge operation to combine several inputs")
    .add< MergeVectors< helper::vector<double> > >(true)
    .add< MergeVectors< helper::vector<int> > >()
    .add< MergeVectors< helper::vector<bool> > >()
    //.add< MergeVectors< helper::vector<std::string> > >()
    .add< MergeVectors< helper::vector<defaulttype::Vec2u> > >()
    .add< MergeVectors< helper::vector<defaulttype::Vec2d> > >()
    .add< MergeVectors< helper::vector<defaulttype::Vec3d> > >()
    .add< MergeVectors< helper::vector<defaulttype::Vec4d> > >()
    .add< MergeVectors< defaulttype::Rigid2Types::VecCoord > >()
    .add< MergeVectors< defaulttype::Rigid2Types::VecDeriv > >()
    .add< MergeVectors< defaulttype::Rigid3Types::VecCoord > >()
    .add< MergeVectors< defaulttype::Rigid3Types::VecDeriv > >()
 
        ;

template class SOFA_GENERAL_ENGINE_API MergeVectors< helper::vector<int> >;
template class SOFA_GENERAL_ENGINE_API MergeVectors< helper::vector<bool> >;
//template class SOFA_GENERAL_ENGINE_API MergeVectors< helper::vector<std::string> >;
template class SOFA_GENERAL_ENGINE_API MergeVectors< helper::vector<defaulttype::Vec2u> >;

template class SOFA_GENERAL_ENGINE_API MergeVectors< helper::vector<double> >;
template class SOFA_GENERAL_ENGINE_API MergeVectors< helper::vector<defaulttype::Vec2d> >;
template class SOFA_GENERAL_ENGINE_API MergeVectors< helper::vector<defaulttype::Vec3d> >;
template class SOFA_GENERAL_ENGINE_API MergeVectors< helper::vector<defaulttype::Vec4d> >;
template class SOFA_GENERAL_ENGINE_API MergeVectors< defaulttype::Rigid2Types::VecCoord >;
template class SOFA_GENERAL_ENGINE_API MergeVectors< defaulttype::Rigid2Types::VecDeriv >;
template class SOFA_GENERAL_ENGINE_API MergeVectors< defaulttype::Rigid3Types::VecCoord >;
template class SOFA_GENERAL_ENGINE_API MergeVectors< defaulttype::Rigid3Types::VecDeriv >;
 


} // namespace constraint

} // namespace component

} // namespace sofa

