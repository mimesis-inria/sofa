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
#ifndef SOFA_COMPONENT_ENGINE_INDICESFROMVALUES_H
#define SOFA_COMPONENT_ENGINE_INDICESFROMVALUES_H
#include "config.h"

#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>


namespace sofa
{

namespace component
{

namespace engine
{

/**
 * This class returns the indices given a list of values.
 */
template <class T>
class IndicesFromValues : public core::DataEngine
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(IndicesFromValues,T),core::DataEngine);
    typedef T Value;
    typedef sofa::helper::vector<T> VecValue;
    typedef unsigned int Index;
    typedef sofa::helper::vector<Index> VecIndex;

protected:

    IndicesFromValues();

    ~IndicesFromValues() override;
public:
    void init() override;

    void reinit() override;

    void doUpdate() override;

    core::objectmodel::Data<VecValue> f_values; ///< input values
    core::objectmodel::Data<VecValue> f_global; ///< Global values, in which the input values are searched
    core::objectmodel::Data<VecIndex> f_indices; ///< Output indices of the given values, searched in global
    core::objectmodel::Data<VecIndex> f_otherIndices; ///< Output indices of the other values, (NOT the given ones) searched in global
    core::objectmodel::Data<bool> f_recursiveSearch; ///< if set to true, output are indices of the "global" data matching with one of the values

    virtual std::string getTemplateName() const override
    {
        return templateName(this);
    }
    static std::string templateName(const IndicesFromValues<T>* = NULL)
    {
        return sofa::defaulttype::DataTypeName<T>::name();
    }
};

#if  !defined(SOFA_COMPONENT_ENGINE_INDICESFROMVALUES_CPP)
extern template class SOFA_GENERAL_ENGINE_API IndicesFromValues<std::string>;
extern template class SOFA_GENERAL_ENGINE_API IndicesFromValues<int>;
extern template class SOFA_GENERAL_ENGINE_API IndicesFromValues<unsigned int>;
extern template class SOFA_GENERAL_ENGINE_API IndicesFromValues< helper::fixed_array<unsigned int, 2> >;
extern template class SOFA_GENERAL_ENGINE_API IndicesFromValues< helper::fixed_array<unsigned int, 3> >;
extern template class SOFA_GENERAL_ENGINE_API IndicesFromValues< helper::fixed_array<unsigned int, 4> >;
extern template class SOFA_GENERAL_ENGINE_API IndicesFromValues< helper::fixed_array<unsigned int, 8> >;
extern template class SOFA_GENERAL_ENGINE_API IndicesFromValues<double>;
extern template class SOFA_GENERAL_ENGINE_API IndicesFromValues<defaulttype::Vec2d>;
extern template class SOFA_GENERAL_ENGINE_API IndicesFromValues<defaulttype::Vec3d>;
// extern template class SOFA_GENERAL_ENGINE_API IndicesFromValues<defaulttype::Rigid2Types::Coord>;
// extern template class SOFA_GENERAL_ENGINE_API IndicesFromValues<defaulttype::Rigid2Types::Deriv>;
// extern template class SOFA_GENERAL_ENGINE_API IndicesFromValues<defaulttype::Rigid3Types::Coord>;
// extern template class SOFA_GENERAL_ENGINE_API IndicesFromValues<defaulttype::Rigid3Types::Deriv>;
 
#endif

} // namespace engine

} // namespace component

} // namespace sofa

#endif
