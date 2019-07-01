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
#ifndef SOFA_COMPONENT_ENGINE_COMPLEMENTARYROI_H
#define SOFA_COMPONENT_ENGINE_COMPLEMENTARYROI_H
#include "config.h"

#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/helper/vectorData.h>

namespace sofa
{

namespace component
{

namespace engine
{

/**
 * Output the positions and their indices in the global mesh not in the specified sets
 *
 * @todo make it general as other ROI (edges, triangles,...)
 *
 * @author Thomas Lemaire @date 2014
 */
template <class DataTypes>
class ComplementaryROI : public core::DataEngine
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(ComplementaryROI, DataTypes), core::DataEngine);

    typedef typename DataTypes::VecCoord VecCoord;
    typedef core::topology::BaseMeshTopology::index_type index_type;
    typedef core::topology::BaseMeshTopology::SetIndex SetIndex;


    ComplementaryROI();
    ~ComplementaryROI() override;

    /// Update
    void doUpdate() override;

    /// Parse the given description to assign values to this object's fields and potentially other parameters
    void parse ( sofa::core::objectmodel::BaseObjectDescription* arg ) override;

    /// Assign the field values stored in the given map of name -> value pairs
    void parseFields ( const std::map<std::string,std::string*>& str ) override;


    void init() override;
    void reinit() override;

    virtual std::string getTemplateName() const override;

    static std::string templateName(const ComplementaryROI<DataTypes>* = NULL);

protected:

    /// inputs
    /// @{
    Data<VecCoord> d_position; ///< input positions
    Data<unsigned int> d_nbSet; ///< number of sets
    helper::vectorData< SetIndex > vd_setIndices; ///< for each set, indices of the included points
    /// @}

    /// outputs
    /// @{
    Data<SetIndex> d_indices; ///< ROI indices
    Data<VecCoord> d_pointsInROI; ///< ROI positions
    /// @}

};

#if  !defined(SOFA_COMPONENT_ENGINE_COMPLEMENTARYROI_CPP)
extern template class SOFA_GENERAL_ENGINE_API ComplementaryROI<defaulttype::Vec3Types>;
 

#endif

} // namespace engine

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_ENGINE_COMPLEMENTARYROI_H
