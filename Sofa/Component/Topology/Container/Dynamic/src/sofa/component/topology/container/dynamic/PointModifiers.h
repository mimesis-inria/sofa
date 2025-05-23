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
#include <sofa/component/topology/container/dynamic/config.h>

#include <sofa/type/vector.h>
#include <sofa/core/topology/BaseTopology.h>

namespace sofa::component::topology::container::dynamic
{

using PointID = core::topology::BaseMeshTopology::PointID;


/**
* This class store all the info to create a new point in the mesh taking into account estimated id
* id of duplicated point if this point will be splitted due to a cut.
* This structure also store all the ancestors and coefficient to efficently add this point into the current mesh.
*/
class SOFA_COMPONENT_TOPOLOGY_CONTAINER_DYNAMIC_API PointToAdd
{
public:
    PointToAdd(PointID uniqueID, PointID idPoint,
        const sofa::type::vector<PointID>& ancestors,
        const sofa::type::vector<SReal>& coefs, SReal snapValue = 1_sreal)
        : m_uniqueID(uniqueID)
        , m_idPoint(idPoint)
        , m_ancestors(ancestors)
        , m_coefs(coefs)
        , m_snapValue(snapValue)
    {
        if (ancestors.size() != coefs.size())
        {
            msg_error("PointToAdd") << "Not the same size of ancestors: " << ancestors.size() << " vs coefs: " << coefs.size() << " given for PointToAdd id: " << uniqueID;
            m_ancestorType = sofa::geometry::ElementType::UNKNOWN;
            return;
        }
    }


    // bool return true if point is snap
    void updatePointIDForDuplication(bool split = true) 
    {
        if (m_ancestorType == sofa::geometry::ElementType::POINT)
        {
            m_idClone = m_idPoint;
            m_idPoint = m_ownerId;
            m_isSnapped = true;
        }
        else if (m_ancestorType == sofa::geometry::ElementType::EDGE)
        {
            m_idClone = m_idPoint + 1;
        }
        else if (m_ancestorType == sofa::geometry::ElementType::TRIANGLE)
        {
            // point in middle of a triangle. It should not be duplicated
        }

        if (split == false)
            m_idClone = sofa::InvalidID;
    }

    PointID getNbrNewPoint()
    {
        PointID nbr = 2;
        if (m_isSnapped)
            nbr--;
        if (m_idClone == sofa::InvalidID)
            nbr--;

        return nbr;
    }


    void printValue()
    {
        msg_info("PointToAdd") << "PTA: " << m_uniqueID << " | idPoint: " << m_idPoint << " | idClone: " << m_idClone << " | m_ancestorType: " << int(m_ancestorType);
        msg_info("PointToAdd") << "PTA: " << m_uniqueID << " | ancestors: " << m_ancestors << " | coefs: " << m_coefs;
    }


    /// Unique ID of this point structure. Will be a code combining ancestors ids
    PointID m_uniqueID;

    /// Future pointID of this pointToAdd
    PointID m_idPoint = sofa::InvalidID;
    /// Future pointID of this pointToAdd if this point is duplicated due to a cut
    PointID m_idClone = sofa::InvalidID;

    sofa::geometry::ElementType m_ancestorType = sofa::geometry::ElementType::UNKNOWN;

    core::topology::BaseMeshTopology::ElemID m_ownerId = sofa::InvalidID;

    /// List of ancestors (existing point ID of the mesh)
    sofa::type::vector<PointID> m_ancestors;
    /// List of corresponding coefficients 
    sofa::type::vector<SReal> m_coefs;

    bool m_isSnapped = false;

    SReal m_snapValue;
};



/// global function to compute and consistently generate a unique ID based on the IDs of its ancestors. The ID is determined by combining both ancestor IDs in ascending order, separated by a unit of 1e6.
inline PointID getUniqueId(PointID ancestor0, PointID ancestor1)
{
    PointID uniqID;
    if (ancestor0 > ancestor1)
        uniqID = 1000000 * ancestor0 + ancestor1;
    else
        uniqID = 1000000 * ancestor1 + ancestor0;

    return uniqID;
}

/// global function to compute and consistently generate a unique ID based on the IDs of its ancestors. The ID is determined by combining the 3 ancestor IDs in ascending order, separated by a unit of 1e6 and 1e12.
inline PointID getUniqueId(PointID ancestor0, PointID ancestor1, PointID ancestor2)
{
    PointID uniqID;
    if (ancestor0 > ancestor1) 
    {
        if (ancestor2 > ancestor0)
            uniqID = 1e12 * ancestor2 + 1e6 * ancestor0 + ancestor1;
        else if (ancestor2 < ancestor1)
            uniqID = 1e12 * ancestor0 + 1e6 * ancestor1 + ancestor2;
        else
            uniqID = 1e12 * ancestor0 + 1e6 * ancestor2 + ancestor1;
    }
    else // ancestor1 > ancestor0
    {
        if (ancestor2 > ancestor1)
            uniqID = 1e12 * ancestor2 + 1e6 * ancestor1 + ancestor0;
        else if (ancestor2 < ancestor0)
            uniqID = 1e12 * ancestor1 + 1e6 * ancestor0 + ancestor2;
        else
            uniqID = 1e12 * ancestor1 + 1e6 * ancestor2 + ancestor0;
    }

    return uniqID;
}


} //namespace sofa::component::topology::container::dynamic
