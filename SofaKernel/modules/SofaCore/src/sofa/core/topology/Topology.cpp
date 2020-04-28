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
#include "Topology.h"
#include <sofa/core/objectmodel/BaseNode.h>

namespace sofa
{

namespace core
{

namespace topology
{

bool Topology::insertInNode( objectmodel::BaseNode* node )
{
    node->addTopology(this);
    Inherit1::insertInNode(node);
    return true;
}

bool Topology::removeInNode( objectmodel::BaseNode* node )
{
    node->removeTopology(this);
    Inherit1::removeInNode(node);
    return true;
}

const sofa::helper::vector<Topology::index_type> Topology::InvalidSet;
const Topology::Edge Topology::InvalidEdge;
const Topology::Triangle Topology::InvalidTriangle;
const Topology::Quad Topology::InvalidQuad;
const Topology::Tetrahedron Topology::InvalidTetrahedron;
const Topology::Pyramid Topology::InvalidPyramid;
const Topology::Pentahedron Topology::InvalidPentahedron;
const Topology::Hexahedron Topology::InvalidHexahedron;


} // namespace topology

} // namespace core

} // namespace sofa

