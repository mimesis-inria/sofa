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
#include <SofaTopologyMapping/Hexa2TetraTopologicalMapping.h>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/core/ObjectFactory.h>

#include <SofaBaseTopology/HexahedronSetTopologyContainer.h>
#include <SofaBaseTopology/HexahedronSetTopologyModifier.h>

#include <SofaBaseTopology/TetrahedronSetTopologyContainer.h>
#include <SofaBaseTopology/TetrahedronSetTopologyModifier.h>

#include <sofa/core/topology/TopologyChange.h>

#include <SofaBaseTopology/GridTopology.h>

#include <sofa/defaulttype/Vec.h>
#include <map>
#include <sofa/defaulttype/VecTypes.h>

namespace sofa
{

namespace component
{

namespace topology
{

using namespace sofa::defaulttype;

using namespace sofa::component::topology;
using namespace sofa::core::topology;

// Register in the Factory
int Hexa2TetraTopologicalMappingClass = core::RegisterObject("Special case of mapping where HexahedronSetTopology is converted to TetrahedronSetTopology")
        .add< Hexa2TetraTopologicalMapping >()

        ;

// Implementation

Hexa2TetraTopologicalMapping::Hexa2TetraTopologicalMapping()
    : swapping(initData(&swapping, false, "swapping","Boolean enabling to swapp hexa-edges\n in order to avoid bias effect"))
{
}

Hexa2TetraTopologicalMapping::~Hexa2TetraTopologicalMapping()
{
}

void Hexa2TetraTopologicalMapping::init()
{    
    // INITIALISATION of TETRAHEDRAL mesh from HEXAHEDRAL mesh :

    // recheck models
    bool modelsOk = true;
    if (!fromModel)
    {
        msg_error() << "Pointer to input topology is invalid.";
        modelsOk = false;
    }

    if (!toModel)
    {
        msg_error() << "Pointer to output topology is invalid.";
        modelsOk = false;
    }
    else
    {
        TetrahedronSetTopologyModifier *to_tstm;
        toModel->getContext()->get(to_tstm);
        if (!to_tstm)
        {
            msg_error() << "No TetrahedronSetTopologyModifier found in the Tetrahedron topology Node.";
            modelsOk = false;
        }
    }

    if (!modelsOk)
    {
        this->m_componentstate = sofa::core::objectmodel::ComponentState::Invalid;
        return;
    }

    TetrahedronSetTopologyContainer *to_tstc;
    toModel->getContext()->get(to_tstc);
    // Clear output topology
    to_tstc->clear();

    // Set the same number of points
    toModel->setNbPoints(fromModel->getNbPoints());

    sofa::helper::vector <unsigned int>& Loc2GlobVec = *(Loc2GlobDataVec.beginEdit());

    Loc2GlobVec.clear();
    Glob2LocMap.clear();

    size_t nbcubes = fromModel->getNbHexahedra();

    // These values are only correct if the mesh is a grid topology
    int nx = 2;
    int ny = 1;
    //int nz = 1;
    {
        topology::GridTopology* grid = dynamic_cast<topology::GridTopology*>(fromModel.get());
        if (grid != NULL)
        {
            nx = grid->getNx()-1;
            ny = grid->getNy()-1;
            //nz = grid->getNz()-1;
        }
    }

    // Tesselation of each cube into 6 tetrahedra
    for (size_t i=0; i<nbcubes; i++)
    {
        core::topology::BaseMeshTopology::Hexa c = fromModel->getHexahedron(i);
#define swap(a,b) { int t = a; a = b; b = t; }
        // TODO : swap indexes where needed (currently crash in TriangleSetContainer)
        bool swapped = false;

        if(swapping.getValue())
        {
            if (!((i%nx)&1))
            {
                // swap all points on the X edges
                swap(c[0],c[1]);
                swap(c[3],c[2]);
                swap(c[4],c[5]);
                swap(c[7],c[6]);
                swapped = !swapped;
            }
            if (((i/nx)%ny)&1)
            {
                // swap all points on the Y edges
                swap(c[0],c[3]);
                swap(c[1],c[2]);
                swap(c[4],c[7]);
                swap(c[5],c[6]);
                swapped = !swapped;
            }
            if ((i/(nx*ny))&1)
            {
                // swap all points on the Z edges
                swap(c[0],c[4]);
                swap(c[1],c[5]);
                swap(c[2],c[6]);
                swap(c[3],c[7]);
                swapped = !swapped;
            }
        }
#undef swap
        if(!swapped)
        {
            to_tstc->addTetra(c[0],c[5],c[1],c[6]);
            to_tstc->addTetra(c[0],c[1],c[3],c[6]);
            to_tstc->addTetra(c[1],c[3],c[6],c[2]);
            to_tstc->addTetra(c[6],c[3],c[0],c[7]);
            to_tstc->addTetra(c[6],c[7],c[0],c[5]);
            to_tstc->addTetra(c[7],c[5],c[4],c[0]);
        }
        else
        {
            to_tstc->addTetra(c[0],c[5],c[6],c[1]);
            to_tstc->addTetra(c[0],c[1],c[6],c[3]);
            to_tstc->addTetra(c[1],c[3],c[2],c[6]);
            to_tstc->addTetra(c[6],c[3],c[7],c[0]);
            to_tstc->addTetra(c[6],c[7],c[5],c[0]);
            to_tstc->addTetra(c[7],c[5],c[0],c[4]);
        }
        for(int j=0; j<6; j++)
            Loc2GlobVec.push_back(i);
        Glob2LocMap[i] = (unsigned int)Loc2GlobVec.size()-1;
    }

    Loc2GlobDataVec.endEdit();

    // Need to fully init the target topology
    toModel->init();

    this->m_componentstate = sofa::core::objectmodel::ComponentState::Valid;
}

unsigned int Hexa2TetraTopologicalMapping::getFromIndex(unsigned int /*ind*/)
{

    return Topology::InvalidID;
}

void Hexa2TetraTopologicalMapping::updateTopologicalMappingTopDown()
{
    msg_warning() << "Method Hexa2TetraTopologicalMapping::updateTopologicalMappingTopDown() not yet implemented!";
// TODO...
}


} // namespace topology

} // namespace component

} // namespace sofa
