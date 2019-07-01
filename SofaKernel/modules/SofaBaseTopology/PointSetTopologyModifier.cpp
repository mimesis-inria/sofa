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
#include <SofaBaseTopology/PointSetTopologyModifier.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/StateChangeVisitor.h>
#include <sofa/simulation/Simulation.h>
#include <sofa/simulation/TopologyChangeVisitor.h>
#include <sofa/core/topology/TopologyChange.h>
#include <SofaBaseTopology/PointSetTopologyContainer.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/AdvancedTimer.h>

namespace sofa
{

namespace component
{

namespace topology
{
int PointSetTopologyModifierClass = core::RegisterObject("Point set topology modifier")
        .add< PointSetTopologyModifier >();

using namespace std;
using namespace sofa::defaulttype;
using namespace sofa::core::behavior;
using namespace sofa::core::topology;

void PointSetTopologyModifier::init()
{
    core::topology::TopologyModifier::init();
    this->getContext()->get(m_container);
}


void PointSetTopologyModifier::swapPoints(const int i1, const int i2)
{
    PointsIndicesSwap *e2 = new PointsIndicesSwap( i1, i2 );
    addStateChange(e2);
    propagateStateChanges();

    PointsIndicesSwap *e = new PointsIndicesSwap( i1, i2 );
    this->addTopologyChange(e);
}


void PointSetTopologyModifier::addPointsProcess(const size_t nPoints)
{
    m_container->addPoints(nPoints);
}

void PointSetTopologyModifier::addPointsWarning(const size_t nPoints, const bool addDOF)
{
    m_container->setPointTopologyToDirty();
    const size_t startIndex = m_container->getNbPoints()-nPoints;
    sofa::helper::vector <PointID> indices; indices.resize(nPoints);
    for(size_t i=0; i<nPoints; ++i)
        indices[i] = PointID(startIndex+i);

    if(addDOF)
    {
        PointsAdded *e2 = new PointsAdded(nPoints,indices);
        addStateChange(e2);
        propagateStateChanges();
    }

    // Warning that vertices just got created
    PointsAdded *e = new PointsAdded(nPoints,indices);
    this->addTopologyChange(e);
}


void PointSetTopologyModifier::addPointsWarning(const size_t nPoints,
        const sofa::helper::vector< sofa::helper::vector< PointID > > &ancestors,
        const sofa::helper::vector< sofa::helper::vector< double       > >& coefs,
        const bool addDOF)
{
    m_container->setPointTopologyToDirty();
    const size_t startIndex = m_container->getNbPoints()-nPoints;
    sofa::helper::vector <PointID> indices; indices.resize(nPoints);
    for(size_t i=0; i<nPoints; ++i)
        indices[i] = PointID(startIndex+i);

    if(addDOF)
    {
        PointsAdded *e2 = new PointsAdded(nPoints, indices, ancestors, coefs);
        addStateChange(e2);
        propagateStateChanges();
    }

    // Warning that vertices just got created
    PointsAdded *e = new PointsAdded(nPoints, indices, ancestors, coefs);
    this->addTopologyChange(e);
}


void PointSetTopologyModifier::addPointsWarning(const size_t nPoints,
        const sofa::helper::vector< core::topology::PointAncestorElem >& ancestorElems,
        const bool addDOF)
{
    using namespace sofa::core::topology;

    assert(ancestorElems.size() == nPoints);

    m_container->setPointTopologyToDirty();

    // Compute standard points construction info based on ancestor points
    // related to topology elems and local coordinates

    const size_t startIndex = m_container->getNbPoints() - nPoints;

    helper::vector< PointID > newPointIndices;
    helper::vector< helper::vector< PointID > > ancestorPointIndices;
    helper::vector< helper::vector< double       > > baryCoefs;

    newPointIndices.resize(nPoints);
    ancestorPointIndices.resize(nPoints);
    baryCoefs.resize(nPoints);

    for(size_t i = 0; i < nPoints; ++i)
    {
        newPointIndices[i] = PointID(startIndex + i);
        PointID ancestorIndex = ancestorElems[i].index;
        // check if this new point has indeed an ancestor.
        if (ancestorIndex != core::topology::BaseMeshTopology::InvalidID )
        {
            core::topology::PointAncestorElem::LocalCoords localCoords = ancestorElems[i].localCoords;
            switch (ancestorElems[i].type)
            {
            case POINT :
                {
                    ancestorPointIndices[i].resize(1);
                    ancestorPointIndices[i][0] = ancestorIndex;

                    baryCoefs[i].resize(1);
                    baryCoefs[i][0] = 1;
                    break;
                }
            case EDGE :
                {
                    const core::topology::Topology::Edge& e = m_container->getEdge(ancestorIndex);
                    ancestorPointIndices[i].resize(2);
                    ancestorPointIndices[i][0] = e[0];
                    ancestorPointIndices[i][1] = e[1];

                    baryCoefs[i].resize(2);
                    baryCoefs[i][0] = 1 - localCoords[0];
                    baryCoefs[i][1] = localCoords[0];
                    break;
                }
            case TRIANGLE :
                {
                    const core::topology::Topology::Triangle& t = m_container->getTriangle(ancestorIndex);
                    ancestorPointIndices[i].resize(3);
                    ancestorPointIndices[i][0] = t[0];
                    ancestorPointIndices[i][1] = t[1];
                    ancestorPointIndices[i][2] = t[2];

                    baryCoefs[i].resize(3);
                    baryCoefs[i][0] = 1 - localCoords[0] - localCoords[1];
                    baryCoefs[i][1] = localCoords[0];
                    baryCoefs[i][2] = localCoords[1];
                    break;
                }
            case TETRAHEDRON :
                {
                    const core::topology::Topology::Tetrahedron& t = m_container->getTetrahedron(ancestorIndex);
                    ancestorPointIndices[i].resize(4);
                    ancestorPointIndices[i][0] = t[0];
                    ancestorPointIndices[i][1] = t[1];
                    ancestorPointIndices[i][2] = t[2];
                    ancestorPointIndices[i][3] = t[3];

                    baryCoefs[i].resize(4);
                    baryCoefs[i][0] = 1 - localCoords[0] - localCoords[1] - localCoords[2];
                    baryCoefs[i][1] = localCoords[0];
                    baryCoefs[i][2] = localCoords[1];
                    baryCoefs[i][3] = localCoords[2];

                    break;
                }
            case QUAD :
                {
                    const core::topology::Topology::Quad& q = m_container->getQuad(ancestorIndex);
                    ancestorPointIndices[i].resize(4);
                    ancestorPointIndices[i][0] = q[0];
                    ancestorPointIndices[i][1] = q[1];
                    ancestorPointIndices[i][2] = q[2];
                    ancestorPointIndices[i][3] = q[3];

                    baryCoefs[i].resize(4);
                    baryCoefs[i][0] = (1 - localCoords[0])*(1 - localCoords[1]);
                    baryCoefs[i][1] = (    localCoords[0])*(1 - localCoords[1]);
                    baryCoefs[i][2] = (    localCoords[0])*(    localCoords[1]);
                    baryCoefs[i][3] = (1 - localCoords[0])*(    localCoords[1]);

                    break;
                }
            case HEXAHEDRON :
                {
                    const core::topology::Topology::Hexahedron& h = m_container->getHexahedron(ancestorIndex);
                    ancestorPointIndices[i].resize(8);
                    ancestorPointIndices[i][0] = h[0];
                    ancestorPointIndices[i][1] = h[1];
                    ancestorPointIndices[i][2] = h[2];
                    ancestorPointIndices[i][3] = h[3];
                    ancestorPointIndices[i][4] = h[4];
                    ancestorPointIndices[i][5] = h[5];
                    ancestorPointIndices[i][6] = h[6];
                    ancestorPointIndices[i][7] = h[7];

                    baryCoefs[i].resize(8);
                    baryCoefs[i][0] = (1 - localCoords[0])*(1 - localCoords[1])*(1 - localCoords[2]);
                    baryCoefs[i][1] = (    localCoords[0])*(1 - localCoords[1])*(1 - localCoords[2]);
                    baryCoefs[i][2] = (    localCoords[0])*(    localCoords[1])*(1 - localCoords[2]);
                    baryCoefs[i][3] = (1 - localCoords[0])*(    localCoords[1])*(1 - localCoords[2]);
                    baryCoefs[i][4] = (1 - localCoords[0])*(1 - localCoords[1])*(    localCoords[2]);
                    baryCoefs[i][5] = (    localCoords[0])*(1 - localCoords[1])*(    localCoords[2]);
                    baryCoefs[i][6] = (    localCoords[0])*(    localCoords[1])*(    localCoords[2]);
                    baryCoefs[i][7] = (1 - localCoords[0])*(    localCoords[1])*(    localCoords[2]);

                    break;
                }
            default :
                msg_error() << "Unsupported ancestor primitive type in addPointsWarning";
                break;
            }
        }
    }

    if(addDOF)
    {
        PointsAdded *e2 = new PointsAdded(nPoints, newPointIndices, ancestorElems,
            ancestorPointIndices, baryCoefs);
        addStateChange(e2);
        propagateStateChanges();
    }

    // Warning that vertices just got created
    PointsAdded *e = new PointsAdded(nPoints, newPointIndices, ancestorElems,
        ancestorPointIndices, baryCoefs);
    this->addTopologyChange(e);
}


void PointSetTopologyModifier::addPoints(const size_t nPoints,
                                         const bool addDOF)
{
    sofa::helper::AdvancedTimer::stepBegin("addPoints");

    sofa::helper::AdvancedTimer::stepBegin("addPointsProcess");
    addPointsProcess(nPoints);

    sofa::helper::AdvancedTimer::stepNext ("addPointsProcess", "addPointsWarning");
    addPointsWarning(nPoints, addDOF);

    sofa::helper::AdvancedTimer::stepNext ("addPointsWarning", "propagateTopologicalChanges");
    propagateTopologicalChanges();
    sofa::helper::AdvancedTimer::stepEnd("propagateTopologicalChanges");

    sofa::helper::AdvancedTimer::stepEnd("addPoints");
}

void PointSetTopologyModifier::addPoints(const size_t nPoints,
     const sofa::helper::vector< sofa::helper::vector< PointID > >& ancestors,
     const sofa::helper::vector< sofa::helper::vector< double> >& coefs,
     const bool addDOF)
{
    sofa::helper::AdvancedTimer::stepBegin("addPoints with ancestors");

    sofa::helper::AdvancedTimer::stepBegin("addPointsProcess");
    addPointsProcess(nPoints);

    sofa::helper::AdvancedTimer::stepNext ("addPointsProcess", "addPointsWarning");
    addPointsWarning(nPoints, ancestors, coefs, addDOF);

    sofa::helper::AdvancedTimer::stepNext ("addPointsWarning", "propagateTopologicalChanges");
    propagateTopologicalChanges();
    sofa::helper::AdvancedTimer::stepEnd("propagateTopologicalChanges");

    sofa::helper::AdvancedTimer::stepEnd("addPoints with ancestors");
}

void PointSetTopologyModifier::addPoints(const size_t nPoints,
     const sofa::helper::vector< core::topology::PointAncestorElem >& srcElems,
     const bool addDOF)
{
    addPointsProcess(nPoints);
    addPointsWarning(nPoints, srcElems, addDOF);
    propagateTopologicalChanges();
}

void PointSetTopologyModifier::movePointsProcess (const sofa::helper::vector <PointID>& id,
        const sofa::helper::vector< sofa::helper::vector< PointID > >& ancestors,
        const sofa::helper::vector< sofa::helper::vector< double > >& coefs,
        const bool moveDOF)
{
    m_container->setPointTopologyToDirty();

    if(moveDOF)
    {
        PointsMoved *ev = new PointsMoved(id, ancestors, coefs);
        addStateChange(ev);
        propagateStateChanges();
    }

    m_container->setPointTopologyToDirty();

    // Warning that vertices just been moved
    PointsMoved *ev2 = new PointsMoved(id, ancestors, coefs);
    this->addTopologyChange(ev2);

}



void PointSetTopologyModifier::removePointsWarning(sofa::helper::vector<PointID> &indices,
        const bool removeDOF)
{
    sofa::helper::AdvancedTimer::stepBegin("removePointsWarning");
    m_container->setPointTopologyToDirty();

    // sort points so that they are removed in a descending order
    std::sort( indices.begin(), indices.end(), std::greater<PointID>() );

    // Warning that these vertices will be deleted
    PointsRemoved *e = new PointsRemoved(indices);
    this->addTopologyChange(e);

    if(removeDOF)
    {
        PointsRemoved *e2 = new PointsRemoved(indices);
        addStateChange(e2);
    }
    sofa::helper::AdvancedTimer::stepEnd("removePointsWarning");
}


void PointSetTopologyModifier::removePointsProcess(const sofa::helper::vector<PointID> & indices,
        const bool removeDOF)
{
    sofa::helper::AdvancedTimer::stepBegin("removePointsProcess");
    if(removeDOF)
    {
        propagateStateChanges();
    }
    m_container->removePoints(indices.size());

    sofa::helper::AdvancedTimer::stepEnd("removePointsProcess");
}


void PointSetTopologyModifier::renumberPointsWarning( const sofa::helper::vector<PointID> &index,
        const sofa::helper::vector<PointID> &inv_index,
        const bool renumberDOF)
{
    // Warning that these vertices will be deleted
    PointsRenumbering *e = new PointsRenumbering(index, inv_index);
    this->addTopologyChange(e);

    if(renumberDOF)
    {
        PointsRenumbering *e2 = new PointsRenumbering(index, inv_index);
        addStateChange(e2);
    }
}


void PointSetTopologyModifier::renumberPointsProcess( const sofa::helper::vector<PointID> &/*index*/,
        const sofa::helper::vector<PointID> &/*inv_index*/,
        const bool renumberDOF)
{
    if(renumberDOF)
    {
        propagateStateChanges();
    }
}

void PointSetTopologyModifier::propagateTopologicalChanges()
{
    if (m_container->beginChange() == m_container->endChange()) return; // nothing to do if no event is stored

    this->propagateTopologicalEngineChanges();
    
    sofa::core::ExecParams* params = sofa::core::ExecParams::defaultInstance();
    sofa::simulation::TopologyChangeVisitor a(params, m_container);

    getContext()->executeVisitor(&a);

    // remove the changes we just propagated, so that we don't send them again next time
    m_container->resetTopologyChangeList();
}

void PointSetTopologyModifier::propagateTopologicalChangesWithoutReset()
{
    if (m_container->beginChange() == m_container->endChange()) return; // nothing to do if no event is stored
    sofa::core::ExecParams* params = sofa::core::ExecParams::defaultInstance();
    sofa::simulation::TopologyChangeVisitor a(params, m_container);

    getContext()->executeVisitor(&a);

    //TODO: temporary code to test topology engine pipeline. Commented by default for the moment
    //this->propagateTopologicalEngineChanges();

}


void PointSetTopologyModifier::propagateTopologicalEngineChanges()
{
    if (m_container->beginChange() == m_container->endChange()) // nothing to do if no event is stored
        return;

    if (!m_container->isPointTopologyDirty()) // triangle Data has not been touched
        return;

    sofa::helper::AdvancedTimer::stepBegin("PointSetTopologyModifier::propagateTopologicalEngineChanges");
    // get directly the list of engines created at init: case of removing.... for the moment
    std::list<sofa::core::topology::TopologyEngine *>::iterator it;

    for ( it = m_container->m_enginesList.begin(); it!=m_container->m_enginesList.end(); ++it)
    {
        // no need to dynamic cast this time? TO BE CHECKED!
        sofa::core::topology::TopologyEngine* topoEngine = (*it);
        if (topoEngine->isDirty())
        {
            topoEngine->update();
        }
    }

    m_container->cleanPointTopologyFromDirty();
    sofa::helper::AdvancedTimer::stepBegin("PointSetTopologyModifier::propagateTopologicalEngineChanges");
}

void PointSetTopologyModifier::propagateStateChanges()
{
    if (m_container->beginStateChange() == m_container->endStateChange()) return; // nothing to do if no event is stored
    sofa::core::ExecParams* params = sofa::core::ExecParams::defaultInstance();
    sofa::simulation::StateChangeVisitor a(params, m_container);
    getContext()->executeVisitor(&a);

    // remove the changes we just propagated, so that we don't send then again next time
    m_container->resetStateChangeList();
}

void PointSetTopologyModifier::notifyEndingEvent()
{
    sofa::core::topology::EndingEvent *e=new sofa::core::topology::EndingEvent();
    m_container->addTopologyChange(e);
}

} // namespace topology

} // namespace component

} // namespace sofa

