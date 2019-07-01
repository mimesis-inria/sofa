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

#include <SofaBaseTopology/QuadSetTopologyContainer.h>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/core/ObjectFactory.h>


namespace sofa
{
namespace component
{
namespace topology
{

using namespace std;
using namespace sofa::defaulttype;

int QuadSetTopologyContainerClass = core::RegisterObject("Quad set topology container")
        .add< QuadSetTopologyContainer >()
        ;

QuadSetTopologyContainer::QuadSetTopologyContainer()
    : EdgeSetTopologyContainer()
    , d_quad(initData(&d_quad, "quads", "List of quad indices"))
{
}


void QuadSetTopologyContainer::addQuad( int a, int b, int c, int d )
{
    helper::WriteAccessor< Data< sofa::helper::vector<Quad> > > m_quad = d_quad;
    m_quad.push_back(Quad(a,b,c,d));
    if (a >= getNbPoints()) setNbPoints( a+1 );
    if (b >= getNbPoints()) setNbPoints( b+1 );
    if (c >= getNbPoints()) setNbPoints( c+1 );
    if (d >= getNbPoints()) setNbPoints( d+1 );
}

void QuadSetTopologyContainer::init()
{
    d_quad.updateIfDirty(); // make sure m_quad is up to date

    helper::ReadAccessor< Data< sofa::helper::vector<Quad> > > m_quads = d_quad;
    if (!m_quads.empty())
    {
        for (size_t i=0; i<m_quads.size(); ++i)
        {
            for(PointID j=0; j<4; ++j)
            {
                int a = m_quads[i][j];
                if (a >= getNbPoints()) setNbPoints(a+1);
            }
        }
    }

    // only init if triangles are present at init.
    if (!m_quads.empty())
        initTopology();
}

void QuadSetTopologyContainer::initTopology()
{
    // Force creation of Edge Neighboordhood buffers.
    EdgeSetTopologyContainer::initTopology();

    // Create triangle cross element buffers.
    createEdgesInQuadArray();
    createQuadsAroundVertexArray();
    createQuadsAroundEdgeArray();
}

void QuadSetTopologyContainer::createQuadSetArray()
{
	if (CHECK_TOPOLOGY)
        msg_error() << "createQuadSetArray method must be implemented by a child topology.";

}

void QuadSetTopologyContainer::createQuadsAroundVertexArray()
{
    // first clear potential previous buffer
    clearQuadsAroundVertex();

    if(!hasQuads()) // this method should only be called when quads exist
        createQuadSetArray();

    if(hasQuadsAroundVertex()) // created by upper topology
        return;
    
    helper::ReadAccessor< Data< sofa::helper::vector<Quad> > > m_quad = d_quad;
    if (m_quad.empty())
    {
        msg_warning() << "QuadsAroundVertex buffer can't be created as no quads are present in this topology.";
        return;
    }

    if (getNbPoints() == 0) // in case only Data have been copied and not going thourgh AddTriangle methods.
        this->setNbPoints(d_initPoints.getValue().size());

    m_quadsAroundVertex.resize(getNbPoints());
    for (size_t i=0; i<m_quad.size(); ++i)
    {
        // adding quad i in the quad shell of all points
        for (size_t j=0; j<4; ++j)
        {
            m_quadsAroundVertex[ m_quad[i][j] ].push_back((QuadID)i);
        }
    }
}

void QuadSetTopologyContainer::createQuadsAroundEdgeArray()
{
    // first clear potential previous buffer
    clearQuadsAroundEdge();

    if(!hasQuads()) // this method should only be called when quads exist
        createQuadSetArray();

    if (hasQuadsAroundEdge()) // created by upper topology
        return;

    const size_t numQuads = getNumberOfQuads();
    if(numQuads == 0) // this method should only be called when quads exist
    {
        msg_warning() << "QuadsAroundEdge buffer can't be created as no quads are present in this topology.";
        return;
    }

    if(!hasEdges()) // this method should only be called when edges exist
        createEdgeSetArray();

    const size_t numEdges = getNumberOfEdges();
    if (numEdges == 0)
    {
        msg_warning() << "QuadsAroundEdge buffer can't be created as no edges are present in this topology.";
        return;
    }

    if(!hasEdgesInQuad())
        createEdgesInQuadArray();

    if (m_edgesInQuad.empty())
    {
        msg_warning() << "TrianglesAroundEdge buffer can't be created as EdgesInQuad buffer creation failed.";
        return;
    }

    m_quadsAroundEdge.resize(numEdges);
    for (size_t i=0; i<numQuads; ++i)
    {
        // adding quad i in the quad shell of all edges
        for (size_t j=0; j<4; ++j)
        {
            m_quadsAroundEdge[ m_edgesInQuad[i][j] ].push_back((QuadID)i);
        }
    }
}

void QuadSetTopologyContainer::createEdgeSetArray()
{
    if(!hasQuads()) // this method should only be called when quads exist
        createQuadSetArray();

    if(hasEdges())
    {
        // clear edges and all shells that depend on edges
        EdgeSetTopologyContainer::clear();

        if(hasEdgesInQuad())
            clearEdgesInQuad();

        if(hasQuadsAroundEdge())
            clearQuadsAroundEdge();
    }

    // create a temporary map to find redundant edges
    std::map<Edge, EdgeID> edgeMap;
    helper::WriteAccessor< Data< sofa::helper::vector<Edge> > > m_edge = d_edge;
    helper::ReadAccessor< Data< sofa::helper::vector<Quad> > > m_quad = d_quad;

    for (size_t i=0; i<m_quad.size(); ++i)
    {
        const Quad &t = m_quad[i];
        for(size_t j=0; j<4; ++j)
        {
            const PointID v1 = t[(j+1)%4];
            const PointID v2 = t[(j+2)%4];

            // sort vertices in lexicographic order
            const Edge e = ((v1<v2) ? Edge(v1,v2) : Edge(v2,v1));

            if(edgeMap.find(e) == edgeMap.end())
            {
                // edge not in edgeMap so create a new one
                const EdgeID edgeIndex = (EdgeID)edgeMap.size();
                edgeMap[e] = edgeIndex;
                m_edge.push_back(e);
            }
        }
    }
}

void QuadSetTopologyContainer::createEdgesInQuadArray()
{
    // first clear potential previous buffer
    clearEdgesInQuad();

    if(!hasQuads()) // this method should only be called when quads exist
        createQuadSetArray();

    if (hasEdgesInQuad()) // created by upper topology
        return;

    if(!hasEdges()) // this method should only be called when edges exist
        createEdgeSetArray();

    if (d_edge.getValue().empty())
    {
        msg_warning() << "EdgesInQuad buffer can't be created as no edges are present in this topology.";
        return;
    }

    const size_t numQuads = getNumberOfQuads();
    m_edgesInQuad.resize( numQuads );
    helper::ReadAccessor< Data< sofa::helper::vector<Quad> > > m_quad = d_quad;

    for(size_t i=0; i<numQuads; ++i)
    {
        const Quad &t = m_quad[i];
        // adding edge i in the edge shell of both points
        for (size_t j=0; j<4; ++j)
        {
            EdgeID edgeIndex = getEdgeIndex(t[(j+1)%4],t[(j+2)%4]);
            assert(edgeIndex != InvalidID);
            m_edgesInQuad[i][j]=edgeIndex;
        }
    }
}

const sofa::helper::vector<QuadSetTopologyContainer::Quad> &QuadSetTopologyContainer::getQuadArray()
{
    if(CHECK_TOPOLOGY)
        msg_warning_when(!hasQuads() && getNbPoints()>0) << "Quad array is empty with " << getNbPoints() << " vertices.";

    return d_quad.getValue();
}

const QuadSetTopologyContainer::Quad QuadSetTopologyContainer::getQuad (QuadID id)
{
    if ((size_t)id >= getNbQuads())
        return Quad(InvalidID, InvalidID, InvalidID, InvalidID);
    else
        return (d_quad.getValue())[id];
}


QuadSetTopologyContainer::QuadID QuadSetTopologyContainer::getQuadIndex(PointID v1, PointID v2, PointID v3, PointID v4)
{
    if(!hasQuadsAroundVertex())
    {
        if(CHECK_TOPOLOGY)
            msg_warning() << "QuadsAroundVertex array is empty with " << getNbPoints() << " vertices.";

        return InvalidID;
    }

    sofa::helper::vector<QuadID> set1 = getQuadsAroundVertex(v1);
    sofa::helper::vector<QuadID> set2 = getQuadsAroundVertex(v2);
    sofa::helper::vector<QuadID> set3 = getQuadsAroundVertex(v3);
    sofa::helper::vector<QuadID> set4 = getQuadsAroundVertex(v4);

    sort(set1.begin(), set1.end());
    sort(set2.begin(), set2.end());
    sort(set3.begin(), set3.end());
    sort(set4.begin(), set4.end());

    // The destination vector must be large enough to contain the result.
    sofa::helper::vector<QuadID> out1(set1.size()+set2.size());
    sofa::helper::vector<QuadID>::iterator result1;
    result1 = std::set_intersection(set1.begin(),set1.end(),set2.begin(),set2.end(),out1.begin());
    out1.erase(result1,out1.end());

    sofa::helper::vector<QuadID> out2(set3.size()+out1.size());
    sofa::helper::vector<QuadID>::iterator result2;
    result2 = std::set_intersection(set3.begin(),set3.end(),out1.begin(),out1.end(),out2.begin());
    out2.erase(result2,out2.end());

    sofa::helper::vector<QuadID> out3(set4.size()+out2.size());
    sofa::helper::vector<QuadID>::iterator result3;
    result3 = std::set_intersection(set4.begin(),set4.end(),out2.begin(),out2.end(),out3.begin());
    out3.erase(result3,out3.end());

    if (CHECK_TOPOLOGY)
        msg_warning_when(out3.size() > 1) << "More than one Quad found for indices: [" << v1 << "; " << v2 << "; " << v3 << "; " << v4 << "]";

    if(out3.size()==1)
        return (int) (out3[0]);

    return InvalidID;
}

size_t QuadSetTopologyContainer::getNumberOfQuads() const
{
    helper::ReadAccessor< Data< sofa::helper::vector<Quad> > > m_quad = d_quad;
    return m_quad.size();
}


size_t QuadSetTopologyContainer::getNumberOfElements() const
{
    return this->getNumberOfQuads();
}


const sofa::helper::vector< QuadSetTopologyContainer::QuadsAroundVertex > &QuadSetTopologyContainer::getQuadsAroundVertexArray()
{
    if(CHECK_TOPOLOGY)	// this method should only be called when the shell array exists
        msg_warning_when(!hasQuadsAroundVertex()) << "QuadsAroundVertex shell array is empty.";

    return m_quadsAroundVertex;
}

const sofa::helper::vector< QuadSetTopologyContainer::QuadsAroundEdge > &QuadSetTopologyContainer::getQuadsAroundEdgeArray()
{
    if(CHECK_TOPOLOGY)	// this method should only be called when the shell array exists
        msg_warning_when(!hasQuadsAroundEdge()) << "QuadsAroundEdge shell array is empty.";

    return m_quadsAroundEdge;
}

const sofa::helper::vector< QuadSetTopologyContainer::EdgesInQuad> &QuadSetTopologyContainer::getEdgesInQuadArray()
{
    if(CHECK_TOPOLOGY)
        msg_warning_when(m_edgesInQuad.empty()) << "EdgesInQuad shell array is empty.";

    return m_edgesInQuad;
}

const QuadSetTopologyContainer::QuadsAroundVertex& QuadSetTopologyContainer::getQuadsAroundVertex(PointID id)
{
    if (id < m_quadsAroundVertex.size())
        return m_quadsAroundVertex[id];
    else if (CHECK_TOPOLOGY)
        msg_error() << "QuadsAroundVertex array access out of bounds: " << id << " >= " << m_quadsAroundVertex.size();

    return InvalidSet;
}

const QuadSetTopologyContainer::QuadsAroundEdge& QuadSetTopologyContainer::getQuadsAroundEdge(EdgeID id)
{
    if (id < m_quadsAroundEdge.size())
        return m_quadsAroundEdge[id];
    else if (CHECK_TOPOLOGY)
        msg_error() << "QuadsAroundEdge array access out of bounds: " << id << " >= " << m_quadsAroundEdge.size();

    return InvalidSet;
}

const QuadSetTopologyContainer::EdgesInQuad &QuadSetTopologyContainer::getEdgesInQuad(QuadID id)
{
    if (id < m_edgesInQuad.size())
        return m_edgesInQuad[id];
    else if (CHECK_TOPOLOGY)
        msg_error() << "EdgesInQuad array access out of bounds: " << id << " >= " << m_edgesInQuad.size();

    return InvalidQuad;
}

int QuadSetTopologyContainer::getVertexIndexInQuad(const Quad &t, PointID vertexIndex) const
{
    if(t[0]==vertexIndex)
        return 0;
    else if(t[1]==vertexIndex)
        return 1;
    else if(t[2]==vertexIndex)
        return 2;
    else if(t[3]==vertexIndex)
        return 3;
    else
        return -1;
}

int QuadSetTopologyContainer::getEdgeIndexInQuad(const EdgesInQuad &t, EdgeID edgeIndex) const
{
    if(t[0]==edgeIndex)
        return 0;
    else if(t[1]==edgeIndex)
        return 1;
    else if(t[2]==edgeIndex)
        return 2;
    else if(t[3]==edgeIndex)
        return 3;
    else
        return -1;
}

QuadSetTopologyContainer::QuadsAroundEdge &QuadSetTopologyContainer::getQuadsAroundEdgeForModification(const EdgeID i)
{
    if(!hasQuadsAroundEdge())	// this method should only be called when the shell array exists
    {
		if (CHECK_TOPOLOGY)
			msg_warning() << "Quad edge shell array is empty.";

        createQuadsAroundEdgeArray();
    }

    if( i >= m_quadsAroundEdge.size())
    {
		if (CHECK_TOPOLOGY)
			msg_error() << "Index out of bounds.";

        createQuadsAroundEdgeArray();
    }

    return m_quadsAroundEdge[i];
}

QuadSetTopologyContainer::QuadsAroundVertex &QuadSetTopologyContainer::getQuadsAroundVertexForModification(const PointID i)
{
    if(!hasQuadsAroundVertex())	// this method should only be called when the shell array exists
    {
		if (CHECK_TOPOLOGY)
			msg_warning() << "Quad vertex shell array is empty.";

        createQuadsAroundVertexArray();
    }

    if( i >= m_quadsAroundVertex.size())
    {
		if (CHECK_TOPOLOGY)
			msg_error() << "Index out of bounds.";

        createQuadsAroundVertexArray();
    }

    return m_quadsAroundVertex[i];
}


bool QuadSetTopologyContainer::checkTopology() const
{
	if (CHECK_TOPOLOGY)
	{
		bool ret = true;
		helper::ReadAccessor< Data< sofa::helper::vector<Quad> > > m_quad = d_quad;

		if (hasQuadsAroundVertex())
		{
			for (size_t i = 0; i < m_quadsAroundVertex.size(); ++i)
			{
				const sofa::helper::vector<QuadID> &tvs = m_quadsAroundVertex[i];
				for (size_t j = 0; j < tvs.size(); ++j)
				{
					if ((m_quad[tvs[j]][0] != i)
						&& (m_quad[tvs[j]][1] != i)
						&& (m_quad[tvs[j]][2] != i)
						&& (m_quad[tvs[j]][3] != i))
					{
						ret = false;
						msg_error() << "*** CHECK FAILED : check_quad_vertex_shell, i = " << i << " , j = " << j;
					}
				}
			}
		}

		if (hasQuadsAroundEdge())
		{
			for (size_t i = 0; i < m_quadsAroundEdge.size(); ++i)
			{
				const sofa::helper::vector<QuadID> &tes = m_quadsAroundEdge[i];
				for (size_t j = 0; j < tes.size(); ++j)
				{
					if ((m_edgesInQuad[tes[j]][0] != i)
						&& (m_edgesInQuad[tes[j]][1] != i)
						&& (m_edgesInQuad[tes[j]][2] != i)
						&& (m_edgesInQuad[tes[j]][3] != i))
					{
						ret = false;
						msg_error() << "*** CHECK FAILED : check_quad_edge_shell, i = " << i << " , j = " << j;
					}
				}
			}
		}

		return ret && EdgeSetTopologyContainer::checkTopology();
	}
	
	return true;
}



/// Get information about connexity of the mesh
/// @{

bool QuadSetTopologyContainer::checkConnexity()
{

    size_t nbr = this->getNbQuads();

    if (nbr == 0)
    {
		if (CHECK_TOPOLOGY)
			msg_error() << "Can't compute connexity as there are no Quads";

        return false;
    }

    VecQuadID elemAll = this->getConnectedElement(0);

    if (elemAll.size() != nbr)
    {
		msg_error() << "In computing connexity, Quads are missings. There is more than one connexe component.";
        return false;
    }

    return true;
}


size_t QuadSetTopologyContainer::getNumberOfConnectedComponent()
{
    size_t nbr = this->getNbQuads();

    if (nbr == 0)
    {
		if (CHECK_TOPOLOGY)
			msg_error() << "Can't compute connexity as there are no Quads";

        return 0;
    }

    VecQuadID elemAll = this->getConnectedElement(0);
    size_t cpt = 1;

    while (elemAll.size() < nbr)
    {
        std::sort(elemAll.begin(), elemAll.end());
        QuadID other_QuadID = (QuadID)elemAll.size();

        for (QuadID i = 0; i<elemAll.size(); ++i)
            if (elemAll[i] != i)
            {
                other_QuadID = i;
                break;
            }

        VecQuadID elemTmp = this->getConnectedElement(other_QuadID);
        cpt++;

        elemAll.insert(elemAll.begin(), elemTmp.begin(), elemTmp.end());
    }

    return cpt;
}


const QuadSetTopologyContainer::VecQuadID QuadSetTopologyContainer::getConnectedElement(QuadID elem)
{
    if(!hasQuadsAroundVertex())	// this method should only be called when the shell array exists
    {
		if (CHECK_TOPOLOGY)
			msg_warning() << "Quad vertex shell array is empty.";

        createQuadsAroundVertexArray();
    }

    VecQuadID elemAll;
    VecQuadID elemOnFront, elemPreviousFront, elemNextFront;
    bool end = false;
    size_t cpt = 0;
    size_t nbr = this->getNbQuads();

    // init algo
    elemAll.push_back(elem);
    elemOnFront.push_back(elem);
    elemPreviousFront.clear();
    cpt++;

    while (!end && cpt < nbr)
    {
        // First Step - Create new region
        elemNextFront = this->getElementAroundElements(elemOnFront); // for each QuadID on the propagation front

        // Second Step - Avoid backward direction
        for (size_t i = 0; i<elemNextFront.size(); ++i)
        {
            bool find = false;
            QuadID id = elemNextFront[i];

            for (size_t j = 0; j<elemAll.size(); ++j)
                if (id == elemAll[j])
                {
                    find = true;
                    break;
                }

            if (!find)
            {
                elemAll.push_back(id);
                elemPreviousFront.push_back(id);
            }
        }

        // cpt for connexity
        cpt +=elemPreviousFront.size();

        if (elemPreviousFront.empty())
        {
            end = true;
			if (CHECK_TOPOLOGY)
				msg_error() << "Loop for computing connexity has reach end.";

        }

        // iterate
        elemOnFront = elemPreviousFront;
        elemPreviousFront.clear();
    }

    return elemAll;
}


const QuadSetTopologyContainer::VecQuadID QuadSetTopologyContainer::getElementAroundElement(QuadID elem)
{
    VecQuadID elems;

    if (!hasQuadsAroundVertex())
    {
		if (CHECK_TOPOLOGY)
			msg_warning() << "Quad vertex shell array is empty.";

        createQuadsAroundVertexArray();
    }

    Quad the_quad = this->getQuad(elem);

    for(size_t i = 0; i<4; ++i) // for each node of the Quad
    {
        QuadsAroundVertex quadAV = this->getQuadsAroundVertex(the_quad[i]);

        for (size_t j = 0; j<quadAV.size(); ++j) // for each Quad around the node
        {
            bool find = false;
            QuadID id = quadAV[j];

            if (id == elem)
                continue;

            for (size_t k = 0; k<elems.size(); ++k) // check no redundancy
                if (id == elems[k])
                {
                    find = true;
                    break;
                }

            if (!find)
                elems.push_back(id);
        }
    }

    return elems;
}


const QuadSetTopologyContainer::VecQuadID QuadSetTopologyContainer::getElementAroundElements(VecQuadID elems)
{
    VecQuadID elemAll;
    VecQuadID elemTmp;

    if (!hasQuadsAroundVertex())
    {
		if (CHECK_TOPOLOGY)
			msg_warning() << "Quad vertex shell array is empty.";

        createQuadsAroundVertexArray();
    }

    for (size_t i = 0; i <elems.size(); ++i) // for each QuadId of input vector
    {
        VecQuadID elemTmp2 = this->getElementAroundElement(elems[i]);

        elemTmp.insert(elemTmp.end(), elemTmp2.begin(), elemTmp2.end());
    }

    for (size_t i = 0; i<elemTmp.size(); ++i) // for each Quad Id found
    {
        bool find = false;
        QuadID id = elemTmp[i];

        for (size_t j = 0; j<elems.size(); ++j) // check no redundancy with input vector
            if (id == elems[j])
            {
                find = true;
                break;
            }

        if (!find)
        {
            for (size_t j = 0; j<elemAll.size(); ++j) // check no redundancy in output vector
                if (id == elemAll[j])
                {
                    find = true;
                    break;
                }
        }

        if (!find)
            elemAll.push_back(id);
    }


    return elemAll;
}

/// @}




bool QuadSetTopologyContainer::hasQuads() const
{
    d_quad.updateIfDirty();
    return !(d_quad.getValue()).empty();
}

bool QuadSetTopologyContainer::hasEdgesInQuad() const
{
    return !m_edgesInQuad.empty();
}

bool QuadSetTopologyContainer::hasQuadsAroundVertex() const
{
    return !m_quadsAroundVertex.empty();
}

bool QuadSetTopologyContainer::hasQuadsAroundEdge() const
{
    return !m_quadsAroundEdge.empty();
}

void QuadSetTopologyContainer::clearQuadsAroundVertex()
{
    for(size_t i=0; i<m_quadsAroundVertex.size(); ++i)
        m_quadsAroundVertex[i].clear();

    m_quadsAroundVertex.clear();
}

void QuadSetTopologyContainer::clearQuadsAroundEdge()
{
    for(size_t i=0; i<m_quadsAroundEdge.size(); ++i)
        m_quadsAroundEdge[i].clear();

    m_quadsAroundEdge.clear();
}

void QuadSetTopologyContainer::clearEdgesInQuad()
{
    m_edgesInQuad.clear();
}

void QuadSetTopologyContainer::clearQuads()
{
    helper::WriteAccessor< Data< sofa::helper::vector<Quad> > > m_quad = d_quad;
    m_quad.clear();
}

void QuadSetTopologyContainer::clear()
{
    clearQuadsAroundVertex();
    clearQuadsAroundEdge();
    clearEdgesInQuad();
    clearQuads();

    EdgeSetTopologyContainer::clear();
}

void QuadSetTopologyContainer::setQuadTopologyToDirty()
{
    // set this container to dirty
    m_quadTopologyDirty = true;

    // set all engines link to this container to dirty
    std::list<sofa::core::topology::TopologyEngine *>::iterator it;
    for (it = m_enginesList.begin(); it!=m_enginesList.end(); ++it)
    {
        sofa::core::topology::TopologyEngine* topoEngine = (*it);
        topoEngine->setDirtyValue();
        if (CHECK_TOPOLOGY)
            msg_info() << "Quad Topology Set dirty engine: " << topoEngine->name;
    }
}

void QuadSetTopologyContainer::cleanQuadTopologyFromDirty()
{
    m_quadTopologyDirty = false;

    // security, clean all engines to avoid loops
    std::list<sofa::core::topology::TopologyEngine *>::iterator it;
    for ( it = m_enginesList.begin(); it!=m_enginesList.end(); ++it)
    {
        if ((*it)->isDirty())
        {
            if (CHECK_TOPOLOGY)
                msg_warning() << "Quad Topology update did not clean engine: " << (*it)->name;
            (*it)->cleanDirty();
        }
    }
}

void QuadSetTopologyContainer::updateTopologyEngineGraph()
{
    // calling real update Data graph function implemented once in PointSetTopologyModifier
    this->updateDataEngineGraph(this->d_quad, this->m_enginesList);

    // will concatenate with edges one:
    EdgeSetTopologyContainer::updateTopologyEngineGraph();
}

} // namespace topology

} // namespace component

} // namespace sofa

