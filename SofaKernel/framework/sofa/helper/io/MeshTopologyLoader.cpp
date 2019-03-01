/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2018 INRIA, USTL, UJF, CNRS, MGH                    *
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
#include <sofa/helper/io/MeshTopologyLoader.h>
#include <sofa/helper/system/FileRepository.h>
#include <sofa/defaulttype/Vec.h>
#include <string.h>

#if defined(WIN32)
#define strcasecmp stricmp
#endif


MSG_REGISTER_CLASS(sofa::helper::io::MeshTopologyLoader, "MeshTopologyLoader")

namespace sofa
{

namespace helper
{

namespace io
{

using namespace sofa::defaulttype;

bool MeshTopologyLoader::addMeshtoTopology()
{
    if (m_mesh == NULL)
        return false;

    setNbPoints((int)m_mesh->getVertices().size());

    const sofa::helper::vector<Vector3>& vertices = m_mesh->getVertices();
    const sofa::helper::vector< Topology::Edge > & edges = m_mesh->getEdges();
    const sofa::helper::vector< Topology::Triangle > & triangles = m_mesh->getTriangles();
    const sofa::helper::vector< Topology::Quad > & quads = m_mesh->getQuads();
    const sofa::helper::vector< Topology::Tetrahedron > & tetra = m_mesh->getTetrahedra();
    const sofa::helper::vector< Topology::Hexahedron > & hexa = m_mesh->getHexahedra();

    for (size_t i = 0; i < vertices.size(); ++i)
        addPoint(vertices[i][0], vertices[i][1], vertices[i][2]);

    for (size_t i = 0; i < edges.size(); ++i)
        addLine(edges[i][0], edges[i][1]);

    for (size_t i = 0; i < triangles.size(); ++i)
        addTriangle(triangles[i][0], triangles[i][1], triangles[i][2]);

    for (size_t i = 0; i < quads.size(); ++i)
        addQuad(quads[i][0], quads[i][1], quads[i][2], quads[i][3]);

    for (size_t i = 0; i < tetra.size(); ++i)
        addTetra(tetra[i][0], tetra[i][1], tetra[i][2], tetra[i][3]);

    for (size_t i = 0; i < hexa.size(); ++i)
        addCube(hexa[i][0], hexa[i][1], hexa[i][2], hexa[i][3],
            hexa[i][4], hexa[i][5], hexa[i][6], hexa[i][7]);

    return true;
}

bool MeshTopologyLoader::loadObj(const char *filename)
{
    m_mesh = helper::io::Mesh::Create(filename);
    if (m_mesh ==NULL)
        return false;

    setNbPoints((int)m_mesh->getVertices().size());
    for (size_t i=0; i<m_mesh->getVertices().size(); i++)
    {
        addPoint((SReal)m_mesh->getVertices()[i][0],
                (SReal)m_mesh->getVertices()[i][1],
                (SReal)m_mesh->getVertices()[i][2]);
    }

    const vector< vector < vector <int> > > & facets = m_mesh->getFacets();
    std::set< std::pair<int,int> > edges;
    for (size_t i=0; i<facets.size(); i++)
    {
        const vector<int>& facet = facets[i][0];
        if (facet.size()==2)
        {
            // Line
            if (facet[0]<facet[1])
                addLine(facet[0],facet[1]);
            else
                addLine(facet[1],facet[0]);
        }
        else if (facet.size()==4)
        {
            // Quad
            addQuad(facet[0],facet[1],facet[2],facet[3]);
        }
        else
        {
            // Triangularize
            for (size_t j=2; j<facet.size(); j++)
                addTriangle(facet[0],facet[j-1],facet[j]);
        }
#if 0
        // Add edges
        if (facet.size()>2)
        {
            for (size_t j=0; j<facet.size(); j++)
            {
                int i1 = facet[j];
                int i2 = facet[(j+1)%facet.size()];
                if (edges.count(std::make_pair(i1,i2))!=0)
                {
                }
                else if (edges.count(std::make_pair(i2,i1))==0)
                {
                    if (i1>i2)
                        addLine(i1,i2);
                    else
                        addLine(i2,i1);
                    edges.insert(std::make_pair(i1,i2));
                }
            }
        }
#endif
    }

    /// delete m_mesh;
    return true;
}

bool MeshTopologyLoader::loadGmsh(const char *filename)
{
    m_mesh = helper::io::Mesh::Create("gmsh", filename);      
    return addMeshtoTopology();
}

bool MeshTopologyLoader::loadMesh(std::ifstream &file)
{
    return false;

    std::string cmd;
    int npoints = 0;
    int nlines = 0;
    int ntris = 0;
    int nquads = 0;
    int ntetrahedra = 0;
    int ncubes = 0;

    while (!file.eof())
    {
        file >> cmd;
        if (cmd=="line")
        {
            int p1,p2;
            file >> p1 >> p2;
            addLine(p1, p2);
            ++nlines;
        }
        else if (cmd=="triangle")
        {
            int p1,p2,p3;
            file >> p1 >> p2 >> p3;
            addTriangle(p1, p2, p3);
            ++ntris;
        }
        else if (cmd=="quad")
        {
            int p1,p2,p3,p4;
            file >> p1 >> p2 >> p3 >> p4;
            addQuad(p1, p2, p3, p4);
            ++nquads;
        }
        else if (cmd=="tetra")
        {
            int p1,p2,p3,p4;
            file >> p1 >> p2 >> p3 >> p4;
            addTetra(p1, p2, p3, p4);
            ++ntetrahedra;
        }
        else if (cmd=="cube")
        {
            int p1,p2,p3,p4,p5,p6,p7,p8;
            file >> p1 >> p2 >> p3 >> p4 >> p5 >> p6 >> p7 >> p8;
            addCube(p1, p2, p3, p4, p5, p6, p7, p8);
            ++ncubes;
        }
        else if (cmd=="point")
        {
            double px,py,pz;
            file >> px >> py >> pz;
            addPoint(px, py, pz);
            ++npoints;
        }
        else if (cmd=="v")
        {
            double px,py,pz;
            file >> px >> py >> pz;
            addPoint(px, py, pz);
            ++npoints;
        }
        else if (cmd=="f")
        {
            int p1,p2,p3,p4=0;
            file >> p1 >> p2 >> p3 >> p4;
            if (p4)
            {
                addQuad(p1-1, p2-1, p3-1, p4-1);
                ++nquads;
            }
            else
            {
                addTriangle(p1-1, p2-1, p3-1);
                ++ntris;
            }
        }
        else if (cmd=="mass")
        {
            int index;
            char location;
            double px,py,pz,vx,vy,vz,mass=0.0,elastic=0.0;
            file >> index >> location >> px >> py >> pz >> vx >> vy >> vz >> mass >> elastic;
            addPoint(px, py, pz);
            ++npoints;
        }
        else if (cmd=="lspg")
        {
            int	index;
            int m1,m2;
            double ks=0.0,kd=0.0,initpos=-1;
            file >> index >> m1 >> m2 >> ks >> kd >> initpos;
            --m1;
            --m2;
            addLine(m1,m2);
            ++nlines;
        }
        else if (cmd[0] == '#')	// it's a comment
        {
        }
        else		// it's an unknown keyword
        {
            msg_error() << "Unknown Mesh keyword:" << cmd;
            return false;
        }
    }

    return true;
}


bool MeshTopologyLoader::loadVtk(const char *filename)
{
    m_mesh = helper::io::Mesh::Create("vtu", filename);
    return addMeshtoTopology();
}

bool MeshTopologyLoader::load(const char *filename)
{
	std::string fname(filename);
	if (!sofa::helper::system::DataRepository.findFile(fname))
	{
		msg_error() << "Cannot find file: " << filename;
		return false;
	}

	bool fileLoaded;

	// check the extension of the filename
	if ((strlen(filename) > 4 && !strcmp(filename + strlen(filename) - 4, ".obj"))
		|| (strlen(filename) > 6 && !strcmp(filename + strlen(filename) - 6, ".trian")))
		fileLoaded = loadObj(fname.c_str());
	else if (strlen(filename) > 4 && !strcmp(filename + strlen(filename) - 4, ".vtk"))
		fileLoaded = loadVtk(fname.c_str());
	else if (strlen(filename) > 9 && !strcmp(filename + strlen(filename) - 9, ".vtk_swap"))
		fileLoaded = loadVtk(fname.c_str());
    else if (strlen(filename) > 4 && !strcmp(filename + strlen(filename) - 4, ".msh"))
        fileLoaded = loadGmsh(fname.c_str());
	else
	{
		std::ifstream file(filename);
		if (!file.good()) return false;
		msg_error() << "This file format: " << filename << " will not be supported anymore in sofa release 18.06.";
		fileLoaded = loadMesh(file);
		file.close();
	}
       
    if(!fileLoaded)
        msg_error() << "Unable to load mesh file '" << fname << "'" ;

    return fileLoaded;
}

} // namespace io

} // namespace helper

} // namespace sofa

