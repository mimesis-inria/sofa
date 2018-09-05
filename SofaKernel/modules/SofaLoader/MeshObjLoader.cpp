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
#include <sofa/core/ObjectFactory.h>
#include <SofaLoader/MeshObjLoader.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/system/SetDirectory.h>
#include <fstream>

namespace sofa
{

namespace component
{

namespace loader
{

using namespace sofa::defaulttype;
using namespace sofa::core::loader;
using namespace sofa::helper::types;

SOFA_DECL_CLASS(MeshObjLoader)

int MeshObjLoaderClass = core::RegisterObject("Specific mesh loader for Obj file format.")
        .add< MeshObjLoader >()
        ;



MeshObjLoader::MeshObjLoader()
    : MeshLoader()
    , d_handleSeams(initData(&d_handleSeams, (bool)false, "handleSeams", "Preserve UV and normal seams information (vertices with multiple UV and/or normals)"))
    , loadMaterial(initData(&loadMaterial, (bool) true, "loadMaterial", "Load the related MTL file or use a default one?"))    
    , faceType(MeshObjLoader::TRIANGLE)
    , d_material(initData(&d_material,"material","Default material") )
    , materials(initData(&materials,"materials","List of materials") )
    , faceList(initData(&faceList,"faceList","List of face definitions.") )
    , texIndexList(initData(&texIndexList,"texcoordsIndex","Indices of textures coordinates used in faces definition."))
    , positionsList(initData(&positionsList,"positionsDefinition", "Vertex positions definition"))
    , texCoordsList(initData(&texCoordsList,"texcoordsDefinition", "Texture coordinates definition"))
    , normalsIndexList(initData(&normalsIndexList,"normalsIndex","List of normals of elements of the mesh loaded."))
    , normalsList(initData(&normalsList,"normalsDefinition","Normals definition"))
    , texCoords(initData(&texCoords,"texcoords","Texture coordinates of all faces, to be used as the parent data of a VisualModel texcoords data"))
    , computeMaterialFaces(initData(&computeMaterialFaces, false, "computeMaterialFaces", "True to activate export of Data instances containing list of face indices for each material"))
    , d_vertPosIdx      (initData   (&d_vertPosIdx, "vertPosIdx", "If vertices have multiple normals/texcoords stores vertices position indices"))
    , d_vertNormIdx     (initData   (&d_vertNormIdx, "vertNormIdx", "If vertices have multiple normals/texcoords stores vertices normal indices"))
{
    faceList.setGroup("OBJ");
    texIndexList.setGroup("OBJ");
    texCoordsList.setGroup("OBJ");
    normalsIndexList.setGroup("OBJ");
    normalsList.setGroup("OBJ");
    positionsList.setGroup("OBJ");
    //BUGFIX: data loaded from OBJ file should not be saved to XML
    faceList.setPersistent(false);
    texIndexList.setPersistent(false);
    texCoordsList.setPersistent(false);
    normalsIndexList.setPersistent(false);
    normalsList.setPersistent(false);
    positionsList.setPersistent(false);
    texCoords.setPersistent(false);
    d_positions.setPersistent(false);
    d_normals.setPersistent(false);
    d_edges.setPersistent(false);
    d_triangles.setPersistent(false);
    d_quads.setPersistent(false);
    d_edgesGroups.setPersistent(false);
    d_trianglesGroups.setPersistent(false);
    d_quadsGroups.setPersistent(false);
    texCoords.setPersistent(false);
    d_vertPosIdx.setPersistent(false);
    d_vertNormIdx.setPersistent(false);
}


MeshObjLoader::~MeshObjLoader()
{

}


bool MeshObjLoader::load()
{
    dmsg_info() << "Loading OBJ file: " << m_filename;

    bool fileRead = false;

    // -- Loading file
    const char* filename = m_filename.getFullPath().c_str();
    std::ifstream file(filename);

    if (!file.good())
    {
        msg_error() << "Error: MeshObjLoader: Cannot read file '" << m_filename << "'.";
        return false;
    }

    // -- Reading file
    fileRead = this->readOBJ (file,filename);
    file.close();

    return fileRead;
}


void MeshObjLoader::addGroup (const PrimitiveGroup& g)
{
    helper::vector< PrimitiveGroup>& my_edgesGroups = *(d_edgesGroups.beginEdit());
    helper::vector< PrimitiveGroup>& my_trianglesGroups = *(d_trianglesGroups.beginEdit());
    helper::vector< PrimitiveGroup>& my_quadsGroups = *(d_quadsGroups.beginEdit());

    switch (faceType)
    {
    case MeshObjLoader::EDGE:
        my_edgesGroups.push_back(g);
        break;
    case MeshObjLoader::TRIANGLE:
        my_trianglesGroups.push_back(g);
        break;
    case MeshObjLoader::QUAD:
        my_quadsGroups.push_back(g);
        break;
    default: break;
    }

    d_edgesGroups.endEdit();
    d_trianglesGroups.endEdit();
    d_quadsGroups.endEdit();
}

bool MeshObjLoader::readOBJ (std::ifstream &file, const char* filename)
{
 
    const bool handleSeams = d_handleSeams.getValue();
    helper::vector<sofa::defaulttype::Vector3>& my_positions = *(d_positions.beginEdit());
    helper::vector<sofa::defaulttype::Vector2>& my_texCoords = *(texCoordsList.beginEdit());
    helper::vector<sofa::defaulttype::Vector3>& my_normals   = *(normalsList.beginEdit());

    Material& material = *(d_material.beginEdit());
    helper::vector<Material>& my_materials = *(materials.beginEdit());
    helper::SVector< helper::SVector <int> >& my_faceList = *(faceList.beginEdit() );
    helper::SVector< helper::SVector <int> >& my_normalsList = *(normalsIndexList.beginEdit());
    helper::SVector< helper::SVector <int> >& my_texturesList   = *(texIndexList.beginEdit());
    helper::vector<int> nodes, nIndices, tIndices;

    helper::vector<Edge >& my_edges = *(d_edges.beginEdit());
    helper::vector<Triangle >& my_triangles = *(d_triangles.beginEdit());
    helper::vector<Quad >& my_quads = *(d_quads.beginEdit());

    //BUGFIX: clear pre-existing data before loading the file
    my_positions.clear();
    material.activated = false;
    my_texCoords.clear();
    my_normals.clear();
    my_materials.clear();
    my_faceList.clear();
    my_normalsList.clear();
    my_texturesList.clear();
    my_edges.clear();
    my_triangles.clear();
    my_quads.clear();
    d_edgesGroups.beginEdit()->clear(); d_edgesGroups.endEdit();
    d_trianglesGroups.beginEdit()->clear(); d_trianglesGroups.endEdit();
    d_quadsGroups.beginEdit()->clear(); d_quadsGroups.endEdit();

    int vtn[3];
    Vector3 result;
    helper::WriteAccessor<Data<helper::vector< PrimitiveGroup> > > my_faceGroups[NBFACETYPE] =
    {
        d_edgesGroups,
        d_trianglesGroups,
        d_quadsGroups
    };
    std::string curGroupName = "Default_Group";
    std::string curMaterialName;
    int curMaterialId = -1;
    int nbFaces[NBFACETYPE] = {0}; // number of edges, triangles, quads
    int groupF0[NBFACETYPE] = {0}; // first primitives indices in current group for edges, triangles, quads
    std::string line;
    while( std::getline(file,line) )
    {
        if (line.empty()) continue;
        std::istringstream values(line);
        std::string token;
        values >> token;

        if (token == "#")
        {
            // comment
        }
        else if (token == "v")
        {
            // vertex
            values >> result[0] >> result[1] >> result[2];
            my_positions.push_back(Vector3(result[0],result[1], result[2]));
        }
        else if (token == "vn")
        {
            // normal 
            values >> result[0] >> result[1] >> result[2];
            my_normals.push_back(Vector3(result[0],result[1], result[2]));
        }
        else if (token == "vt")
        {
            // texcoord
            values >> result[0] >> result[1];
            my_texCoords.push_back(Vector2(result[0],result[1]));
        }
        else if ((token == "mtllib") && loadMaterial.getValue())
        {
            while (!values.eof())
            {
                std::string materialLibaryName;
                values >> materialLibaryName;
                std::string mtlfile = sofa::helper::system::SetDirectory::GetRelativeFromFile(materialLibaryName.c_str(), filename);
                this->readMTL(mtlfile.c_str(), my_materials);
            }
        }
        else if (token == "usemtl" || token == "g")
        {
            // end of current group
            //curGroup.nbp = nbf - curGroup.p0;
            for (int ft = 0; ft < NBFACETYPE; ++ft)
                if (nbFaces[ft] > groupF0[ft])
                {
                    my_faceGroups[ft].push_back(PrimitiveGroup(groupF0[ft], nbFaces[ft]-groupF0[ft], curMaterialName, curGroupName, curMaterialId));
                    groupF0[ft] = nbFaces[ft];
                }
            if (token == "usemtl")
            {
                values >> curMaterialName;
                curMaterialId = -1;
                helper::vector<Material>::iterator it = my_materials.begin();
                helper::vector<Material>::iterator itEnd = my_materials.end();
                for (; it != itEnd; ++it)
                {
                    if (it->name == curMaterialName)
                    {
                        (*it).activated = true;
                        if (!material.activated)
                            material = *it;
                        curMaterialId = it - my_materials.begin();
                        break;
                    }
                }
            }
            else if (token == "g")
            {
                curGroupName.clear();
                while (!values.eof())
                {
                    std::string g;
                    values >> g;
                    if (!curGroupName.empty())
                        curGroupName += " ";
                    curGroupName += g;
                }
            }
        }
        else if (token == "l" || token == "f")
        {
            // face
            nodes.clear();
            nIndices.clear();
            tIndices.clear();

            while (!values.eof())
            {
                std::string face;
                values >> face;
                if (face.empty()) continue;
                for (int j = 0; j < 3; j++)
                {
                    vtn[j] = -1;
                    std::string::size_type pos = face.find('/');
                    std::string tmp = face.substr(0, pos);
                    if (pos == std::string::npos)
                        face = "";
                    else
                    {
                        face = face.substr(pos + 1);
                    }

                    if (!tmp.empty())
                    {
                        vtn[j] = atoi(tmp.c_str());
                        if (vtn[j] >= 1)
                            vtn[j] -=1; // -1 because the numerotation begins at 1 and a vector begins at 0
                        else if (vtn[j] < 0)
                            vtn[j] += (j==0) ? my_positions.size() : (j==1) ? my_texCoords.size() : my_normals.size();
                        else
                        {
                            msg_error() << "Invalid index " << tmp;
                            vtn[j] = -1;
                        }
                    }
                }

                nodes.push_back(vtn[0]);
                tIndices.push_back(vtn[1]);
                nIndices.push_back(vtn[2]);
            }

            my_faceList.push_back(nodes);
            my_normalsList.push_back(nIndices);
            my_texturesList.push_back(tIndices);

            if (nodes.size() == 2) // Edge
            {
                if (!handleSeams) // we have to wait for renumbering vertices if we handle seams
                {
                    if (nodes[0]<nodes[1])
                        addEdge(&my_edges, Edge(nodes[0], nodes[1]));
                    else
                        addEdge(&my_edges, Edge(nodes[1], nodes[0]));
                }
                ++nbFaces[MeshObjLoader::EDGE];
                faceType = MeshObjLoader::EDGE;
            }
            else if (nodes.size()==4 && !this->d_triangulate.getValue()) // Quad
            {
                if (!handleSeams) // we have to wait for renumbering vertices if we handle seams
                {
                    addQuad(&my_quads, Quad(nodes[0], nodes[1], nodes[2], nodes[3]));
                }
                ++nbFaces[MeshObjLoader::QUAD];
                faceType = MeshObjLoader::QUAD;
            }
            else // Triangulate
            {
                if (!handleSeams) // we have to wait for renumbering vertices if we handle seams
                {
                    for (size_t j=2; j<nodes.size(); j++)
                        addTriangle(&my_triangles, Triangle(nodes[0], nodes[j-1], nodes[j]));
                }
                ++nbFaces[MeshObjLoader::TRIANGLE];
                faceType = MeshObjLoader::TRIANGLE;
            }

        }
        else
        {
            // std::cerr << "readObj : Unknown token for line " << line << std::endl;
        }
    }

    // end of current group
    for (size_t ft = 0; ft < NBFACETYPE; ++ft)
        if (nbFaces[ft] > groupF0[ft])
        {
            my_faceGroups[ft].push_back(PrimitiveGroup(groupF0[ft], nbFaces[ft]-groupF0[ft], curMaterialName, curGroupName, curMaterialId));
            groupF0[ft] = nbFaces[ft];
        }

    if (!d_handleSeams.getValue())
    { // default mode, vertices are never duplicated, only one texcoord and normal is used per vertex
        helper::vector<sofa::defaulttype::Vector2>& vTexCoords = *texCoords.beginEdit();
        helper::vector<sofa::defaulttype::Vector3>& vNormals   = *d_normals.beginEdit();
        helper::vector<sofa::defaulttype::Vector3>& vVertices  = *d_positions.beginEdit();
        vVertices = my_positions;
        size_t vertexCount = my_positions.size();
        if( my_texCoords.size() > 0 )
        {
            vTexCoords.resize(vertexCount);
        }
        else
        {
            vTexCoords.resize(0);
        }
        if( my_normals.size() > 0 )
        {
            vNormals.resize(vertexCount);
        }
        else
        {
            vNormals.resize(0);
        }
        for (size_t fi=0; fi<my_faceList.size(); ++fi)
        {
            const helper::SVector<int>& nodes = my_faceList[fi];
            const helper::SVector<int>& nIndices = my_normalsList[fi];
            const helper::SVector<int>& tIndices = my_texturesList[fi];
            for (size_t i = 0; i < nodes.size(); ++i)
            {
                unsigned int pi = nodes[i];
                unsigned int ni = nIndices[i];
                unsigned int ti = tIndices[i];
                if (pi >= vertexCount) continue;
                if (ti < my_texCoords.size() && (vTexCoords[pi] == sofa::defaulttype::Vector2() ||
                                                 (my_texCoords[ti]-vTexCoords[pi])*sofa::defaulttype::Vector2(-1,1) > 0))
                    vTexCoords[pi] = my_texCoords[ti];
                if (ni < my_normals.size())
                    vNormals[pi] += my_normals[ni];
            }
        }
        for (size_t i=0; i<vNormals.size(); ++i)
        {
            vNormals[i].normalize();
        }
    }
    else
    { // handleSeam mode : vertices are duplicated in case they have different texcoords and/or normals
        // This code was initially in VisualModelImpl::setMesh()
        
        int nbVIn = (int)my_positions.size();
        // First we compute for each point how many pair of normal/texcoord indices are used
        // The map store the final index of each combinaison
        std::vector< std::map< std::pair<int,int>, int > > vertTexNormMap;
        vertTexNormMap.resize(nbVIn);
        for (size_t fi=0; fi<my_faceList.size(); ++fi)
        {
            const helper::SVector<int>& nodes = my_faceList[fi];
            const helper::SVector<int>& nIndices = my_normalsList[fi];
            const helper::SVector<int>& tIndices = my_texturesList[fi];
            for (size_t i = 0; i < nodes.size(); ++i)
            {
                unsigned int pi = nodes[i];
                unsigned int ni = nIndices[i];
                unsigned int ti = tIndices[i];
                vertTexNormMap[pi][std::make_pair(ti, ni)] = 0;
            }
        }

        // Then we can compute how many vertices are created
        int nbVOut = 0;
        bool vsplit = false;
        for (int i = 0; i < nbVIn; i++)
        {
            int s = vertTexNormMap[i].size();
            nbVOut += s;
        }

		dmsg_info() << nbVIn << " input positions, " << nbVOut << " final vertices.";

        if (nbVIn != nbVOut)
            vsplit = true;

        // Then we can create the final arrays
        helper::WriteAccessor<Data<helper::vector<sofa::defaulttype::Vector3> > > vertices2 = d_positions;
        helper::WriteAccessor<Data<helper::vector<sofa::defaulttype::Vector3> > > vnormals = d_normals;
        helper::WriteAccessor<Data<helper::vector<sofa::defaulttype::Vector2> > > vtexcoords = texCoords;
        helper::WriteAccessor<Data<helper::vector<int> > > vertPosIdx = d_vertPosIdx;
        helper::WriteAccessor<Data<helper::vector<int> > > vertNormIdx = d_vertNormIdx;

        vertices2.resize(nbVOut);
        vnormals.resize(nbVOut);
        vtexcoords.resize(nbVOut);
        if (vsplit)
        {
            vertPosIdx.resize(nbVOut);
            vertNormIdx.resize(nbVOut);
        }

        int nbNOut = 0; /// Number of different normals
        for (int i = 0, j = 0; i < nbVIn; i++)
        {
            std::map<int, int> normMap;
            for (std::map<std::pair<int, int>, int>::iterator it = vertTexNormMap[i].begin();
                 it != vertTexNormMap[i].end(); ++it)
            {
                int t = it->first.first;
                int n = it->first.second;
                if ( (unsigned)n < my_normals.size())
                    vnormals[j] = my_normals[n];
                if ((unsigned)t < my_texCoords.size())
                    vtexcoords[j] = my_texCoords[t];
                
                vertices2[j] = my_positions[i];
                if (vsplit)
                {
                    vertPosIdx[j] = i;
                    if (normMap.count(n))
                        vertNormIdx[j] = normMap[n];
                    else
                    {
                        vertNormIdx[j] = nbNOut;
                        normMap[n] = nbNOut++;
                    }
                }
                it->second = j++;
            }
        }

        if( vsplit && nbNOut == nbVOut )
            vertNormIdx.resize(0);

        // Then we create the triangles and quads
        
        for (size_t fi=0; fi<my_faceList.size(); ++fi)
        {
            const helper::SVector<int>& verts = my_faceList[fi];
            const helper::SVector<int>& nIndices = my_normalsList[fi];
            const helper::SVector<int>& tIndices = my_texturesList[fi];
            std::vector<int> nodes;
            nodes.resize(verts.size());
            for (size_t i = 0; i < verts.size(); ++i)
            {
                unsigned int pi = verts[i];
                unsigned int ni = nIndices[i];
                unsigned int ti = tIndices[i];
                nodes[i] = vertTexNormMap[pi][std::make_pair(ti, ni)];
                if ((unsigned)nodes[i] >= (unsigned)nbVOut)
                {
                    msg_error() << this->getPathName()<<" index "<<nodes[i]<<" out of range";
                    nodes[i] = 0;
                }
            }

            if (nodes.size() == 2) // Edge
            {
                if (nodes[0]<nodes[1])
                    addEdge(&my_edges, Edge(nodes[0], nodes[1]));
                else
                    addEdge(&my_edges, Edge(nodes[1], nodes[0]));
            }
            else if (nodes.size()==4 && !this->d_triangulate.getValue()) // Quad
            {
                addQuad(&my_quads, Quad(nodes[0], nodes[1], nodes[2], nodes[3]));
            }
            else // Triangulate
            {
                for (size_t j=2; j<nodes.size(); j++)
                    addTriangle(&my_triangles, Triangle(nodes[0], nodes[j-1], nodes[j]));
            }
        }
        for (size_t i=0; i<vnormals.size(); ++i)
        {
            vnormals[i].normalize();
        }
    }


    if (computeMaterialFaces.getValue())
    {
        // create subset lists
        std::map< std::string, helper::vector<unsigned int> > materialFaces[NBFACETYPE];
        for (int ft = 0; ft < NBFACETYPE; ++ft)
        {
            for (size_t gi=0; gi<my_faceGroups[ft].size(); ++gi)
            {
                PrimitiveGroup g = my_faceGroups[ft][gi];
                helper::vector<unsigned int>& out = materialFaces[ft][g.materialName];
                for (int f=g.p0; f<g.p0+g.nbp; ++f)
                    out.push_back(f);
            }
        }
        for (int ft = 0; ft < NBFACETYPE; ++ft)
        {
            std::string fname;
            switch (faceType)
            {
            case MeshObjLoader::EDGE:     fname = "edge"; break;
            case MeshObjLoader::TRIANGLE: fname = "triangle"; break;
            case MeshObjLoader::QUAD:     fname = "quad"; break;
            default: break;
            }
            for (std::map< std::string, helper::vector<unsigned int> >::const_iterator it = materialFaces[ft].begin(), itend = materialFaces[ft].end(); it != itend; ++it)
            {
                std::string materialName = it->first;
                const helper::vector<unsigned>& faces = it->second;
                if (faces.empty()) continue;
                std::ostringstream oname;
                oname << "material_" << materialName << "_" << fname << "Indices";
                Data< helper::vector<unsigned int> >* dOut = new Data< helper::vector<unsigned int> >("list of face indices corresponding to a given material");
                dOut->setName(oname.str());

                this->addData(dOut);
                dOut->setGroup("Materials");
                dOut->setValue(faces);
                subsets_indices.push_back(dOut);
            }
        }
    }

    d_edgesGroups.endEdit();
    d_trianglesGroups.endEdit();
    d_quadsGroups.endEdit();
    d_positions.endEdit();
    d_edges.endEdit();
    d_triangles.endEdit();
    d_quads.endEdit();
    normalsList.endEdit();
    normalsIndexList.endEdit();
    d_material.endEdit();
    materials.endEdit();
    texIndexList.endEdit();
    texCoordsList.endEdit();
    texCoords.endEdit();
    faceList.endEdit();
    //vertices.endEdit();
    d_normals.endEdit();
    return true;
}



// -----------------------------------------------------
// readMTL: read a wavefront material library file
//
//    model - properly initialized GLMmodel structure
//    name  - name of the material library
// -----------------------------------------------------
bool MeshObjLoader::readMTL(const char* filename, helper::vector <Material>& materials)
{
    FILE* file;
    char buf[128]; // Note: in the strings below, 127 is sizeof(buf)-1
    const char *single_string_format = "%127s"; // Better than "%s" for scanf
    const char *double_string_format = "%127s %127s"; // Better than "%s %s"

    file = fopen(filename, "r");
    Material *mat = NULL;
    if (!file);//serr << "readMTL() failed: can't open material file " << filename << sendl;
    else
    {
        /* now, read in the data */
        while (fscanf(file, single_string_format, buf) != EOF)
        {

            switch (buf[0])
            {
            case '#':
                /* comment */
                /* eat up rest of line */
                if ( fgets(buf, sizeof(buf), file) == NULL)
                {
					if (feof(file))
						msg_error() << "Error: MeshObjLoader: fgets function has encounter end of file. case #.";
                    else
                        msg_error() << "Error: MeshObjLoader: fgets function has encounter an error. case #.";
                }
                break;
            case 'n':
                /* newmtl */
                if (mat != NULL)
                {
                    materials.push_back(*mat);
                    delete mat;
                    mat = NULL;
                }
                mat = new Material();
                if ( fgets(buf, sizeof(buf), file) == NULL)
                {
                    if (feof (file) )
                        msg_error() << "Error: MeshObjLoader: fgets function has encounter end of file. case n.";
                    else
                        msg_error() << "Error: MeshObjLoader: fgets function has encounter an error. case n.";
                }
                sscanf(buf, double_string_format, buf, buf);
                mat->name = buf;
                break;
            case 'N':
                switch (buf[1])
                {
                case 'i':
                {
                    float optical_density;
                    if ( fscanf(file, "%f", &optical_density) == EOF )
                        msg_error() << "Error: MeshObjLoader: fscanf function has encounter an error. case N i.";
                    break;
                }
                case 's':
                    if (fscanf(file, "%f", &mat->shininess) == EOF )
                        msg_error() << "Error: MeshObjLoader: fscanf function has encounter an error. case N s.";

                    mat->useShininess = true;
                    break;
                default:
                    /* eat up rest of line */
                    if ( fgets(buf, sizeof(buf), file) == NULL)
                    {
                        if (feof (file) )
                            msg_error() << "Error: MeshObjLoader: fgets function has encounter end of file. case N.";
                        else
                            msg_error() << "Error: MeshObjLoader: fgets function has encounter an error. case N.";
                    }
                    break;
                }
                break;
            case 'K':
                switch (buf[1])
                {
                case 'd':
                    if ( fscanf(file, "%f %f %f", &mat->diffuse[0], &mat->diffuse[1], &mat->diffuse[2]) == EOF)
                        msg_error() << "Error: MeshObjLoader: fscanf function has encounter an error. case K d.";
                    mat->useDiffuse = true;
                    break;
                case 's':
                    if ( fscanf(file, "%f %f %f", &mat->specular[0], &mat->specular[1], &mat->specular[2]) == EOF)
                        msg_error() << "Error: MeshObjLoader: fscanf function has encounter an error. case K s.";
                    mat->useSpecular = true;
                    break;
                case 'a':
                    if ( fscanf(file, "%f %f %f", &mat->ambient[0], &mat->ambient[1], &mat->ambient[2]) == EOF)
                        msg_error() << "Error: MeshObjLoader: fscanf function has encounter an error. case K a.";
                    mat->useAmbient = true;
                    break;
                default:
                    /* eat up rest of line */
                    if ( fgets(buf, sizeof(buf), file) == NULL)
                    {
                        if (feof (file) )
                            msg_error() << "Error: MeshObjLoader: fgets function has encounter end of file. case K.";
                        else
                            msg_error() << "Error: MeshObjLoader: fgets function has encounter an error. case K.";
                    }
                    break;
                }
                break;
            case 'd':
            case 'T':
				if (!mat)
				{
					msg_error("MeshOBJ") << "readMTL 'T' no material";
					break;
				}
                // transparency value
                if ( fscanf(file, "%f", &mat->diffuse[3]) == EOF)
                    msg_error() << "Error: MeshObjLoader: fscanf function has encounter an error. case T i.";
                break;

			case 'm':
			{
				if (!mat)
				{
					msg_error("MeshOBJ") << "readMTL 'm' no material";
					break;
				}
				//texture map
				char charFilename[128] = { 0 };
				if (fgets(charFilename, sizeof(charFilename), file) == NULL)
				{
					msg_error("MeshOBJ") << "fgets has encountered an error";
				}
				else
				{
					mat->useTexture = true;

					//store the filename of the texture map in the material
					std::string stringFilename(charFilename);

					//delete carriage return from the string assuming the next property of the .mtl file is at the next line
					stringFilename.erase(stringFilename.end() - 1, stringFilename.end());
					stringFilename.erase(stringFilename.begin(), stringFilename.begin() + 1);
					mat->textureFilename = stringFilename;
				}

				break;
			}
			case 'b':
			{
				if (!mat)
				{
					msg_error("MeshOBJ") << "readMTL 'b' no material";
					break;
				}
				//bump mapping texture map
				char charFilename[128] = { 0 };
				if (fgets(charFilename, sizeof(charFilename), file) == NULL)
				{
					msg_error("MeshOBJ") << "fgets has encountered an error";
				}
				else
				{
					mat->useBumpMapping = true;

					//store the filename of the texture map in the material
					std::string stringFilename(charFilename);

					//delete carriage return from the string assuming the next property of the .mtl file is at the next line
					stringFilename.erase(stringFilename.end() - 1, stringFilename.end());
					stringFilename.erase(stringFilename.begin(), stringFilename.begin() + 1);
					mat->bumpTextureFilename = stringFilename;
				}

				break;
			}
            default:
                /* eat up rest of line */
                if ( fgets(buf, sizeof(buf), file) == NULL)
                {
                    if (feof (file) )
                        msg_error() << "Error: MeshObjLoader: fgets function has encounter end of file. case default.";
                    else
                        msg_error() << "Error: MeshObjLoader: fgets function has encounter an error. case default.";
                }
                break;
            }

        }
        fclose(file);
    }
    if (mat != NULL)
    {
        materials.push_back(*mat);
        delete mat;
        mat = NULL;
    }

    return true;
}



} // namespace loader

} // namespace component

} // namespace sofa

