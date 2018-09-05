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
#ifndef SOFA_COMPONENT_LOADER_MESHOBJLOADER_H
#define SOFA_COMPONENT_LOADER_MESHOBJLOADER_H
#include "config.h"

#include <sofa/core/loader/MeshLoader.h>
#include <sofa/helper/SVector.h>
#include <sofa/helper/types/Material.h>

namespace sofa
{

namespace component
{

namespace loader
{

class SOFA_LOADER_API MeshObjLoader : public sofa::core::loader::MeshLoader
{
public:
    enum FaceType { EDGE, TRIANGLE, QUAD, NBFACETYPE };

    SOFA_CLASS(MeshObjLoader,sofa::core::loader::MeshLoader);
protected:
    MeshObjLoader();
    virtual ~MeshObjLoader();

public:
    virtual bool load() override;

    template <class T>
    static bool canCreate ( T*& obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg )
    {
        return BaseLoader::canCreate (obj, context, arg);
    }

protected:
    bool readOBJ (std::ifstream &file, const char* filename);
    bool readMTL (const char* filename, helper::vector <sofa::helper::types::Material>& materials);
    void addGroup (const sofa::core::loader::PrimitiveGroup& g);

    Data<bool> d_handleSeams;
    Data<bool> loadMaterial;
    std::string textureName;
    FaceType faceType;

public:
    Data<sofa::helper::types::Material> d_material;
    Data <helper::vector <sofa::helper::types::Material> > materials;
    Data <helper::SVector <helper::SVector <int> > > faceList;
    Data <helper::SVector <helper::SVector <int> > > texIndexList;
    Data <helper::vector<sofa::defaulttype::Vector3> > positionsList;
    Data< helper::vector<sofa::defaulttype::Vector2> > texCoordsList;
    Data <helper::SVector<helper::SVector<int> > > normalsIndexList;
    Data <helper::vector<sofa::defaulttype::Vector3> > normalsList;
    Data< helper::vector<sofa::defaulttype::Vector2> > texCoords;
    Data< bool > computeMaterialFaces;
    helper::vector< Data <helper::vector <unsigned int> >* > subsets_indices;

    /// If vertices have multiple normals/texcoords, then we need to separate them
    /// This vector store which input position is used for each vertex
    /// If it is empty then each vertex correspond to one position
    Data< helper::vector<int> > d_vertPosIdx;

    /// Similarly this vector store which input normal is used for each vertex
    /// If it is empty then each vertex correspond to one normal
    Data< helper::vector<int> > d_vertNormIdx;

    virtual std::string type() { return "The format of this mesh is OBJ."; }
};




} // namespace loader

} // namespace component

} // namespace sofa

#endif
