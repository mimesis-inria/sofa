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
#ifndef SOFA_COMPONENT_ENGINE_EXTRUDEQUADSANDGENERATEHEXAS_H
#define SOFA_COMPONENT_ENGINE_EXTRUDEQUADSANDGENERATEHEXAS_H
#include "config.h"

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

#include <sofa/defaulttype/Vec.h>
#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/defaulttype/Vec3Types.h>

namespace sofa
{

namespace component
{

namespace engine
{

/**
 * This class extrudes a quad surface into a set of hexahedra
 */
template <class DataTypes>
class ExtrudeQuadsAndGenerateHexas : public core::DataEngine
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(ExtrudeQuadsAndGenerateHexas,DataTypes),core::DataEngine);

    typedef typename DataTypes::Coord     Coord;
    typedef typename DataTypes::VecCoord  VecCoord;
    typedef typename DataTypes::Real      Real;
    typedef defaulttype::Vec<3,Real>      Vec3;

protected:

    ExtrudeQuadsAndGenerateHexas();

    ~ExtrudeQuadsAndGenerateHexas() {}
public:
    void init() override;

    void reinit() override;

    void update() override;

    void draw( const core::visual::VisualParams* ) override;

    virtual std::string getTemplateName() const override
    {
        return templateName(this);
    }

    static std::string templateName(const ExtrudeQuadsAndGenerateHexas<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

    bool                                             initialized;
    Data<bool>                                       isVisible; ///< is Visible ?
    Data<Coord>                                      f_scale; ///< Apply a scaling factor to the extruded mesh
    Data<Real>                                       f_thickness;
    Data<Real>                                       f_thicknessIn; ///< Thickness of the extruded volume in the opposite direction of the normals
    Data<Real>                                       f_thicknessOut; ///< Thickness of the extruded volume in the direction of the normals
    Data<int>                                        f_numberOfSlices; ///< Number of slices / steps in the extrusion
    Data<VecCoord>                                   f_surfaceVertices; ///< Position coordinates of the surface
    Data< helper::vector<sofa::core::topology::BaseMeshTopology::Quad> >   f_surfaceQuads; ///< Indices of the quads of the surface to extrude
    Data<VecCoord>                                   f_extrudedVertices; ///< Coordinates of the extruded vertices
    Data< helper::vector<sofa::core::topology::BaseMeshTopology::Quad> >   f_extrudedSurfaceQuads; ///< List of new surface quads generated during the extrusion
    Data< helper::vector<sofa::core::topology::BaseMeshTopology::Quad> >   f_extrudedQuads; ///< List of all quads generated during the extrusion
    Data< helper::vector<sofa::core::topology::BaseMeshTopology::Hexa> >   f_extrudedHexas; ///< List of hexahedra generated during the extrusion
};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_ENGINE_EXTRUDEQUADSANDGENERATEHEXAS_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_GENERAL_ENGINE_API ExtrudeQuadsAndGenerateHexas<defaulttype::Vec3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
extern template class SOFA_GENERAL_ENGINE_API ExtrudeQuadsAndGenerateHexas<defaulttype::Vec3fTypes>;
#endif //SOFA_DOUBLE
#endif

} // namespace engine

} // namespace component

} // namespace sofa

#endif
