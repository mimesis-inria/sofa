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
#ifndef SOFA_COMPONENT_LOADER_SPHERELOADER_H
#define SOFA_COMPONENT_LOADER_SPHERELOADER_H
#include "config.h"

#include <sofa/core/loader/BaseLoader.h>

namespace sofa
{
namespace component
{
namespace loader
{

class SphereLoader : public sofa::core::loader::BaseLoader
{
public:
    SOFA_CLASS(SphereLoader,sofa::core::loader::BaseLoader);
protected:
    SphereLoader();
public:
    // Point coordinates in 3D in double.
    Data< helper::vector<sofa::defaulttype::Vec<3,SReal> > > positions; ///< Sphere centers
    Data< helper::vector<SReal> > radius; ///< Radius of each sphere
    Data< defaulttype::Vector3 > d_scale; ///< Scale applied to sphere positions
    Data< defaulttype::Vector3 > d_translation; ///< Translation applied to sphere positions
    bool load() override;
};

} //loader
} //component
} //sofa

#endif // SOFA_COMPONENT_LOADER_SPHERELOADER_H
