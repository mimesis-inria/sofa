/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2018 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU General Public License as published by the Free  *
* Software Foundation; either version 2 of the License, or (at your option)   *
* any later version.                                                          *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    *
* more details.                                                               *
*                                                                             *
* You should have received a copy of the GNU General Public License along     *
* with this program. If not, see <http://www.gnu.org/licenses/>.              *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
/*
 * ToolTracker.h
 *
 *  Created on: 8 sept. 2009
 *      Author: froy
 */

#ifndef SOFAVRPNCLIENT_TOOLTRACKER_H_
#define SOFAVRPNCLIENT_TOOLTRACKER_H_

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/DataEngine.h>
#include <sofa/defaulttype/Quat.h>
#include <VRPNDevice.h>

#include <vrpn/vrpn_Analog.h>

namespace sofavrpn
{

namespace client
{

/*
 * Find a tool given various parameters...
 *
 */

template<class DataTypes>
class ToolTracker : public virtual sofa::core::objectmodel::BaseObject, public virtual sofa::core::DataEngine
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(ToolTracker, DataTypes), sofa::core::objectmodel::BaseObject);

    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Point;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::VecCoord VecCoord;

    typedef typename sofa::defaulttype::RigidTypes::Coord RPoint;
    typedef typename sofa::defaulttype::RigidTypes::Coord RCoord;

    //input
    sofa::core::objectmodel::Data<VecCoord > f_points; ///< Incoming 3D Points
    //distances between each point for the given tool
    sofa::core::objectmodel::Data<sofa::helper::vector<double> > f_distances; ///< Distances between each point


    //output
    sofa::core::objectmodel::Data<Coord> f_center; ///< Tool's center
    sofa::core::objectmodel::Data<sofa::defaulttype::Quat> f_orientation; ///< Tool's orientation
    //the same...
    sofa::core::objectmodel::Data<RCoord> f_rigidCenter; ///< Rigid center of the tool

    //parameters
    sofa::core::objectmodel::Data<bool> f_drawTool; ///< Draw tool's contour

    ToolTracker();
    virtual ~ToolTracker();

//	void init();
//	void reinit();
    void update();
    void draw();

private:

};

#if  !defined(SOFAVRPNCLIENT_TOOLTRACKER_CPP_)
extern template class SOFA_SOFAVRPNCLIENT_API ToolTracker<defaulttype::Vec3Types>;
 
#endif

}

}

#endif /* SOFAVRPNCLIENT_TOOLTRACKER_H_ */
