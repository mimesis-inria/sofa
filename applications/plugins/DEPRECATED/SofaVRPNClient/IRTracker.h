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
 * IRTracker.h
 *
 *  Created on: 8 sept. 2009
 *      Author: froy
 */

#ifndef SOFAVRPNCLIENT_IRTRACKER_H_
#define SOFAVRPNCLIENT_IRTRACKER_H_

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/DataEngine.h>

#include <VRPNDevice.h>

#include <vrpn/vrpn_Analog.h>

namespace sofavrpn
{

namespace client
{

template<class DataTypes>
class IRTracker : public virtual sofa::core::objectmodel::BaseObject, public virtual sofa::core::DataEngine
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(IRTracker, DataTypes), sofa::core::objectmodel::BaseObject);

    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Point;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::VecCoord VecCoord;

    static const double WIIMOTE_X_RESOLUTION;
    static const double WIIMOTE_Y_RESOLUTION;
    static const double WIIMOTE_X_ANGLE;
    static const double WIIMOTE_Y_ANGLE;

    IRTracker();
    virtual ~IRTracker();


    //input
    sofa::core::objectmodel::Data<VecCoord > f_leftDots; ///< IR dots from L camera
    sofa::core::objectmodel::Data<VecCoord > f_rightDots; ///< IR dots from R camera
    sofa::core::objectmodel::Data<Real> f_distance; ///< Distance between the 2 cameras
    sofa::core::objectmodel::Data<Real> f_distanceSide; ///< Distance of a side
    sofa::core::objectmodel::Data<Real> f_scale; ///< Scale

    //output
    sofa::core::objectmodel::Data<VecCoord> f_points; ///< Computed 3D Points

    //Parameters
    sofa::core::objectmodel::Data<double> p_yErrorCoeff; ///< Y Error Coefficient
    sofa::core::objectmodel::Data<double> p_sideErrorCoeff; ///< Side Error Coefficient
    sofa::core::objectmodel::Data<double> p_realSideErrorCoeff; ///< Real Side Error Coefficient

    void update();

    Coord get3DPoint(double lx, double ly, double rx, double ry);
private:

};

#if  !defined(SOFAVRPNCLIENT_IRTRACKER_CPP_)
extern template class SOFA_SOFAVRPNCLIENT_API IRTracker<defaulttype::Vec3Types>;
 
#endif

}

}

#endif /* SOFAVRPNCLIENT_IRTRACKER_H_ */
