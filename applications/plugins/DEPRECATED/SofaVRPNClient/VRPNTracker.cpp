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
 * VRPNTracker.cpp
 *
 *  Created on: 8 sept. 2009
 *      Author: froy
 */

#define SOFAVRPNCLIENT_VRPNTRACKER_CPP_

#include "VRPNTracker.inl"

#include <sofa/core/ObjectFactory.h>
#include <vrpnclient_config.h>

namespace sofavrpn
{

namespace client
{

using namespace sofa::defaulttype;
using namespace sofavrpn::client;

int VRPNTrackerClass = sofa::core::RegisterObject("VRPN Tracker")
        .add< VRPNTracker<Vec3Types> >()
        .add< VRPNTracker<Rigid3Types> >()
 
        ;

template class SOFA_SOFAVRPNCLIENT_API VRPNTracker<Vec3Types>;
template class SOFA_SOFAVRPNCLIENT_API VRPNTracker<Rigid3Types>;
 

}

}
