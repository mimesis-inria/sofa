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
#include <SofaHaptics/ForceFeedback.h>

namespace sofa
{

namespace component
{

namespace controller
{

ForceFeedback::ForceFeedback():
    d_activate(initData(&d_activate, false, "activate", "boolean to activate or deactivate the forcefeedback"))
  , d_indice(initData(&d_indice, 0, "indice", "Tool indice in the OmniDriver"))
{
}

void ForceFeedback::init()
{
    context = dynamic_cast<simulation::Node *>(this->getContext());
}

void ForceFeedback::setReferencePosition(sofa::defaulttype::SolidTypes<SReal>::Transform& referencePosition)
{
    SOFA_UNUSED(referencePosition);
}

bool ForceFeedback::isEnabled() {
    return this->getContext()->isActive();
}

} // namespace controller

} // namespace component

} // namespace sofa
