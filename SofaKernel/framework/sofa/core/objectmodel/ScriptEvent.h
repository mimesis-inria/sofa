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
#ifndef SOFA_CORE_OBJECTMODEL_SCRIPTEVENT_H
#define SOFA_CORE_OBJECTMODEL_SCRIPTEVENT_H

#include <sofa/core/objectmodel/Event.h>
#include <string>
#include <sofa/simulation/Node.h>


namespace sofa
{

namespace core
{

namespace objectmodel
{

/**
* @brief Generic Event class to send a message through the simulation graph.
*/
class SOFA_CORE_API ScriptEvent : public sofa::core::objectmodel::Event
{
public:

    SOFA_EVENT_H( ScriptEvent )

    /**
     * @brief Constructor.
     */
    ScriptEvent(sofa::simulation::Node::SPtr sender, const char* eventName);

    /**
     * @brief Destructor.
     */
    ~ScriptEvent() override;

    /**
     * @brief Get the sender name
     */
    const sofa::simulation::Node::SPtr getSender(void) const {return m_sender;}

    /**
     * @brief Get the event name
     */
    const std::string getEventName(void) const {return m_eventName;}

    const char* getClassName() const override { return "ScriptEvent"; }
private:

    sofa::simulation::Node::SPtr m_sender;
    std::string m_eventName;

};

} // namespace objectmodel

} // namespace core

} // namespace sofa

#endif // SOFA_CORE_OBJECTMODEL_SCRIPTEVENT_H
