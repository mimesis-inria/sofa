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
#ifndef SOFA_SIMULATION_SCENECHECKMISSINGREQUIREDPLUGIN_H
#define SOFA_SIMULATION_SCENECHECKMISSINGREQUIREDPLUGIN_H

#include "config.h"
#include "SceneCheck.h"
#include <map>
#include <vector>

namespace sofa {
namespace simulation {
        class Node;
} // namespace sofa
} // namespace simulation

namespace sofa
{
namespace simulation
{
namespace _scenechecking_
{

class SOFA_GRAPH_COMPONENT_API SceneCheckMissingRequiredPlugin : public SceneCheck
{
public:
    typedef std::shared_ptr<SceneCheckMissingRequiredPlugin> SPtr;
    static SPtr newSPtr() { return SPtr(new SceneCheckMissingRequiredPlugin()); }
    virtual const std::string getName() override;
    virtual const std::string getDesc() override;
    void doInit(Node* node) override;
    void doCheckOn(Node* node) override;
    void doPrintSummary() override;

private:    
    std::map<std::string, bool > m_loadedPlugins;
    std::map<std::string, std::vector<std::string> > m_requiredPlugins;
};

} // namespace _scenechecking_

namespace scenechecking
{
    using _scenechecking_::SceneCheckMissingRequiredPlugin;
}

} // namespace simulation
} // namespace sofa

#endif // SOFA_SIMULATION_SCENECHECKMISSINGREQUIREDPLUGIN_H
