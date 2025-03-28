/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
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
#include <MultiThreading/ParallelImplementationsRegistry.h>
#include <MultiThreading/initMultiThreading.h>
#include <gtest/gtest.h>
#include <sofa/Modules.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/simpleapi/SimpleApi.h>
#include <sofa/testing/ScopedPlugin.h>

namespace multithreading
{

TEST(ParallelImplementationsRegistry, existInObjectFactory)
{
    // sequential versions will be added to the ObjectFactory
    const auto plugins = sofa::testing::makeScopedPlugin({
        Sofa.Component.LinearSolver.Iterative,
        Sofa.Component.Collision.Detection.Algorithm,
        Sofa.Component.SolidMechanics.FEM.Elastic,
        Sofa.Component.Mapping.Linear,
        "MultiThreading"
    });

    const auto implementations = ParallelImplementationsRegistry::getImplementations();

    for (const auto& [seq, par] : implementations)
    {
        ASSERT_FALSE(seq.empty());
        ASSERT_FALSE(par.empty());

        EXPECT_TRUE(sofa::core::ObjectFactory::getInstance()->hasCreator(seq)) << seq;
        EXPECT_TRUE(sofa::core::ObjectFactory::getInstance()->hasCreator(par)) << par;
    }
}

}
