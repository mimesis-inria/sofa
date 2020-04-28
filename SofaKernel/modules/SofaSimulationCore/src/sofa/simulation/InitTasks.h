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
#ifndef InitTasks_h__
#define InitTasks_h__

#include <sofa/simulation/TaskScheduler.h>

namespace sofa
{
    namespace simulation
    {
        
        using namespace sofa;
        
        
        
        class SOFA_SIMULATION_CORE_API InitPerThreadDataTask : public CpuTask
        {
            
        public:
            
            InitPerThreadDataTask(std::atomic<int>* atomicCounter, std::mutex* mutex, CpuTask::Status* status);
            
            ~InitPerThreadDataTask() override;
            
            MemoryAlloc run() override;
            
        private:
            
            std::mutex*	 IdFactorygetIDMutex;
            std::atomic<int>* _atomicCounter;
        };
        
        
        // thread storage initialization
        SOFA_SIMULATION_CORE_API void initThreadLocalData();
        
        
        
    } // namespace simulation

} // namespace sofa

#endif // InitTasks_h__
