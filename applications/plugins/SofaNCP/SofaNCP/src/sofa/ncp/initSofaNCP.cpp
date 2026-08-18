#include <sofa/ncp/config.h>
#include <sofa/ncp/initSofaNCP.h>

#include <sofa/core/ObjectFactory.h>
using sofa::core::ObjectFactory;
#include <sofa/helper/system/PluginManager.h>

namespace sofa::ncp
{
    extern void registerNCPStaticSolver(sofa::core::ObjectFactory* factory);
    extern void registerNCPDebugNewtonRaphsonSolver(sofa::core::ObjectFactory* factory);
    extern void registerNCPMechanicalContinuationNewtonRaphsonSolver(sofa::core::ObjectFactory* factory);

    extern void registerStaticComplianceProvider(sofa::core::ObjectFactory* factory);

    extern void registerPlaneNCPContactForceField(sofa::core::ObjectFactory* factory);
    extern void registerCylinderNCPContactForceField(sofa::core::ObjectFactory* factory);
    extern void registerSignedDistanceFieldNCPContactForceField(sofa::core::ObjectFactory* factory);
    extern void registerMeshNCPContactForceField(sofa::core::ObjectFactory* factory);

    extern "C" {
        SOFANCP_API void initExternalModule();
        SOFANCP_API const char* getModuleName();
        SOFANCP_API const char* getModuleVersion();
        SOFANCP_API const char* getModuleLicense();
        SOFANCP_API const char* getModuleDescription();
        SOFANCP_API void registerObjects(sofa::core::ObjectFactory* factory);
    }

    void initExternalModule()
    {
        init();
    }

    void init()
    {
        static bool first = true;
        if (first)
        {
            // make sure that this plugin is registered into the PluginManager
            sofa::helper::system::PluginManager::getInstance().registerPlugin(ModuleName);

            first = false;
        }
    }

    const char* getModuleName()
    {
        return sofa::ncp::ModuleName;
    }

    const char* getModuleVersion()
    {
        return sofa::ncp::ModuleVersion;
    }

    const char* getModuleLicense()
    {
        return "LGPL-2.1-or-later";
    }

    const char* getModuleDescription()
    {
        return "Transactional nonlinear solvers and Fischer-Burmeister NCP contact for SOFA.";
    }

    void registerObjects(sofa::core::ObjectFactory* factory)
    {
        registerNCPStaticSolver(factory);
        registerNCPDebugNewtonRaphsonSolver(factory);
        registerNCPMechanicalContinuationNewtonRaphsonSolver(factory);

        registerStaticComplianceProvider(factory);
        
        registerPlaneNCPContactForceField(factory);
        registerCylinderNCPContactForceField(factory);
        registerSignedDistanceFieldNCPContactForceField(factory);
        registerMeshNCPContactForceField(factory);
    }

}
