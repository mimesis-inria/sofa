cmake_minimum_required(VERSION 3.22)

sofa_add_subdirectory(application SceneChecking SceneChecking ON)
sofa_add_subdirectory(application Modeler Modeler OFF)

sofa_add_subdirectory(application getDeprecatedComponents getDeprecatedComponents OFF)

sofa_add_subdirectory(application GenerateRigid GenerateRigid)

sofa_add_subdirectory(application SofaPhysicsAPI SofaPhysicsAPI)
sofa_add_subdirectory(application SofaGuiGlut SofaGuiGlut OFF)

sofa_add_subdirectory(directory SofaGLFW SofaGLFW EXTERNAL GIT_REF master)
sofa_add_subdirectory(plugin Sofa.Qt Sofa.Qt EXTERNAL GIT_REF master OFF)
sofa_add_subdirectory(application runSofa runSofa ON)
sofa_add_subdirectory(application sofaOPENCL sofaOPENCL OFF)

sofa_add_subdirectory(directory Regression Regression EXTERNAL GIT_REF master)
sofa_add_subdirectory(application sofaProjectExample sofaProjectExample)
sofa_add_subdirectory(application sofaInfo sofaInfo)
