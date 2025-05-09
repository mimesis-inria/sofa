cmake_minimum_required(VERSION 3.22)
project(Sofa.Component.LinearSolver.Iterative LANGUAGES CXX)

set(SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR "src/sofa/component/linearsolver/iterative")

set(HEADER_FILES
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/config.h.in
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/init.h
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/CGLinearSolver.h
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/CGLinearSolver.inl
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/GraphScatteredTypes.h
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/LinearSystemData[GraphScattered].h
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/MatrixFreeSystem[GraphScattered].h
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/MatrixLinearSolver.h
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/MatrixLinearSolver.inl
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/MatrixLinearSystem[GraphScattered].h
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/MinResLinearSolver.h
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/MinResLinearSolver.inl
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/ShewchukPCGLinearSolver.h
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/ShewchukPCGLinearSolver.inl
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/PCGLinearSolver.h
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/PCGLinearSolver.inl
)

set(SOURCE_FILES
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/init.cpp
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/CGLinearSolver.cpp
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/GraphScatteredTypes.cpp
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/LinearSystemData[GraphScattered].cpp
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/MatrixFreeSystem[GraphScattered].cpp
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/MatrixLinearSolver.cpp
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/MatrixLinearSystem[GraphScattered].cpp
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/MinResLinearSolver.cpp
    ${SOFACOMPONENTLINEARSOLVERITERATIVE_SOURCE_DIR}/PCGLinearSolver.cpp
)

sofa_find_package(Sofa.Simulation.Core REQUIRED)
sofa_find_package(Sofa.Component.LinearSystem REQUIRED)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${WRAPPER_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC Sofa.Simulation.Core Sofa.Component.LinearSystem)

sofa_create_package_with_targets(
    PACKAGE_NAME ${PROJECT_NAME}
    PACKAGE_VERSION ${Sofa_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR "${PROJECT_NAME}"
)
