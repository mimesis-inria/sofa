cmake_minimum_required(VERSION 3.12)
project(SofaNCP VERSION 0.1.0 LANGUAGES CXX)

if(NOT TARGET Sofa.Framework)
    find_package(Sofa.Framework REQUIRED)
endif()

sofa_find_package(Sofa.Simulation.Core REQUIRED)
sofa_find_package(Sofa.Component.LinearSystem REQUIRED)
sofa_find_package(Sofa.Component.ODESolver.Backward REQUIRED)
sofa_find_package(Sofa.Component.Constraint.Projective REQUIRED)
sofa_find_package(Sofa.Component.SolidMechanics.FEM.Elastic REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)

set(SOFANCP_SRC_DIR src/sofa/ncp)

set(HEADER_FILES
    ${SOFANCP_SRC_DIR}/config.h.in
    ${SOFANCP_SRC_DIR}/initSofaNCP.h

    ${SOFANCP_SRC_DIR}/solvers/NCPStaticSolver.h
    ${SOFANCP_SRC_DIR}/solvers/NCPDebugNewtonRaphsonSolver.h
    ${SOFANCP_SRC_DIR}/solvers/NCPMechanicalContinuationNewtonRaphsonSolver.h

    ${SOFANCP_SRC_DIR}/contact/FischerBurmeisterContactForceField.h
    ${SOFANCP_SRC_DIR}/contact/FischerBurmeisterContactForceField.inl
    ${SOFANCP_SRC_DIR}/contact/PlaneNCPContactForceField.h
    ${SOFANCP_SRC_DIR}/contact/PlaneNCPContactForceField.inl
    ${SOFANCP_SRC_DIR}/contact/CylinderNCPContactForceField.h
    ${SOFANCP_SRC_DIR}/contact/CylinderNCPContactForceField.inl
    ${SOFANCP_SRC_DIR}/contact/SignedDistanceFieldNCPContactForceField.h
    ${SOFANCP_SRC_DIR}/contact/SignedDistanceFieldNCPContactForceField.inl
    ${SOFANCP_SRC_DIR}/contact/MeshNCPContactForceField.h
    ${SOFANCP_SRC_DIR}/contact/MeshNCPContactForceField.inl

    ${SOFANCP_SRC_DIR}/preconditioners/StaticComplianceProvider.h
    ${SOFANCP_SRC_DIR}/preconditioners/StaticComplianceProvider.inl
)
set(HEADER_FILES_TO_MOC
)

set(SOURCE_FILES
    ${SOFANCP_SRC_DIR}/initSofaNCP.cpp

    ${SOFANCP_SRC_DIR}/solvers/NCPStaticSolver.cpp
    ${SOFANCP_SRC_DIR}/solvers/NCPDebugNewtonRaphsonSolver.cpp
    ${SOFANCP_SRC_DIR}/solvers/NCPMechanicalContinuationNewtonRaphsonSolver.cpp

    ${SOFANCP_SRC_DIR}/contact/PlaneNCPContactForceField.cpp
    ${SOFANCP_SRC_DIR}/contact/CylinderNCPContactForceField.cpp
    ${SOFANCP_SRC_DIR}/contact/SignedDistanceFieldNCPContactForceField.cpp
    ${SOFANCP_SRC_DIR}/contact/MeshNCPContactForceField.cpp
    
    ${SOFANCP_SRC_DIR}/preconditioners/StaticComplianceProvider.cpp
)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES_TO_MOC} ${MOCCED_HEADER_FILES} ${HEADER_FILES} ${SOURCE_FILES} ${README_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC Sofa.Framework Sofa.Simulation.Core Sofa.Component.LinearSystem Sofa.Component.ODESolver.Backward Sofa.Component.Constraint.Projective Sofa.Component.SolidMechanics.FEM.Elastic Eigen3::Eigen)

sofa_create_package_with_targets(
    PACKAGE_NAME ${PROJECT_NAME}
    PACKAGE_VERSION ${PROJECT_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR ${PROJECT_NAME}
    RELOCATABLE "plugins"
)
