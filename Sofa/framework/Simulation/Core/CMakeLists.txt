cmake_minimum_required(VERSION 3.22)
project(Sofa.Simulation.Core)

set(SRC_ROOT "src/sofa/simulation")

set(HEADER_FILES
    ${SRC_ROOT}/config.h.in
    ${SRC_ROOT}/fwd.h
    ${SRC_ROOT}/AnimateBeginEvent.h
    ${SRC_ROOT}/AnimateEndEvent.h
    ${SRC_ROOT}/AnimateVisitor.h
    ${SRC_ROOT}/BaseMechanicalVisitor.h
    ${SRC_ROOT}/BehaviorUpdatePositionVisitor.h
    ${SRC_ROOT}/CleanupVisitor.h
    ${SRC_ROOT}/CollisionAnimationLoop.h
    ${SRC_ROOT}/CollisionBeginEvent.h
    ${SRC_ROOT}/CollisionEndEvent.h
    ${SRC_ROOT}/CollisionVisitor.h
    ${SRC_ROOT}/Colors.h
    ${SRC_ROOT}/CpuTask.h
    ${SRC_ROOT}/CpuTaskStatus.h
    ${SRC_ROOT}/DeactivatedNodeVisitor.h
    ${SRC_ROOT}/DefaultAnimationLoop.h
    ${SRC_ROOT}/DefaultVisualManagerLoop.h
    ${SRC_ROOT}/DeleteVisitor.h
    ${SRC_ROOT}/ExportDotVisitor.h
    ${SRC_ROOT}/ExportGnuplotVisitor.h
    ${SRC_ROOT}/ExportVisualModelOBJVisitor.h
    ${SRC_ROOT}/InitVisitor.h
    ${SRC_ROOT}/IntegrateBeginEvent.h
    ${SRC_ROOT}/IntegrateEndEvent.h
    ${SRC_ROOT}/MechanicalOperations.h
    ${SRC_ROOT}/MechanicalVPrintVisitor.h
    ${SRC_ROOT}/MechanicalVisitor.h
    ${SRC_ROOT}/MutationListener.h
    ${SRC_ROOT}/Node.h
    ${SRC_ROOT}/Node.inl
    ${SRC_ROOT}/ParallelForEach.h
    ${SRC_ROOT}/ParallelSparseMatrixProduct.h
    ${SRC_ROOT}/PauseEvent.h
    ${SRC_ROOT}/PipelineImpl.h
    ${SRC_ROOT}/PositionEvent.h
    ${SRC_ROOT}/PrintVisitor.h
    ${SRC_ROOT}/PropagateEventVisitor.h
    ${SRC_ROOT}/RequiredPlugin.h
    ${SRC_ROOT}/ResetVisitor.h
    ${SRC_ROOT}/SceneLoaderFactory.h
    ${SRC_ROOT}/Simulation.h
    ${SRC_ROOT}/SimulationLoop.h
    ${SRC_ROOT}/SolveVisitor.h
    ${SRC_ROOT}/StateChangeVisitor.h
    ${SRC_ROOT}/TopologyChangeVisitor.h
    ${SRC_ROOT}/UpdateBoundingBoxVisitor.h
    ${SRC_ROOT}/UpdateContextVisitor.h
    ${SRC_ROOT}/UpdateInternalDataVisitor.h
    ${SRC_ROOT}/UpdateLinksVisitor.h
    ${SRC_ROOT}/UpdateMappingEndEvent.h
    ${SRC_ROOT}/UpdateMappingVisitor.h
    ${SRC_ROOT}/VectorOperations.h
    ${SRC_ROOT}/VelocityThresholdVisitor.h
    ${SRC_ROOT}/Visitor.h
    ${SRC_ROOT}/VisitorExecuteFunc.h
    ${SRC_ROOT}/VisualVisitor.h
    ${SRC_ROOT}/WriteStateVisitor.h
    ${SRC_ROOT}/XMLPrintVisitor.h
    ${SRC_ROOT}/init.h
    ${SRC_ROOT}/BaseSimulationExporter.h
    ${SRC_ROOT}/TaskScheduler.h
    ${SRC_ROOT}/TaskSchedulerFactory.h
    ${SRC_ROOT}/TaskSchedulerRegistry.h
    ${SRC_ROOT}/DefaultTaskScheduler.h
    ${SRC_ROOT}/Task.h
    ${SRC_ROOT}/InitTasks.h
    ${SRC_ROOT}/Locks.h
    ${SRC_ROOT}/MainTaskSchedulerFactory.h
    ${SRC_ROOT}/MainTaskSchedulerRegistry.h
    ${SRC_ROOT}/SceneCheck.h
    ${SRC_ROOT}/SceneCheckRegistry.h
    ${SRC_ROOT}/SceneCheckMainRegistry.h
    ${SRC_ROOT}/WorkerThread.h
    ${SRC_ROOT}/events/BuildConstraintSystemEndEvent.h
    ${SRC_ROOT}/events/SimulationInitDoneEvent.h
    ${SRC_ROOT}/events/SimulationInitStartEvent.h
    ${SRC_ROOT}/events/SimulationInitTexturesDoneEvent.h
    ${SRC_ROOT}/events/SimulationStartEvent.h
    ${SRC_ROOT}/events/SimulationStopEvent.h
    ${SRC_ROOT}/events/SolveConstraintSystemEndEvent.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalAccFromFVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalAccumulateJacobian.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalAccumulateJacobian.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalAccumulateMatrixDeriv.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalAddMBK_ToMatrixVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalAddMBKdxVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalAddMDxVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalAddSeparateGravityVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalApplyConstraintsVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalApplyProjectiveConstraint_ToMatrixVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalBeginIntegrationVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalBuildConstraintMatrix.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalComputeContactForceVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalComputeDfVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalComputeEnergyVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalComputeForceVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalComputeGeometricStiffness.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalEndIntegrationVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalGetConstraintInfoVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalGetConstraintJacobianVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalGetMatrixDimensionVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalGetMomentumVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalGetNonDiagonalMassesCountVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalIdentityBlocksInJacobianVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalIntegrateConstraintVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalIntegrationVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalMultiVectorFromBaseVectorVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalMultiVectorPeqBaseVectorVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalMultiVectorToBaseVectorVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalPickParticlesVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalPickParticlesWithTagsVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalProjectJacobianMatrixVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalProjectPositionAndVelocityVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalProjectPositionVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalProjectVelocityVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalPropagateDxAndResetForceVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalPropagateDxVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalPropagateOnlyPositionAndResetForceVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalPropagateOnlyPositionAndVelocityVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalPropagateOnlyPositionVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalPropagateOnlyVelocityVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalResetConstraintVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalResetForceVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalSetPositionAndVelocityVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalVAllocVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalVAvailVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalVDotVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalVFreeVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalVInitVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalVMultiOpVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalVNormVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalVOpVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalVReallocVisitor.h
    ${SRC_ROOT}/mechanicalvisitor/MechanicalVSizeVisitor.h
)

set(SOURCE_FILES
    ${SRC_ROOT}/AnimateBeginEvent.cpp
    ${SRC_ROOT}/AnimateEndEvent.cpp
    ${SRC_ROOT}/AnimateVisitor.cpp
    ${SRC_ROOT}/BaseMechanicalVisitor.cpp
    ${SRC_ROOT}/BehaviorUpdatePositionVisitor.cpp
    ${SRC_ROOT}/CleanupVisitor.cpp
    ${SRC_ROOT}/CollisionAnimationLoop.cpp
    ${SRC_ROOT}/CollisionBeginEvent.cpp
    ${SRC_ROOT}/CollisionEndEvent.cpp
    ${SRC_ROOT}/CollisionVisitor.cpp
    ${SRC_ROOT}/CpuTask.cpp
    ${SRC_ROOT}/CpuTaskStatus.cpp
    ${SRC_ROOT}/DeactivatedNodeVisitor.cpp
    ${SRC_ROOT}/DefaultAnimationLoop.cpp
    ${SRC_ROOT}/DefaultVisualManagerLoop.cpp
    ${SRC_ROOT}/DeleteVisitor.cpp
    ${SRC_ROOT}/ExportDotVisitor.cpp
    ${SRC_ROOT}/ExportGnuplotVisitor.cpp
    ${SRC_ROOT}/ExportVisualModelOBJVisitor.cpp
    ${SRC_ROOT}/InitVisitor.cpp
    ${SRC_ROOT}/IntegrateBeginEvent.cpp
    ${SRC_ROOT}/IntegrateEndEvent.cpp
    ${SRC_ROOT}/MainTaskSchedulerRegistry.cpp
    ${SRC_ROOT}/MainTaskSchedulerFactory.cpp
    ${SRC_ROOT}/MechanicalOperations.cpp
    ${SRC_ROOT}/MechanicalVPrintVisitor.cpp
    ${SRC_ROOT}/MechanicalVisitor.cpp
    ${SRC_ROOT}/MutationListener.cpp
    ${SRC_ROOT}/Node.cpp
    ${SRC_ROOT}/PauseEvent.cpp
    ${SRC_ROOT}/PipelineImpl.cpp
    ${SRC_ROOT}/PositionEvent.cpp
    ${SRC_ROOT}/PrintVisitor.cpp
    ${SRC_ROOT}/PropagateEventVisitor.cpp
    ${SRC_ROOT}/RequiredPlugin.cpp
    ${SRC_ROOT}/ResetVisitor.cpp
    ${SRC_ROOT}/SceneLoaderFactory.cpp
    ${SRC_ROOT}/SceneCheck.cpp
    ${SRC_ROOT}/SceneCheckRegistry.cpp
    ${SRC_ROOT}/SceneCheckMainRegistry.cpp
    ${SRC_ROOT}/Simulation.cpp
    ${SRC_ROOT}/SolveVisitor.cpp
    ${SRC_ROOT}/StateChangeVisitor.cpp
    ${SRC_ROOT}/TopologyChangeVisitor.cpp
    ${SRC_ROOT}/UpdateBoundingBoxVisitor.cpp
    ${SRC_ROOT}/UpdateContextVisitor.cpp
    ${SRC_ROOT}/UpdateInternalDataVisitor.cpp
    ${SRC_ROOT}/UpdateLinksVisitor.cpp
    ${SRC_ROOT}/UpdateMappingEndEvent.cpp
    ${SRC_ROOT}/UpdateMappingVisitor.cpp
    ${SRC_ROOT}/VectorOperations.cpp
    ${SRC_ROOT}/VelocityThresholdVisitor.cpp
    ${SRC_ROOT}/Visitor.cpp
    ${SRC_ROOT}/VisualVisitor.cpp
    ${SRC_ROOT}/WriteStateVisitor.cpp
    ${SRC_ROOT}/XMLPrintVisitor.cpp
    ${SRC_ROOT}/init.cpp
    ${SRC_ROOT}/fwd.cpp
    ${SRC_ROOT}/BaseSimulationExporter.cpp
    ${SRC_ROOT}/TaskScheduler.cpp
    ${SRC_ROOT}/TaskSchedulerFactory.cpp
    ${SRC_ROOT}/TaskSchedulerRegistry.cpp
    ${SRC_ROOT}/DefaultTaskScheduler.cpp
    ${SRC_ROOT}/Task.cpp
    ${SRC_ROOT}/InitTasks.cpp
    ${SRC_ROOT}/WorkerThread.cpp
    ${SRC_ROOT}/events/BuildConstraintSystemEndEvent.cpp
    ${SRC_ROOT}/events/SimulationInitDoneEvent.cpp
    ${SRC_ROOT}/events/SimulationInitStartEvent.cpp
    ${SRC_ROOT}/events/SimulationInitTexturesDoneEvent.cpp
    ${SRC_ROOT}/events/SimulationStartEvent.cpp
    ${SRC_ROOT}/events/SimulationStopEvent.cpp
    ${SRC_ROOT}/events/SolveConstraintSystemEndEvent.cpp

    ${SRC_ROOT}/mechanicalvisitor/MechanicalAccFromFVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalAccumulateJacobian.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalAccumulateJacobian.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalAccumulateMatrixDeriv.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalAddMBK_ToMatrixVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalAddMBKdxVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalAddMDxVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalAddSeparateGravityVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalApplyConstraintsVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalApplyProjectiveConstraint_ToMatrixVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalBeginIntegrationVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalBuildConstraintMatrix.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalComputeContactForceVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalComputeDfVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalComputeEnergyVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalComputeForceVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalComputeGeometricStiffness.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalEndIntegrationVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalGetConstraintInfoVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalGetConstraintJacobianVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalGetMatrixDimensionVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalGetMomentumVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalGetNonDiagonalMassesCountVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalIdentityBlocksInJacobianVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalIntegrateConstraintVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalIntegrationVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalMultiVectorFromBaseVectorVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalMultiVectorPeqBaseVectorVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalMultiVectorToBaseVectorVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalPickParticlesVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalPickParticlesWithTagsVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalProjectJacobianMatrixVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalProjectPositionAndVelocityVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalProjectPositionVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalProjectVelocityVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalPropagateDxAndResetForceVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalPropagateDxVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalPropagateOnlyPositionAndResetForceVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalPropagateOnlyPositionAndVelocityVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalPropagateOnlyPositionVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalPropagateOnlyVelocityVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalResetConstraintVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalResetForceVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalSetPositionAndVelocityVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalVAllocVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalVAvailVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalVDotVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalVFreeVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalVInitVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalVMultiOpVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalVNormVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalVOpVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalVReallocVisitor.cpp
    ${SRC_ROOT}/mechanicalvisitor/MechanicalVSizeVisitor.cpp
)

sofa_find_package(Sofa.Core REQUIRED)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})

target_link_libraries(${PROJECT_NAME} PUBLIC Sofa.Core)

# is this a compiler/linker version specific problem?
if(CMAKE_SYSTEM_NAME STREQUAL Linux)
    # target_link_options is  available from CMake 3.13
    # target_link_options(${PROJECT_NAME} PUBLIC "-pthread")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread" )
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER Sofa.Framework) # IDE folder

sofa_create_package_with_targets(
    PACKAGE_NAME ${PROJECT_NAME}
    PACKAGE_VERSION ${Sofa_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR "${PROJECT_NAME}"
)

# Tests
# If SOFA_BUILD_TESTS exists and is OFF, then these tests will be auto-disabled
cmake_dependent_option(SOFA_SIMULATION_CORE_BUILD_TESTS "Compile the automatic tests" ON "SOFA_BUILD_TESTS OR NOT DEFINED SOFA_BUILD_TESTS" OFF)
if(SOFA_SIMULATION_CORE_BUILD_TESTS)
    add_subdirectory(test)
    add_subdirectory(simutest)
endif()
