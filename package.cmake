######################
# Wrapper macro to set boolean value to a variable
macro(setSofaOption name value)
    set(${name} "${value}" CACHE BOOL "" FORCE)
    message("${name} ${value}")
endmacro()

macro(setSofaPath name value)
    set(${name} "${value}" CACHE PATH "" FORCE)
    message("${name} ${value}")
endmacro()

macro(setSofaString name value)
    set(${name} "${value}" CACHE STRING "" FORCE)
    message("${name} ${value}")
endmacro()

macro(setSofaFilePath name value)
    set(${name} "${value}" CACHE FILEPATH "" FORCE)
    message("${name} ${value}")
endmacro()
######################

setSofaString(CMAKE_BUILD_TYPE Release)

setSofaOption(APPLICATION_RUNSOFA ON)
setSofaOption(APPLICATION_MODELER OFF)

setSofaOption(SOFA_USE_MASK OFF)

setSofaOption(SOFA_BUILD_TESTS OFF)
setSofaOption(SOFA_BUILD_TUTORIALS OFF)

# Set all plugins/modules OFF
get_cmake_property(_variableNames VARIABLES)
list (SORT _variableNames)
foreach (_variableName ${_variableNames})
    if(_variableName MATCHES "^PLUGIN_" OR _variableName MATCHES "^MODULE_")
        setSofaOption(${_variableName} OFF)
    endif()
endforeach()

# Set some plugins/modules ON
setSofaOption(PLUGIN_SOFAALLCOMMONCOMPONENTS ON)
setSofaOption(PLUGIN_CIMGPLUGIN ON)
setSofaOption(PLUGIN_SOFAPYTHON ON)
setSofaOption(PLUGIN_SOFAMISCCOLLISION ON)
setSofaOption(MODULE_SOFASPARSESOLVER ON)
setSofaOption(MODULE_SOFAPRECONDITIONER ON)

# Copy resources files (etc/, share/, examples/) when installing 
setSofaOption(SOFA_INSTALL_RESOURCES_FILES ON)

# install GTest even if SOFA_BUILD_TESTS=OFF
add_subdirectory(extlibs/gtest)
