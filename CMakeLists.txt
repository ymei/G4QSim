#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(G4QSim)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT REQUIRED)
if(NOT ROOT_FOUND)
  message(STATUS "ROOT package not found.")
  return()
endif()
if(NOT GCCXML)
   message(STATUS "GCCXML not found.")
#  return()
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIRS})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Generate dictionaries, add ROOT libraries properties
#
# REFLEX_GENERATE_DICTIONARY(ExP01Classes include/ExP01Classes.hh SELECTION xml/selection.xml)
# add_library(ExP01ClassesDict SHARED ${sources} ExP01Classes.cpp)
# set(libsuffix .so)
# set(ROOT_LIBRARY_PROPERTIES ${ROOT_LIBRARY_PROPERTIES} SUFFIX ${libsuffix})
# set_target_properties(ExP01ClassesDict PROPERTIES ${ROOT_LIBRARY_PROPERTIES})
# target_link_libraries(ExP01ClassesDict  ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} Cintex Reflex)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(G4QSim G4QSim.cc ${sources} ${headers})
target_link_libraries(G4QSim ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build G04. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(G4QSim_SCRIPTS
    )

foreach(_script ${G4QSim_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
# add_custom_target(G4QSim DEPENDS G4QSim)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS G4QSim DESTINATION bin)
