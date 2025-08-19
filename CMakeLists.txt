cmake_minimum_required(VERSION 3.16)

# Declare the project and the languages it uses.
project(zin LANGUAGES CXX)

######################################################################
# Compiler settings
######################################################################

# Enforce C++17 across all targets.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(ZIN_BUILD_PYTHON_BINDINGS)
  # Shared objects such as Python modules need position independent code.
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

######################################################################
# Options
######################################################################

option(ZIN_TEST_PYTHON_SCRIPTS                    "" OFF)
option(ZIN_EXPERIMENT_SCALE_MASSIVE               "" OFF)
option(ZIN_EXPERIMENT_SCALE_MINIMUM               "" OFF)
option(ZIN_EXPERIMENT_USE_LAST_SELECTION_STRATEGY "" ON )
option(ZIN_BUILD_APP                              "" ON )
option(ZIN_BUILD_TESTS                            "" ON )
option(ZIN_BUILD_PYTHON_BINDINGS                  "" OFF)

# This option enables the "satisfied" button and sets the maximum number of iterations to 15.
option(ZIN_APP_FOR_STUDY "Custom the app behavior for user study" OFF)

if(ZIN_EXPERIMENT_SCALE_MASSIVE AND ZIN_EXPERIMENT_SCALE_MINIMUM)
  message(FATAL_ERROR "ZIN_EXPERIMENT_SCALE_MASSIVE and ZIN_EXPERIMENT_SCALE_MINIMUM cannot be ON at the same time.")
endif()

######################################################################
# Prerequisites
######################################################################

# Eigen3 provides linear algebra utilities.
find_package(Eigen3 REQUIRED)
if((NOT TARGET Eigen3::Eigen) AND (DEFINED EIGEN3_INCLUDE_DIR))
  add_library(AliasEigen3 INTERFACE)
  target_include_directories(AliasEigen3 INTERFACE ${EIGEN3_INCLUDE_DIR})
  add_library(Eigen3::Eigen ALIAS AliasEigen3)
endif()

if(ZIN_BUILD_APP)
  if(APPLE AND EXISTS /usr/local/opt/qt)
    list(APPEND CMAKE_PREFIX_PATH "/usr/local/opt/qt")
  endif()
  if(APPLE AND EXISTS /usr/local/opt/qt5)
    list(APPEND CMAKE_PREFIX_PATH "/usr/local/opt/qt5")
  endif()
  if(APPLE AND EXISTS /opt/homebrew/opt/qt)
    list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/opt/qt")
  endif()
  if(APPLE AND EXISTS /opt/homebrew/opt/qt5)
    list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/opt/qt5")
  endif()
  if(APPLE AND EXISTS /opt/homebrew/opt/qt@5)
    list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/opt/qt@5")
  endif()
  find_package(Qt5 COMPONENTS Concurrent Widgets REQUIRED)
endif()

######################################################################
# Third parties
######################################################################

set(SEQUENTIAL_LINE_SEARCH_BUILD_COMMAND_DEMOS                 OFF CACHE INTERNAL "")
set(SEQUENTIAL_LINE_SEARCH_BUILD_PHOTO_DEMOS                   OFF CACHE INTERNAL "")
set(SEQUENTIAL_LINE_SEARCH_BUILD_PYTHON_BINDING                OFF CACHE INTERNAL "")
set(SEQUENTIAL_LINE_SEARCH_BUILD_VISUAL_DEMOS                  OFF CACHE INTERNAL "")
set(SEQUENTIAL_LINE_SEARCH_USE_NOISELESS_FORMULATION           ON  CACHE INTERNAL "")
set(SEQUENTIAL_LINE_SEARCH_USE_PARALLELIZED_MULTI_START_SEARCH ON  CACHE INTERNAL "")
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/sequential-line-search)

if(ZIN_BUILD_APP)
  set(ENHANCER_BUILD_QT_TESTS          OFF CACHE INTERNAL "")
  set(ENHANCER_USE_QT_FEATURES         ON  CACHE INTERNAL "")
  set(ENHANCER_USE_ADVANCED_PARAMETERS ON  CACHE BOOL     "")
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/enhancer)

  set(VISOPTSLIDER_BUILD_TESTS OFF CACHE INTERNAL "")
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/visoptslider)

  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/string-util)
endif()

if(ZIN_BUILD_TESTS)
  set(OTF_WITH_EIGEN ON CACHE INTERNAL "")
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/optimization-test-functions)
endif()

######################################################################
# Core library
######################################################################

# Collect library sources. CONFIGURE_DEPENDS lets CMake track new files.
file(GLOB CONFIGURE_DEPENDS lib_headers ${CMAKE_CURRENT_SOURCE_DIR}/sps/include/sps/*.hpp)
file(GLOB CONFIGURE_DEPENDS lib_sources  ${CMAKE_CURRENT_SOURCE_DIR}/sps/src/*.cpp)

add_library(sps STATIC ${lib_headers} ${lib_sources})
target_include_directories(sps PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/sps/include)
target_link_libraries(sps PUBLIC Eigen3::Eigen SequentialLineSearch)

######################################################################
# Main app
######################################################################

if(ZIN_BUILD_APP)
  # Application sources.
  file(
    GLOB CONFIGURE_DEPENDS
    app_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/app/*.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/app/*.cpp
  )

  add_executable(sequential-gallery ${app_sources})
  target_link_libraries(sequential-gallery Eigen3::Eigen enhancer Qt5::Concurrent Qt5::Widgets sps string-util visoptslider)
  add_custom_command(TARGET sequential-gallery POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/resources $<TARGET_FILE_DIR:sequential-gallery>/resources)
  if(ZIN_APP_FOR_STUDY)
    target_compile_definitions(sequential-gallery PUBLIC ZIN_APP_FOR_STUDY)
  endif()
endif()

######################################################################
# Experiment executables
######################################################################

if(ZIN_BUILD_TESTS)
  set(exp_source ${CMAKE_CURRENT_SOURCE_DIR}/exp/main.cpp)
  add_executable(exp-main ${exp_source})
  target_link_libraries(exp-main Eigen3::Eigen optimization-test-functions sps)

  if(ZIN_EXPERIMENT_SCALE_MASSIVE)
    target_compile_definitions(exp-main PUBLIC ZIN_EXPERIMENT_SCALE_MASSIVE)
  endif()
  if(ZIN_EXPERIMENT_SCALE_MINIMUM)
    target_compile_definitions(exp-main PUBLIC ZIN_EXPERIMENT_SCALE_MINIMUM)
  endif()
  if(ZIN_EXPERIMENT_USE_LAST_SELECTION_STRATEGY)
    target_compile_definitions(exp-main PUBLIC ZIN_EXPERIMENT_USE_LAST_SELECTION_STRATEGY)
  endif()
endif()

######################################################################
# Python bindings
######################################################################

if(ZIN_BUILD_PYTHON_BINDINGS)
	# pybind11
	set(PYBIND11_INSTALL OFF CACHE INTERNAL "" FORCE)
	set(PYBIND11_TEST    OFF CACHE INTERNAL "" FORCE)
	add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/sequential-line-search/external/pybind11)

	add_subdirectory(python-bindings)
endif()

######################################################################
# Tests
######################################################################

if(ZIN_BUILD_TESTS)
  enable_testing()

  ##############################################################################
  # Unit tests
  ##############################################################################

  set(BUILD_GMOCK   OFF CACHE INTERNAL "" FORCE)
  set(INSTALL_GTEST OFF CACHE INTERNAL "" FORCE)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/googletest)

  add_executable(unit-test ${CMAKE_CURRENT_SOURCE_DIR}/test/unit-test.cpp)
  target_link_libraries(unit-test sps gtest_main)

  add_test(NAME unit-test COMMAND $<TARGET_FILE:unit-test>)

  ##############################################################################
  # Other tests
  ##############################################################################

  add_test(NAME exp-main COMMAND $<TARGET_FILE:exp-main>)
  if(ZIN_TEST_PYTHON_SCRIPTS)
    add_test(NAME exp-main-python COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/exp/visualize.py $<TARGET_FILE_DIR:exp-main>)
    add_test(NAME exp-main-python-additional COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/exp/visualize-together.py $<TARGET_FILE_DIR:exp-main>)
  endif()
endif()
