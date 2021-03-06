cmake_minimum_required (VERSION 3.14)

###############################################################################
#  Get the Current Version Number so that it can be embedded in files later
###############################################################################
file(STRINGS "${CMAKE_CURRENT_LIST_DIR}/VERSION" _VERSION_DATA REGEX ^[0-9]+_[0-9]+_[0-9]+)
string(REPLACE "_" ";" _VERSION_LIST ${_VERSION_DATA})
list(GET _VERSION_LIST 0 MAJOR_VERSION)
list(GET _VERSION_LIST 1 MINOR_VERSION)
list(GET _VERSION_LIST 2 MICRO_VERSION)

set(CMAKE_PROJECT_VERSION ${MAJOR_VERSION}.${MINOR_VERSION}.${MICRO_VERSION})


###############################################################################
##  From:
##    https://cliutils.gitlab.io/modern-cmake/chapters/projects/submodule.html
##
##  Use git to add, initialize, and update submodules for the user
###############################################################################
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    # Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "** Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

################################################################################
# TODO:  Add licensing and authorship information
# TODO:  Test with FEtk build
# TODO:  Handle special mac dependencies
#        (configure.ac:1306)
################################################################################
set(CMAKE_VERBOSE_MAKEFILE true)
cmake_policy(SET CMP0054 OLD)
cmake_policy(SET CMP0042 NEW)
set(CMAKE_MACOSX_RPATH 1)

################################################################################
# Set up basic project stuff                                                   #
################################################################################
include(CheckIncludeFiles)
include(CheckFunctionExists)
include(ExternalProject)
include(FeatureSummary)
include(InstallRequiredSystemLibraries)

set(APBS_VERSION     "${MAJOR_VERSION}.${MINOR_VERSION}.${MICRO_VERSION}")
set(PACKAGE_STRING   "APBS ${APBS_VERSION}")
set(CMAKE_BUILD_TYPE "RELWITHDEBINFO")

################################################################################
#  NOTE: Overwrite the RPATH information during the install from
#  https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling
################################################################################
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
#include(GNUInstallDirs)
#file(RELATIVE_PATH relDir
#     ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
#     ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}
#)
#set(CMAKE_INSTALL_RPATH $ORIGIN $ORIGIN/${relDir})

project("APBS" VERSION ${APBS_VERSION})
set(CMAKE_PROJECT_DESCRIPTION "Adaptive Poisson-Boltzmann Solver")
set(CMAKE_PROJECT_HOMEPAGE_URL "http://www.poissonboltzmann.org/")

set_package_properties(APBS PROPERTIES
    URL "${CMAKE_PROJECT_HOMEPAGE_URL}"
    DESCRIPTION "${CMAKE_PROJECT_DESCRIPTION}"
    PURPOSE "APBS solves the equations of continuum electrostatics for large \
             biomolecular assemblages. This software was designed 'from the \
             ground up' using modern design principles to ensure its ability \
             to interface with other computational packages and evolve as \
             methods and applications change over time. The APBS code is \
             accompanied by extensive documentation for both users and \
             programmers and is supported by a variety of utilities for \
             preparing calculations and analyzing results. Finally, the free, \
             open-source APBS license ensures its accessibility to the entire \
             biomedical community.")

###############################################################################
##  NOTE: For Debugging and generally nice information to have
###############################################################################
string(TIMESTAMP COMPILE_TIME)
string(TIMESTAMP COMPILE_YEAR "%Y")

message(STATUS "#############################################################")
message(STATUS "** PROJECT NAME:  ${CMAKE_PROJECT_NAME}")
message(STATUS "** PROJECT DESC:  ${CMAKE_PROJECT_DESCRIPTION}")
message(STATUS "** PROJECT URL:   ${CMAKE_PROJECT_HOMEPAGE_URL}")
message(STATUS "** VERSION_DATA:  ${CMAKE_PROJECT_VERSION}")
message(STATUS "** CURRENT_TIME:  ${COMPILE_TIME}")
message(STATUS "** CURRENT_YEAR:  ${COMPILE_YEAR}")

OPTION(BUILD_SHARED_LIBS "Build shared libraries." OFF)

################################################################################
# Set project paths                                                            #
################################################################################
message(STATUS "Setting project paths")

set(CMAKE_CXX_STANDARD 17)
if(WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -W4")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:100000000")
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -fpermissive -fPIC")
endif()
if(APPLE AND BUILD_SHARED_LIBS AND (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU"))
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -undefined dynamic_lookup")
endif()

set(APBS_ROOT ${PROJECT_SOURCE_DIR})
set(APBS_BUILD ${PROJECT_BINARY_DIR})
set(EXECUTABLE_OUTPUT_PATH ${APBS_BUILD}/bin)
set(LIBRARY_OUTPUT_PATH ${APBS_BUILD}/lib)
set(TOOLS_PATH ${APBS_BUILD}/tools)
set(APBS_BINARY ${EXECUTABLE_OUTPUT_PATH}/apbs)

set(LIBRARY_INSTALL_PATH lib)
set(HEADER_INSTALL_PATH include/apbs)
set(EXECUTABLE_INSTALL_PATH bin)
set(SHARE_INSTALL_PATH share/apbs)

set(APBS_LIBS)
set(APBS_LIB_DIRS)

find_file(CONTRIB_PATH "contrib"
          PATHS "${APBS_ROOT}"
          DOC "The path to contributed modules for apbs")
find_file(EXTERNALS_PATH "externals"
          PATHS "${APBS_ROOT}"
          DOC "The path to the external git submodules")

set(SYS_LIBPATHS /usr/lib64)
#list(APPEND SYS_LIBPATHS /usr/lib64)


###############################################################################
#####  Boilerplate Setup
###############################################################################
# include(GNUInstallDirs)
# set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
# message(STATUS "** Project will be installed to ${CMAKE_INSTALL_PREFIX}")


###############################################################################
#####  Offer the builder the choice of overriding the installation directories
###############################################################################
# set(INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR} CACHE PATH "Installation directory for libraries")
# set(INSTALL_BINDIR ${CMAKE_INSTALL_BINDIR} CACHE PATH "Installation directory for executables")
# set(INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR} CACHE PATH "Installation directory for header files")
# set(INSTALL_CMAKEDIR ${DEF_INSTALL_CMAKEDIR} CACHE PATH "Installation directory for CMake files")


###############################################################################
##### Report to builder
###############################################################################
# foreach(p LIB BIN INCLUDE CMAKE)
#   file(TO_NATIVE_PATH ${CMAKE_INSTALL_PREFIX}/${INSTALL_${p}DIR} _path )
#   message(STATUS "** Installing ${p} components to ${_path}")
#   unset(_path)
# endforeach()


################################################################################
# Debian/Ubuntu hack                                                           #
################################################################################
if(EXISTS "/etc/debian_version" AND
  "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed")
endif()


################################################################################
# Set up temporary files and directories                                       #
################################################################################
file(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}/temp)


################################################################################
# Set the lookup paths for external library dependencies                       #
################################################################################
message(STATUS "Setting lookup paths for headers and libraries")

set(CMAKE_INCLUDE_PATH "${CMAKE_INCLUDE_PATH}")
list(APPEND CMAKE_INCLUDE_PATH /usr/include)
list(APPEND CMAKE_INCLUDE_PATH /usr/local/include)


################################################################################
# Enable ansi pedantic compiling                                               #
################################################################################
option(ENABLE_PEDANTIC "Enable the pedantic ansi compilation" OFF)

if(ENABLE_PEDANTIC)
    ADD_DEFINITIONS("-Wall -pedantic -ansi")
    message(STATUS "Pedantic compilation enabled")
endif()


################################################################################
# Determine Machine Epsilon values                                             #
################################################################################
OPTION(CHECK_EPSILON "Compute machine epsilon values" YES)

if(CHECK_EPSILON)
    message(STATUS "Computing machine epsilon values")
    try_run(
        FLOAT_EPSILON_COMPILED
        FLOAT_EPSILON_COMPUTED
        ${EXECUTABLE_OUTPUT_PATH}/temp
        ${APBS_ROOT}/src/.config/float_epsilon.c
        RUN_OUTPUT_VARIABLE FLOAT_EPSILON_OUTPUT
    )

    if(FLOAT_EPSILON_COMPUTED)
        message(STATUS "Floating point epsilon is ${FLOAT_EPSILON_OUTPUT}")
        set(FLOAT_EPSILON ${FLOAT_EPSILON_OUTPUT})
    else()
        message(FATAL_ERROR "Couldn't compute floating point machine epsilon")
    endif()

    try_run(
        DOUBLE_EPSILON_COMPILED
        DOUBLE_EPSILON_COMPUTED
        ${EXECUTABLE_OUTPUT_PATH}/temp
        ${APBS_ROOT}/src/.config/double_epsilon.c
        RUN_OUTPUT_VARIABLE DOUBLE_EPSILON_OUTPUT
    )

    if(DOUBLE_EPSILON_COMPUTED)
        message(STATUS "Double precision epsilon is ${DOUBLE_EPSILON_OUTPUT}")
        set(DOUBLE_EPSILON ${DOUBLE_EPSILON_OUTPUT})
    else()
        message(FATAL_ERROR "Couldn't compute double precision machine epsilon")
    endif()
else()
  set(FLOAT_EPSILON "1.19209290e-7")
  set(DOUBLE_EPSILON "2.2204460492503131e-16")
endif()



################################################################################
# Check for a few required functions                                           #
################################################################################
CHECK_FUNCTION_EXISTS(time HAVE_TIME_FUNC)
if(NOT HAVE_TIME_FUNC)
    message(FATAL_ERROR "Required time function not found")
endif()


CHECK_FUNCTION_EXISTS(rand HAVE_RAND_FUNC)
if(NOT HAVE_RAND_FUNC)
    message(FATAL_ERROR "Required rand function not found")
endif()


CHECK_FUNCTION_EXISTS(srand HAVE_SRAND_FUNC)
if(NOT HAVE_SRAND_FUNC)
    message(FATAL_ERROR "Required srand function not found")
endif()


################################################################################
# Find some libraries; Windows finds these automatically                       #
################################################################################
if(NOT WIN32)
    find_library(MATH_LIBRARY "m")
    list(APPEND APBS_LIBS m stdc++)
endif()


################################################################################
# Create APBS configuration file
################################################################################
configure_file(
    src/.config/apbscfg.h.in
    ${APBS_BUILD}/src/apbscfg.h
    @ONLY
)
#file(READ ${APBS_BUILD}/src/apbscfg.h APBS_CFG)
#message(STATUS "** APBS_CFG")
#message(FATAL_ERROR ${APBS_CFG})


################################################################################
# Optionally copy NanoShaper executable                                        #
# the actual grab is later for BEM enable                                      #
################################################################################
option(GET_NanoShaper "Put NanoShaper executable in the build directory" OFF)


################################################################################
# Optionally build BEM: TABI-PB method                                         #
################################################################################
option(ENABLE_BEM "Boundary element method using TABIPB" OFF)

if(ENABLE_BEM)
    message(STATUS "Building TABIPB")
    add_definitions(-DENABLE_BEM)
    set(TABI bem)
    set(TABI_PATH ${EXTERNALS_PATH}/${TABI})

    set(GET_NanoShaper ON) #CACHE BOOL "Grabbing NanoShaper path"

    if(WIN32)
        ExternalProject_add(${TABI}
            SOURCE_DIR ${TABI_PATH}
            CMAKE_ARGS -DENABLE_TABIPB_APBS=ON -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                       -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            BUILD_COMMAND "${CMAKE_COMMAND}" --build . --config $<CONFIG> --target TABIPBlib
            INSTALL_COMMAND "")
    else()
        ExternalProject_add(${TABI}
            SOURCE_DIR ${TABI_PATH}
            CMAKE_ARGS -DENABLE_TABIPB_APBS=ON -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                       -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}
                       -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            BUILD_COMMAND "${CMAKE_COMMAND}" --build . --target TABIPBlib
            INSTALL_COMMAND "")
    endif()

    ExternalProject_Get_Property(${TABI} binary_dir)
    list(APPEND APBS_LIB_DIRS "${binary_dir}/lib")
    list(APPEND APBS_LIBS TABIPBlib)
    include_directories(${EXTERNALS_PATH}/${TABI})
    include_directories(${EXTERNALS_PATH}/${TABI}/src)
    include_directories(${EXTERNALS_PATH}/${TABI}/src/tabipb_wrap)

    if(BUILD_SHARED_LIBS)
        install(FILES ${binary_dir}/lib/libTABIPBlib${CMAKE_SHARED_LIBRARY_SUFFIX}
                DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
    endif()
endif() # ENABLE_BEM


################################################################################
# Getting nanoshaper binary
################################################################################
set(MESH_ROOT_URL_PATH "https://github.com/lwwilson1/mesh_routines/releases/download/v1.6")

if(GET_NanoShaper)
    set(NANOSHAPER_FULL_URL_PATH)
    set(NANOSHAPER_RAW_BINARY_NAME "NanoShaper")
    message(STATUS "Copying NanoShaper executable")
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")      # Mac OS X specific code
        set(NANOSHAPER_FULL_URL_PATH "${MESH_ROOT_URL_PATH}/NanoShaper_OSX")
    elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")   # Linux specific code
        set(NANOSHAPER_FULL_URL_PATH "${MESH_ROOT_URL_PATH}/NanoShaper_Linux64")
    elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows") # Windows specific code
        set(NANOSHAPER_FULL_URL_PATH "${MESH_ROOT_URL_PATH}/NanoShaper32.exe")
        set(NANOSHAPER_RAW_BINARY_NAME "NanoShaper.exe")
    endif()

    if(NOT EXISTS "${EXECUTABLE_OUTPUT_PATH}/${NANOSHAPER_RAW_BINARY_NAME}")
        file(DOWNLOAD ${NANOSHAPER_FULL_URL_PATH} ${EXECUTABLE_OUTPUT_PATH}/${NANOSHAPER_RAW_BINARY_NAME})
        message(STATUS "NanoShaper Download: ${NANOSHAPER_FULL_URL_PATH} ${EXECUTABLE_OUTPUT_PATH}/${NANOSHAPER_RAW_BINARY_NAME}")
    endif()

    install(PROGRAMS ${EXECUTABLE_OUTPUT_PATH}/${NANOSHAPER_RAW_BINARY_NAME}
            DESTINATION ${EXECUTABLE_INSTALL_PATH})
    message(STATUS "NanoShaper install: ${EXECUTABLE_OUTPUT_PATH}/${NANOSHAPER_RAW_BINARY_NAME} ${EXECUTABLE_INSTALL_PATH}")
endif()


################################################################################
# Handle conditional building with MPI Support                                 #
# There may be a better way, but for now this comes before FETK because FETK   #
# depends on having the libraries located.                                     #
################################################################################
option(ENABLE_MPI "Enable MPI parallelism" OFF)

if(ENABLE_MPI)
    if(NOT ENABLE_DEBUG)
        message(STATUS "Checking for MPI")
        find_package(MPI)
        if(MPI_FOUND)
            message(STATUS "MPI support enabled")
            set(HAVE_MPI_H 1)
            include_directories(${MPI_C_INCLUDE_PATH})
            list(APPEND APBS_LIBS ${MPI_C_LIBRARIES})
        else()
            message(WARNING "MPI was not found; disabling")
        endif()
    else()
        message(WARNING "MPI may not be enabled in debugging mode")
    endif()
endif()


################################################################################
# Handle the finite element solver dependencies                                #
################################################################################
option(ENABLE_FETK "Enable the finite element solver" OFF)


#  HACK START - This used to be in tools/windows/configuration headers/*.h
#file(COPY tools/apbscfg.h DESTINATION ${APBS_BUILD}/src)
file(COPY tools/maloccf.h DESTINATION ${APBS_BUILD}/src)
#  if should be done something like the following:
#configure_file(
#     ${EXTERNALS_PATH}/fetk/maloc/src/config/maloccf.h.in
#     DESTINATION ${APBS_BUILD}/src/maloccf.h
#)
#  HACK END   - This used to be in tools/windows/configuration headers/*.h


################################################################################
# Deal with MALOC                                                              #
# Note that we depend exclusively on the MALOC that is in our external FETK    #
# repository.  Also note, that in order for this to work on Windows, we've had #
# to add CMake support to MALOC, and that FETK has no such support, in general.#
################################################################################
# add_library(apbs_fetk_libraries INTERFACE)

if(ENABLE_FETK)
    if(WIN32)
        message(FATAL_ERROR "FETK will not built on Windows because it uses autotools!")
    else()
        message(STATUS "Building FETK components.")
        set(CONFIG_OPTS --prefix=<INSTALL_DIR> --without-doxygen --with-doxygen=no)

        set(MPI_ENV "")
        if(HAVE_MPI_H)
            list(APPEND CONFIG_OPTS --enable-mpi)

            # MALOC just cares about libmpi
            foreach(LIB_PATH ${MPI_C_LIBRARIES})
                get_filename_component(LIB ${LIB_PATH} NAME_WE)
                if(LIB STREQUAL "libmpi")
                    get_filename_component(MPI_LIB_DIR ${LIB_PATH} DIRECTORY)
                endif()
            endforeach()

            set(MPI_ENV "FETK_MPI_INCLUDE=${MPI_C_INCLUDE_PATH}" "FETK_MPI_LIBRARY=${MPI_LIB_DIR}")
        endif()

        ExternalProject_Add(fetk-maloc
            PREFIX fetk
            SOURCE_DIR ${EXTERNALS_PATH}/fetk/maloc/src
            CONFIGURE_COMMAND ${MPI_ENV} ${EXTERNALS_PATH}/fetk/maloc/configure ${CONFIG_OPTS}
            BUILD_COMMAND ${MAKE})

        ExternalProject_Add(fetk-punc
            PREFIX fetk
            DEPENDS fetk-maloc
            SOURCE_DIR ${EXTERNALS_PATH}/fetk/punc/src
            CONFIGURE_COMMAND ${MPI_ENV} FETK_INCLUDE=<INSTALL_DIR>/include FETK_LIBRARY=<INSTALL_DIR>/lib ${EXTERNALS_PATH}/fetk/punc/configure ${CONFIG_OPTS}
            BUILD_COMMAND ${MAKE})

        ExternalProject_Add(fetk-gamer
            PREFIX fetk
            DEPENDS fetk-maloc
            SOURCE_DIR ${EXTERNALS_PATH}/fetk/gamer/src
            CONFIGURE_COMMAND ${MPI_ENV} FETK_INCLUDE=<INSTALL_DIR>/include FETK_LIBRARY=<INSTALL_DIR>/lib ${EXTERNALS_PATH}/fetk/gamer/configure ${CONFIG_OPTS}
            BUILD_COMMAND ${MAKE})

        ExternalProject_Add(fetk-mc
            PREFIX fetk
            DEPENDS fetk-maloc fetk-punc
            SOURCE_DIR ${EXTERNALS_PATH}/fetk/mc/src
            CONFIGURE_COMMAND ${MPI_ENV} FETK_INCLUDE=<INSTALL_DIR>/include FETK_LIBRARY=<INSTALL_DIR>/lib ${EXTERNALS_PATH}/fetk/mc/configure ${CONFIG_OPTS}
            BUILD_COMMAND ${MAKE})

        set(FETK_ENABLED 1)

        ExternalProject_Get_Property(fetk-maloc install_dir)
        include_directories(${install_dir}/include)

        list(APPEND APBS_LIB_DIRS "${install_dir}/lib")
        # target_link_libraries(apbs_fetk_libraries INTERFACE maloc amd punc mc gamer superlu umfpack blas vf2c tetgen triangle readline)
        list(APPEND APBS_LIBS
                                  maloc
                                  amd
                                  punc
                                  mc
                                  gamer
                                  superlu
                                  umfpack
                                  blas
                                  vf2c
                                  tetgen
                                  triangle
                                  readline)

        if(BUILD_SHARED_LIBS)
            install(DIRECTORY ${install_dir}/lib/
                    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
        endif()

        SET(HAVE_MC_H 1)
        SET(HAVE_PUNC_H 1)
        SET(HAVE_BIOM_H 1)
    endif()
else()
    # We'll use autotools to build MALOC if FETK is enabled.
    add_subdirectory(${EXTERNALS_PATH}/fetk/maloc)
    include_directories(${EXTERNALS_PATH}/fetk/maloc/src/base
                        ${EXTERNALS_PATH}/fetk/maloc/src/vsys
                        ${EXTERNALS_PATH}/fetk/maloc/src/vsh
                        ${EXTERNALS_PATH}/fetk/maloc/src/psh
                        ${EXTERNALS_PATH}/fetk/mc/src/aaa_inc)
    list(APPEND APBS_LIB_DIRS ${EXTERNALS_PATH}/fetk/maloc/lib)
    list(APPEND APBS_LIBS maloc)
    if (WIN32)
      list(APPEND MALOC_LIBS wsock32 WS2_32)
    endif()
endif()


################################################################################
# Optionally build geoflow (git submodule)
################################################################################
option(ENABLE_GEOFLOW "geoflow method" OFF)

if(ENABLE_GEOFLOW)
    message(STATUS "Building GEOFLOW")
    add_definitions(-DGEOFLOW_APBS)

    if(WIN32)
        ExternalProject_Add(geoflow
            SOURCE_DIR ${EXTERNALS_PATH}/geoflow_c
            CMAKE_ARGS -DENABLE_GEOFLOW_APBS=ON -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
            BUILD_COMMAND "${CMAKE_COMMAND}" --build . --config $<CONFIG> --target GeometricFlowLib
            INSTALL_COMMAND "")
    else()
        ExternalProject_Add(geoflow
            SOURCE_DIR ${EXTERNALS_PATH}/geoflow_c
            CMAKE_ARGS -DENABLE_GEOFLOW_APBS=ON -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
            BUILD_COMMAND "${CMAKE_COMMAND}" --build . --target GeometricFlowLib
            INSTALL_COMMAND "")
    endif()

    ExternalProject_Get_Property(geoflow binary_dir)
    add_definitions(-DENABLE_GEOFLOW)
    list(APPEND APBS_LIB_DIRS "${binary_dir}/lib")
    list(APPEND APBS_LIBS GeometricFlowLib)
    include_directories(${EXTERNALS_PATH}/geoflow_c/src)

    if(BUILD_SHARED_LIBS)
        install(FILES ${binary_dir}/lib/libGeometricFlowLib${CMAKE_SHARED_LIBRARY_SUFFIX}
                DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
    endif()
endif()


################################################################################
# Optionally build pbam (git submodule)
################################################################################
option(ENABLE_PBAM "pbam method" OFF)

if(ENABLE_PBAM)
    message(STATUS "Building PBAM")
    add_definitions(-DPBAM_APBS)

    if(WIN32)
        ExternalProject_Add(pbam
            SOURCE_DIR ${EXTERNALS_PATH}/pb_s_am
            CMAKE_ARGS -DENABLE_PBAM_APBS=ON -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
            BUILD_COMMAND "${CMAKE_COMMAND}" --build . --config $<CONFIG> --target PBAMLib
            INSTALL_COMMAND "")
    else()
        ExternalProject_Add(pbam
            SOURCE_DIR ${EXTERNALS_PATH}/pb_s_am
            CMAKE_ARGS -DENABLE_PBAM_APBS=ON -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
            BUILD_COMMAND "${CMAKE_COMMAND}" --build . --target PBAMLib
            INSTALL_COMMAND "")
    endif()

    ExternalProject_Get_Property(pbam binary_dir)
    add_definitions(-DENABLE_PBAM)
    list(APPEND APBS_LIB_DIRS "${binary_dir}/pbam/src")
    if (WIN32)
        list(APPEND APBS_LIB_DIRS "${binary_dir}/pbam/src/${CMAKE_BUILD_TYPE}")
    endif()
    list(APPEND APBS_LIBS PBAMLib)
    include_directories(${EXTERNALS_PATH}/pb_s_am/pbam/src)
    include_directories(${EXTERNALS_PATH}/pb_s_am/pb_shared/src)
    include_directories(${EXTERNALS_PATH}/pb_s_am/pb_wrap/src)

    if(BUILD_SHARED_LIBS)
        install(FILES ${binary_dir}/pbam/src/libPBAMLib${CMAKE_SHARED_LIBRARY_SUFFIX}
                DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
    endif()
endif()


################################################################################
# Optionally build pbsam (git submodule)
################################################################################
option(ENABLE_PBSAM "pbsam method" OFF)

if(ENABLE_PBSAM)
    message(STATUS "Building PBSAM")
    add_definitions(-DPBSAM_APBS)

    list(APPEND CMAKE_MODULE_PATH ${EXTERNALS_PATH}/pb_s_am/cmake/modules)
    include(${EXTERNALS_PATH}/pb_s_am/cmake/Utils.cmake)
    include(${EXTERNALS_PATH}/pb_s_am/cmake/Dependencies.cmake)
    include(${EXTERNALS_PATH}/pb_s_am/cmake/Summary.cmake)

    if(WIN32)
        ExternalProject_Add(pbsam
            SOURCE_DIR ${EXTERNALS_PATH}/pb_s_am
            CMAKE_ARGS -DENABLE_PBSAM_APBS=ON -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
            BUILD_COMMAND "${CMAKE_COMMAND}" --build . --config $<CONFIG> --target PBSAMLib
            INSTALL_COMMAND "")
    else()
        ExternalProject_Add(pbsam
            SOURCE_DIR ${EXTERNALS_PATH}/pb_s_am
            CMAKE_ARGS -DENABLE_PBSAM_APBS=ON -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
            BUILD_COMMAND "${CMAKE_COMMAND}" --build . --target PBSAMLib
            INSTALL_COMMAND "")
    endif()

    ExternalProject_Get_Property(pbsam binary_dir)
    add_definitions(-DENABLE_PBSAM)
    list(APPEND APBS_LIB_DIRS "${binary_dir}/pbsam/src")

    if(WIN32)
        list(APPEND APBS_LIBS PBSAMLib ${PBSAM_LINKER_LIBS})
    else()
        list(APPEND APBS_LIBS stdc++ PBSAMLib ${PBSAM_LINKER_LIBS})
    endif()

    message(STATUS "Adding in " ${PBSAM_LINKER_LIBS} )
    include_directories(${EXTERNALS_PATH}/pb_s_am/pbsam/src)
    include_directories(${EXTERNALS_PATH}/pb_s_am/pb_shared/src)
    include_directories(${EXTERNALS_PATH}/pb_s_am/pb_wrap/src)

    if(BUILD_SHARED_LIBS)
        install(FILES ${binary_dir}/pbsam/src/libPBSAMLib${CMAKE_SHARED_LIBRARY_SUFFIX}
                DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
    endif()
endif()


################################################################################
# Handle conditional fast mode                                                 #
################################################################################
option(ENABLE_FAST "Enable fast mode" OFF)

if(ENABLE_FAST)
    set(APBS_FAST 1)
    message(STATUS "Fast mode enabled")
endif()


################################################################################
# Handle conditional TINKER support                                            #
################################################################################
option(ENABLE_TINKER "Enable TINKER support" OFF)

if(ENABLE_TINKER)
    set(WITH_TINKER 1)
    message(STATUS "Tinker enabled")
endif()


################################################################################
# Handle conditional availability of macro embedding                           #
################################################################################
try_compile(
    HAVE_EMBED
    ${APBS_BUILD}/build
    ${APBS_ROOT}/src/.config/embed_test.c
)

# TODO: Determine if the EMBED macro is even used


################################################################################
# Handle conditional debug building                                            #
################################################################################
option(ENABLE_DEBUG "Enable debugging mode" OFF)

if(ENABLE_DEBUG)
    set(CMAKE_BUILD_TYPE "Debug")
    set(DEBUG 1)
    message(STATUS "Debugging compilation enabled")
endif()


################################################################################
# Enable inline functions conditionally dependent on debug mode                #
################################################################################
option(ENABLE_INLINE "Enable inline functions" ON)

if(ENABLE_INLINE)
    set(APBS_NOINLINE 1)
    message(STATUS "Inline functions enabled")
endif()


################################################################################
# Handle conditional building with quiet mode                                  #
################################################################################
option(ENABLE_QUIET "Enable quiet mode" OFF)

if(ENABLE_QUIET)
    set(VAPBSQUIET 1)
    message(STATUS "Quiet mode enabled")
endif()


################################################################################
# Handle conditional building with verbose debugging information printouts     #
################################################################################
option(ENABLE_VERBOSE_DEBUG "Enable verbose debugging mode" ON)

if(ENABLE_VERBOSE_DEBUG)
    set(VERBOSE_DEBUG 1)
    message(STATUS "Verbose debugging mode enabled")
endif()


################################################################################
# Configure Python                                                             #
# Buidling the Python goodness happens in the tool/python directory.           #
################################################################################
option(ENABLE_PYTHON "Enable python support" OFF)

#If(ENABLE_PYTHON AND UNIX AND NOT APPLE AND NOT BUILD_SHARED_LIBS)
#    message(FATAL_ERROR "In order to build the APBS Python shared library, BUILD_SHARED_LIBS must be enabled.")
#elseif(ENABLE_PYTHON AND APPLE AND BUILD_SHARED_LIBS)
#             message(FATAL_ERROR "ENABLE_PYTHON option in OSX systems requires the variable BUILD_SHARED_LIBS to be set to off.")
#endif()


################################################################################
# Handle conditional building with OpenMP                                      #
################################################################################
option(ENABLE_OPENMP "Enable OpenMP parallelism" OFF)

if(ENABLE_OPENMP)
    if(NOT ENABLE_DEBUG)
        message(STATUS "Checking for OpenMP")
        find_package(OpenMP)
        if(OPENMP_FOUND)
            message(STATUS "OpenMP support enabled")
            add_definitions("${OpenMP_C_FLAGS}")
            list(APPEND APBS_LIBS ${OpenMP_C_FLAGS})
        else()
            message(FATAL_ERROR "OpenMP was not found.  OpenMP support disabled")
        endif()
    else()
        message(WARNING "OpenMP may not be enabled in debugging mode")
    endif()
endif()


################################################################################
# Handle library checks for embedded unix environments in windows              #
################################################################################
if(MINGW)
    message(STATUS "Checking for wsock32 in MinGW environment")
    find_library(MINGW_WSOCK32
                 NAMES wsock32
                 PATHS ${SYS_LIBPATHS}
                 DOC   "The wsock32 library")

    if(MINGW_WSOCK32)
        message(STATUS "The wsock32 library was found: ${MINGW_HAS_WSOCK32}")
    else()
        message(FATAL_ERROR "The wsock32 library was not fond")
    endif()
endif()

if(CYGWIN)
    message(STATUS "Checking for wsock32 in Cygwin environment")
    find_library(CYGWIN_WSOCK32
                 NAMES wsock32
                 PATHS ${SYS_LIBPATHS}
                 DOC   "The wsock32 library")

    if(CYGWIN_WSOCK32)
        message(STATUS "The wsock32 library was found: ${CYGWIN_WSOCK32}")
        list(APPEND APBS_LIBS ${CYGWIN_WSOCK32})
    else()
        message(FATAL_ERROR "The wsock32 library was not fond")
    endif()

    set(HAVE_CYGWIN 1)
endif()

if(NOT CYGWIN AND NOT MINGW AND WIN32)
    list(APPEND APBS_LIBS wsock32 WS2_32)
    ADD_DEFINITIONS("/D _CRT_SECURE_NO_WARNINGS")
endif()


################################################################################
# Build APBS sources                                                           #
################################################################################
link_directories(${APBS_LIB_DIRS})
include_directories(${APBS_ROOT}/src ${APBS_BUILD}/src ${APBS_ROOT}/include)
add_subdirectory(src)


################################################################################
# Build APBS documentation                                                     #
################################################################################
option(BUILD_DOC "Build/Rebuild documentation" ON)

if(BUILD_DOC)
    message(STATUS "Building documentation enabled")
    add_subdirectory(docs)
endif()


################################################################################
# Handle conditional building with verbose debugging information printouts     #
################################################################################
option(BUILD_TOOLS "Build supplemental tools" ON)

if(BUILD_TOOLS)
    message(STATUS "Supplemental tools enabled")
    add_subdirectory(tools)
endif()


################################################################################
# Set up additional directories to install                                     #
################################################################################
install(DIRECTORY ${APBS_ROOT}/docs
        DESTINATION ${SHARE_INSTALL_PATH}
        PATTERN "programmer" EXCLUDE
        PATTERN "CMakeLists.txt" EXCLUDE)

install(DIRECTORY ${APBS_ROOT}/examples
        DESTINATION ${SHARE_INSTALL_PATH})

install(DIRECTORY ${APBS_ROOT}/tests
        DESTINATION ${SHARE_INSTALL_PATH}
        FILES_MATCHING
        PATTERN "*.py"
        PATTERN "README")

if(BUILD_TOOLS)
    install(DIRECTORY ${APBS_ROOT}/tools
            DESTINATION ${SHARE_INSTALL_PATH}
            USE_SOURCE_PERMISSIONS
            PATTERN "CMakeLists.txt"
            PATTERN "matlab"      EXCLUDE
            PATTERN "__pycache__" EXCLUDE
            PATTERN "*.py[co]"    EXCLUDE)

    if(NOT "${APBS_ROOT}" STREQUAL "${APBS_BUILD}")
    install(DIRECTORY ${APBS_BUILD}/tools/bin
            DESTINATION ${SHARE_INSTALL_PATH}/tools
            USE_SOURCE_PERMISSIONS)
    endif()
endif()

################################################################################
# Optionally build iAPBS interface                                             #
################################################################################
option(ENABLE_iAPBS "Enable iAPBS" OFF)

if(ENABLE_iAPBS)
    message(STATUS "Building of iAPBS interface enabled")
    add_subdirectory(contrib/iapbs/src)
    list(APPEND APBS_LIBS ${iAPBS_LIBRARY})
endif()


################################################################################
# Clean up temporary files and directories                                     #
################################################################################
file(REMOVE_RECURSE ${EXECUTABLE_OUTPUT_PATH}/temp)


###############################################################################
#####  Do Testing
###############################################################################
option(ENABLE_TESTS "Enable tests" ON)
if(${ENABLE_TESTS})
  enable_testing()
  add_subdirectory(tests)
endif()

###############################################################################
#####  Generate Packages
###############################################################################
# string(TOLOWER ${PROJECT_NAME}            PACKAGE_NAME)
set(PACKAGE_NAME                       "${PROJECT_NAME}")
set(CPACK_PACKAGE_NAME                 "${PACKAGE_NAME}")
set(CPACK_PACKAGE_DESCRIPTION          "APBS - Adaptive Poisson Boltzmann Solver")

set(CPACK_RESOURCE_FILE_LICENSE        "${PROJECT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_DESCRIPTION_FILE     "${PROJECT_SOURCE_DIR}/README.md")

set(CPACK_SOURCE_IGNORE_FILES          "${PROJECT_BINARY_DIR};/.git/;.gitignore;/tools/matlab/")
set(CPACK_PACKAGE_VERSION_MAJOR        "${CMAKE_MAJOR_VERSION}")
set(CPACK_PACKAGE_VERSION_MINOR        "${CMAKE_MINOR_VERSION}")
set(CPACK_PACKAGE_VERSION_PATCH        "${CMAKE_MICRO_VERSION}")
#  NOTE: The following is tempting but does not work!
# set(FULL_PACKAGE_NAME                  "${PACKAGE_NAME}-${MAJOR_VERSION}.${MINOR_VERSION}.${MICRO_VERSION}")
set(CPACK_PACKAGE_VENDOR               "PNNL")
set(CPACK_PACKAGE_CONTACT              "APBS - https://github.com/Electrostatics/apbs")
set(CPACK_PACKAGE_FILE_NAME            "${PACKAGE_NAME}-${APBS_VERSION}.${CMAKE_HOST_SYSTEM_NAME}")
set(CPACK_STRIP_FILES                   True)

list(APPEND CPACK_GENERATOR            "ZIP")
set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY    True)

if(CPACK_GENERATOR MATCHES "^(7Z|TBZ2|TGZ|TXZ|TZ|TZST|ZIP)$")
    set(CPACK_TOPLEVEL_TAG ${CPACK_PACKAGE_NAME})
    set(CPACK_SET_DESTDIR  False)
    set(CPACK_COMPONENT_INCLUDE_TOPLEVEL_DIRECTORY False)
endif()

include (CPack)
