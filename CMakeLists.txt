cmake_minimum_required (VERSION 3.0)

project ( HADMOLEE LANGUAGES CXX Fortran)

set(CMAKE_CXX_FLAGS "-fPIC -O3") 
set(CMAKE_BUILD_TYPE "Release")

# Set output lib and bin locations
set(CMAKE_INSTALL_BINDIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_INSTALL_BINDIR})

# Make sure executables are always built with the install RPATH 
# else things get linked to the .so in the build directory and not the one in lib
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_CURRENT_SOURCE_DIR}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_CURRENT_SOURCE_DIR}/lib")
endif("${isSystemDir}" STREQUAL "-1")

# Optional toggle to turn off generating cling executable
option(BUILD_CLING "Build cling executable and ROOT Dictionary" "FALSE")

##-----------------------------------------------------------------------
## ROOT finding and connecting

# Make sure we can find ROOT 
execute_process(COMMAND root-config --prefix 
                COMMAND tr -d \\n
                RESULT_VARIABLE RC_RESULT 
                OUTPUT_VARIABLE ROOT_PREFIX)
list(APPEND CMAKE_MODULE_PATH "${ROOT_PREFIX}/share/root/cmake")

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS MathMore)
if (ROOT_FOUND)
    message("-- ROOT found! (${ROOT_USE_FILE})")
    #---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
    include(${ROOT_USE_FILE})
    list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}) 
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
    include_directories(${ROOT_INCLUDE_DIRS})
else()
    message(SEND_ERROR "-- Root not found")
endif()

##-----------------------------------------------------------------------
## Look for jpacStyle to build the cling executable 

# Find the jpacStyle library
find_library(JSTYLELIB NAMES jpacStyle libjpacStyle 
                       HINTS "$ENV{JPACSTYLE}/lib")
if (JSTYLELIB)
    include_directories($ENV{JPACSTYLE}/include)
    message("-- jpacStyle found! (${JSTYLELIB})")
else()
    message("-- jpacStyle not found!")
endif()

##-----------------------------------------------------------------------
# Find the cubature  library
find_library(CUBATURELIB NAMES cubature libcubature
                         HINTS "$ENV{CUBATURE}/lib")
if (CUBATURELIB) 
    include_directories($ENV{CUBATURE}/lib)
    message("-- cubature found! (${CUBATURELIB})")
else()
    message("-- cubature not found!")
endif()

##-----------------------------------------------------------------------
# Find the LoopTools library
find_library(LOOPLIB NAMES LoopTools libLoopTools
                     HINTS "$ENV{LOOPTOOLS}/lib")
if (LOOPLIB)
    include_directories($ENV{LOOPTOOLS}/include)
    message("-- LoopTools found! (${LOOPLIB})")
else()
    message("-- LoopTools not found!")
endif()

##-----------------------------------------------------------------------
# Build core physics library 

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/core)
if (LOOPLIB)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/LoopTools)
    add_subdirectory(src/LoopTools)
endif()
add_subdirectory(src/core)

if (BUILD_CLING)
    message("-- Building cling executable.")
    # Build running and loading executable
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/cling)
    add_subdirectory(src/cling)
else()
    message("-- Not building cling executable.")
endif()