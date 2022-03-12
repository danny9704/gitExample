## -------------------------------------------------------------------------- ##
#
#   Project: FortranProjectBase
#
#   Description:
#       Project example to begin the computer package and to educate student.
#
#   Author:
#       Young-Myung CHOI
#         Dept. of Naval Arch. and Ocean Eng.,
#         Pusan National University, Republic of Korea
#
## -------------------------------------------------------------------------- ##

## Default Settings --------------------------------------------------------- ##

##... Set the minimum version of cmake
cmake_minimum_required(VERSION 2.9 FATAL_ERROR)

##... CMake Project Name
project(FortranProjectBase)

##... Set Default Fortran Compiler & Eable the Fortran Compilation Environment
enable_language(Fortran)
if (NOT DEFINED CMAKE_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER gfortran)
endif()

## CMake Configuration ------------------------------------------------------ ##

##... Load Default CMake Configuration
include( config/config.cmake )

##... Set Fortran Compiler Flags
include( config/fortranFlag.cmake )

##... Auto Compilation Functions
include( config/Func_AutoCompile.cmake )

## Print CMake Project Information ------------------------------------------ ##
message("")
message("${BoldYellow}Project: ${BoldWhite}${CMAKE_PROJECT_NAME}")
message("")
message("  - ${Green}Projct Path       ${White}: ${PROJECT_COMPILE_PATH} ${Green}")
message("  - ${Green}GIT_BRANCH        ${White}: ${GIT_BRANCH} ${Green}")
message("  - ${Green}GIT_COMMIT_HASH   ${White}: ${GIT_COMMIT_HASH} ${Green}")
message("")
message("  - ${Green}Lib. Install Path ${White}: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${Green}")
message("  - ${Green}Exe. Install Path ${White}: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} ${Green}")
message("")

if( ${FORTRAN_COMPILE_MODE} STREQUAL "Release" )
message("  - ${Green}Compilation Mode  ${White}: Release ${Green}")
elseif( ${FORTRAN_COMPILE_MODE} STREQUAL "Debug" )
message("  - ${Green}Compilation Mode  ${White}: Debug ${Green}")
endif()

message("${Green}")
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
message("  - ${Green}Fortran Compiler  ${White}: Intel Fortran ${Green}")
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
message("  - ${Green}Fortran Compiler  ${White}: GNU Fortran ${Green}")
endif()
message("  - ${Green}Global Flag       ${White}: ${CMAKE_Fortran_FLAGS} ${Green}")
message("${ColourReset}")
## Compile Library ---------------------------------------------------------- ##

#set(COMPILE_DYNAMIC_LIB False)
set(COMPILE_DYNAMIC_LIB True)

#... Library 1
set( LIB_NAME_lib1 lib1)
# set(LIB_PATH_lib1 )
# set(LIB_DEPEND_lib1 )
# set(LIB_EXTLIB_lib1 )
# set(LIB_EXTINC_lib1 )
# set(LIB_EXTFLAG_lib1 )
add_subdirectory( src/lib1 )

#... Library 2
set( LIB_NAME_lib2 lib2)
# set(LIB_PATH_lib2 )
set( LIB_DEPEND_lib2 lib1 )
# set(LIB_EXTLIB_lib2 )
# set(LIB_EXTINC_lib2 )
# set(LIB_EXTFLAG_lib2 )
add_subdirectory(src/lib2)

## Compile Test Programs ---------------------------------------------------- ##

#... Test 1
set( TEST_NAME_testLib1 testLib1)
# set(TEST_PATH_testLib1 )
set( TEST_DEPEND_testLib1 lib1 )
# set(TEST_EXTLIB_testLib1 )
# set(TEST_EXTINC_testLib1 )
# set(TEST_EXTFLAG_testLib1 )
add_subdirectory(tests/testLib1)

#... Test 2
set( TEST_NAME_testLib2 testLib2)
# set(TEST_PATH_testLib2 )
set( TEST_DEPEND_testLib2 lib1 lib2 )
# set(TEST_EXTLIB_testLib2 )
# set(TEST_EXTINC_testLib2 )
# set(TEST_EXTFLAG_testLib2 )
add_subdirectory(tests/testLib2)

## Compile Applications ----------------------------------------------------- ##

#... Application 1
set( APP_NAME_app1 app1)
# set(APP_PATH_app1 )
set( APP_DEPEND_app1 lib1 lib2 )
# set(APP_EXTLIB_app1 )
# set(APP_EXTINC_app1 )
# set(APP_EXTFLAG_app1 )
add_subdirectory(applications/app1)
