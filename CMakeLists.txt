cmake_minimum_required(VERSION 3.18)

project(libjupitermag VERSION 1.5.0 LANGUAGES C CXX)

include(GNUInstallDirs)
include(CTest)
include(FetchContent)

option(LIBJUPITERMAG_BUILD_TESTS "Build sample test executables" ON)
option(LIBJUPITERMAG_BUILD_TIMING_EXE "Build legacy timing executable" OFF)
option(BUILD_SHARED_LIBS "Build shared libraries instead of static" ON)

set(LIBJUPITERMAG_LIBCON2020_GIT_TAG "a55701a7c8a94aa5a2e2a15d6a77e8674539b3ff" CACHE STRING "Git tag/branch for libcon2020 when using FetchContent")
set(LIBJUPITERMAG_LIBSPLINE_GIT_TAG "0d3c05b291a3887541331423e0d2cbbcfe4c7bca" CACHE STRING "Git tag/branch for libspline when using FetchContent")
set(LIBJUPITERMAG_LIBINTERNALFIELD_GIT_TAG "46008b41ab92ce39f789df44dfafe9072a2cc753" CACHE STRING "Git tag/branch for libinternalfield when using FetchContent")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Common include roots used throughout the legacy source tree.
set(LIBJUPITERMAG_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(LIBJUPITERMAG_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(LIBINTERNALFIELD_BUILD_TESTS OFF CACHE BOOL "Disable libinternalfield tests in superbuild" FORCE)

FetchContent_Declare(
    libcon2020
    GIT_REPOSITORY https://github.com/mattkjames7/libcon2020.git
    GIT_TAG ${LIBJUPITERMAG_LIBCON2020_GIT_TAG}
    SOURCE_SUBDIR src
)
FetchContent_Declare(
    libspline
    GIT_REPOSITORY https://github.com/mattkjames7/libspline.git
    GIT_TAG ${LIBJUPITERMAG_LIBSPLINE_GIT_TAG}
    SOURCE_SUBDIR src
)
FetchContent_Declare(
    libinternalfield
    GIT_REPOSITORY https://github.com/mattkjames7/libinternalfield.git
    GIT_TAG ${LIBJUPITERMAG_LIBINTERNALFIELD_GIT_TAG}
)

FetchContent_MakeAvailable(libcon2020 libspline libinternalfield)

set(LIBCON2020_ROOT ${libcon2020_SOURCE_DIR})
set(LIBSPLINE_ROOT ${libspline_SOURCE_DIR})
set(LIBINTERNALFIELD_ROOT ${libinternalfield_SOURCE_DIR})

add_subdirectory(src)

if(LIBJUPITERMAG_BUILD_TESTS)
    add_subdirectory(test)
endif()

install(FILES include/jupitermag.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
