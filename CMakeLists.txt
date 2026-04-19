cmake_minimum_required(VERSION 3.18)

project(libjupitermag VERSION 1.6.2 LANGUAGES C CXX)

include(GNUInstallDirs)
include(CTest)
include(FetchContent)

option(LIBJUPITERMAG_BUILD_TESTS "Build sample test executables" ON)
option(LIBJUPITERMAG_BUILD_TIMING_EXE "Build legacy timing executable" OFF)
option(BUILD_SHARED_LIBS "Build shared libraries instead of static" ON)

set(LIBJUPITERMAG_LIBCON2020_GIT_TAG "187cb26764eac5223d90dc57c69e68a0aee273e2" CACHE STRING "Git tag/branch for libcon2020 when using FetchContent")  # v1.2.1
set(LIBJUPITERMAG_LIBSPLINE_GIT_TAG "38f91bacb6a4e8f2d9598aa460ed94bc0966783e" CACHE STRING "Git tag/branch for libspline when using FetchContent")  # v0.3.1
set(LIBJUPITERMAG_LIBINTERNALFIELD_GIT_TAG "9a6b6e562b57a95146b7d608ee9cdf229b674a63" CACHE STRING "Git tag/branch for libinternalfield when using FetchContent")  # v1.4.1

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Common include roots used throughout the legacy source tree.
set(LIBJUPITERMAG_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(LIBJUPITERMAG_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(LIBINTERNALFIELD_BUILD_TESTS OFF CACHE BOOL "Disable libinternalfield tests in superbuild" FORCE)

# Build all FetchContent dependencies as static libraries to avoid
# platform-specific shared-library runtime lookup issues in CI.
set(LIBJUPITERMAG_LOCAL_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries instead of static" FORCE)

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

# Restore the originally requested build mode for libjupitermag and
# subsequent local targets.
set(BUILD_SHARED_LIBS ${LIBJUPITERMAG_LOCAL_BUILD_SHARED_LIBS} CACHE BOOL "Build shared libraries instead of static" FORCE)

set(LIBCON2020_ROOT ${libcon2020_SOURCE_DIR})
set(LIBSPLINE_ROOT ${libspline_SOURCE_DIR})
set(LIBINTERNALFIELD_ROOT ${libinternalfield_SOURCE_DIR})

add_subdirectory(src)

if(LIBJUPITERMAG_BUILD_TESTS)
    add_subdirectory(test)
endif()

install(
    FILES
        include/jupitermag.h
        ${LIBCON2020_ROOT}/include/con2020.h
        ${LIBSPLINE_ROOT}/include/spline.h
        ${LIBINTERNALFIELD_ROOT}/include/internalfield.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
