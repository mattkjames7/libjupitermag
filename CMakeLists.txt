cmake_minimum_required(VERSION 3.16)

project(libjupitermag VERSION 1.5.0 LANGUAGES C CXX)

include(GNUInstallDirs)
include(CTest)
include(FetchContent)

option(LIBJUPITERMAG_BUILD_TESTS "Build sample test executables" ON)
option(LIBJUPITERMAG_USE_FETCHCONTENT "Fetch dependencies with CMake instead of using local submodules" OFF)
option(LIBJUPITERMAG_BUILD_TIMING_EXE "Build legacy timing executable" OFF)

set(LIBJUPITERMAG_LIBCON2020_GIT_TAG "a55701a7c8a94aa5a2e2a15d6a77e8674539b3ff" CACHE STRING "Git tag/branch for libcon2020 when using FetchContent")
set(LIBJUPITERMAG_LIBSPLINE_GIT_TAG "0d3c05b291a3887541331423e0d2cbbcfe4c7bca" CACHE STRING "Git tag/branch for libspline when using FetchContent")
set(LIBJUPITERMAG_LIBINTERNALFIELD_GIT_TAG "93733365e2b18f6fc24db731d5765ddda62b510d" CACHE STRING "Git tag/branch for libinternalfield when using FetchContent")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Common include roots used throughout the legacy source tree.
set(LIBJUPITERMAG_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(LIBJUPITERMAG_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

if(LIBJUPITERMAG_USE_FETCHCONTENT)
    set(LIBINTERNALFIELD_BUILD_TESTS OFF CACHE BOOL "Disable libinternalfield tests in superbuild" FORCE)

    FetchContent_Declare(
        libcon2020
        GIT_REPOSITORY https://github.com/mattkjames7/libcon2020.git
        GIT_TAG ${LIBJUPITERMAG_LIBCON2020_GIT_TAG}
    )
    FetchContent_Declare(
        libspline
        GIT_REPOSITORY https://github.com/mattkjames7/libspline.git
        GIT_TAG ${LIBJUPITERMAG_LIBSPLINE_GIT_TAG}
    )
    FetchContent_Declare(
        libinternalfield
        GIT_REPOSITORY https://github.com/mattkjames7/libinternalfield.git
        GIT_TAG ${LIBJUPITERMAG_LIBINTERNALFIELD_GIT_TAG}
    )

    FetchContent_GetProperties(libcon2020)
    if(NOT libcon2020_POPULATED)
        FetchContent_Populate(libcon2020)
    endif()
    add_subdirectory(${libcon2020_SOURCE_DIR}/src ${libcon2020_BINARY_DIR}/src EXCLUDE_FROM_ALL)

    FetchContent_GetProperties(libspline)
    if(NOT libspline_POPULATED)
        FetchContent_Populate(libspline)
    endif()
    add_subdirectory(${libspline_SOURCE_DIR}/src ${libspline_BINARY_DIR}/src EXCLUDE_FROM_ALL)

    FetchContent_GetProperties(libinternalfield)
    if(NOT libinternalfield_POPULATED)
        FetchContent_Populate(libinternalfield)
    endif()
    add_subdirectory(${libinternalfield_SOURCE_DIR} ${libinternalfield_BINARY_DIR} EXCLUDE_FROM_ALL)

    set(LIBCON2020_ROOT ${libcon2020_SOURCE_DIR})
    set(LIBSPLINE_ROOT ${libspline_SOURCE_DIR})
    set(LIBINTERNALFIELD_ROOT ${libinternalfield_SOURCE_DIR})
else()
    set(LIBINTERNALFIELD_BUILD_TESTS OFF CACHE BOOL "Disable libinternalfield tests in superbuild" FORCE)

    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/lib/libcon2020/CMakeLists.txt")
        message(FATAL_ERROR "Missing dependency at lib/libcon2020. Initialize submodules or set LIBJUPITERMAG_USE_FETCHCONTENT=ON.")
    endif()
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/lib/libspline/CMakeLists.txt")
        message(FATAL_ERROR "Missing dependency at lib/libspline. Initialize submodules or set LIBJUPITERMAG_USE_FETCHCONTENT=ON.")
    endif()
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/lib/libinternalfield/CMakeLists.txt")
        message(FATAL_ERROR "Missing dependency at lib/libinternalfield. Initialize submodules or set LIBJUPITERMAG_USE_FETCHCONTENT=ON.")
    endif()

    add_subdirectory(lib/libcon2020/src lib/libcon2020/build-cmake/src EXCLUDE_FROM_ALL)
    add_subdirectory(lib/libspline/src lib/libspline/build-cmake/src EXCLUDE_FROM_ALL)
    add_subdirectory(lib/libinternalfield lib/libinternalfield/build-cmake EXCLUDE_FROM_ALL)

    set(LIBCON2020_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/lib/libcon2020)
    set(LIBSPLINE_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/lib/libspline)
    set(LIBINTERNALFIELD_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/lib/libinternalfield)
endif()

set(LIBJUPITERMAG_SOURCES
    src/model.cc
    src/trace.cc
    src/interptraceclosestpos.cc
    src/footprint.cc
    src/coordconv.cc
    src/libjupitermag.cc
)

add_library(jupitermag_shared SHARED ${LIBJUPITERMAG_SOURCES})
add_library(jupitermag_static STATIC ${LIBJUPITERMAG_SOURCES})

set_target_properties(jupitermag_shared PROPERTIES OUTPUT_NAME jupitermag)
set_target_properties(jupitermag_static PROPERTIES OUTPUT_NAME jupitermag)

foreach(target jupitermag_shared jupitermag_static)
    target_include_directories(${target}
        PUBLIC
            $<BUILD_INTERFACE:${LIBJUPITERMAG_INCLUDE_DIR}>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        PRIVATE
            ${LIBCON2020_ROOT}/include
            ${LIBCON2020_ROOT}/src
            ${LIBINTERNALFIELD_ROOT}/include
            ${LIBSPLINE_ROOT}/include
            ${LIBSPLINE_ROOT}/src
            ${LIBJUPITERMAG_SRC_DIR}
    )

    target_compile_options(${target} PRIVATE -Wextra)

    target_link_libraries(${target}
        PUBLIC
            con2020::con2020
            internalfield
            spline::spline
            m
    )
endforeach()

if(LIBJUPITERMAG_BUILD_TESTS)
    add_executable(cpptest test/cpptest.cc)
    target_link_libraries(cpptest PRIVATE jupitermag_shared)
    target_include_directories(cpptest PRIVATE ${LIBJUPITERMAG_INCLUDE_DIR})

    add_executable(ctest test/ctest.c)
    set_target_properties(ctest PROPERTIES LINKER_LANGUAGE CXX)
    target_link_libraries(ctest PRIVATE jupitermag_shared m)
    target_include_directories(ctest PRIVATE ${LIBJUPITERMAG_INCLUDE_DIR})

    if(LIBJUPITERMAG_BUILD_TIMING_EXE)
        add_executable(timing test/timing.cc test/timing.h)
        target_link_libraries(timing PRIVATE jupitermag_shared)
        target_include_directories(timing PRIVATE ${LIBJUPITERMAG_INCLUDE_DIR})
    endif()

    add_test(NAME ctest_sample COMMAND ctest)
    add_test(NAME cpptest_sample COMMAND cpptest)
endif()

install(TARGETS jupitermag_shared jupitermag_static
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(FILES include/jupitermag.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
