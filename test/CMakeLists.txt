add_executable(cpptest cpptest.cc)
target_link_libraries(cpptest PRIVATE jupitermag)
target_include_directories(cpptest PRIVATE ${LIBJUPITERMAG_INCLUDE_DIR})

add_executable(ctest ctest.c)
set_target_properties(ctest PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(ctest PRIVATE jupitermag m)
target_include_directories(ctest PRIVATE ${LIBJUPITERMAG_INCLUDE_DIR})

if(LIBJUPITERMAG_BUILD_TIMING_EXE)
    add_executable(timing timing.cc timing.h)
    target_link_libraries(timing PRIVATE jupitermag)
    target_include_directories(timing PRIVATE ${LIBJUPITERMAG_INCLUDE_DIR})
endif()

add_test(NAME ctest_sample COMMAND ctest)
add_test(NAME cpptest_sample COMMAND cpptest)

set(INSTALL_GTEST OFF CACHE BOOL "Disable googletest installation for embedded builds" FORCE)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/external/googletest/CMakeLists.txt")
    add_subdirectory(external/googletest external/googletest-build EXCLUDE_FROM_ALL)
else()
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0
    )
    FetchContent_MakeAvailable(googletest)
endif()

add_subdirectory(unit)
