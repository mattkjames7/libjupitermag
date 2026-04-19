include(GoogleTest)

set(REGRESSION_COMMON_SOURCE regression_common.cc)
set(REGRESSION_TEST_SOURCES
    regression_config_gtest.cc
    regression_field_gtest.cc
    regression_geometry_gtest.cc
    regression_robustness_gtest.cc
    regression_trace_gtest.cc
)

foreach(regression_test_source ${REGRESSION_TEST_SOURCES})
    get_filename_component(regression_test_name ${regression_test_source} NAME_WE)

    add_executable(${regression_test_name}
        ${regression_test_source}
        ${REGRESSION_COMMON_SOURCE}
    )

    target_link_libraries(${regression_test_name}
        PRIVATE
            jupitermag
            GTest::gtest_main
    )

    target_include_directories(${regression_test_name}
        PRIVATE
            ${LIBJUPITERMAG_INCLUDE_DIR}
            ${CMAKE_CURRENT_SOURCE_DIR}
    )

    gtest_discover_tests(
        ${regression_test_name}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..
        DISCOVERY_MODE PRE_TEST
    )
endforeach()
