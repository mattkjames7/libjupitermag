set(LIBJUPITERMAG_SOURCES
    model.cc
    trace.cc
    interptraceclosestpos.cc
    footprint.cc
    coordconv.cc
    libjupitermag.cc
)

add_library(jupitermag ${LIBJUPITERMAG_SOURCES})
add_library(jupitermag::jupitermag ALIAS jupitermag)

set_target_properties(jupitermag PROPERTIES OUTPUT_NAME jupitermag)

if(WIN32)
    set_target_properties(jupitermag PROPERTIES
        PREFIX "lib"
        IMPORT_PREFIX "lib"
    )
endif()


# Only use .def file for MSVC, not MinGW/MSYS2
if(WIN32 AND BUILD_SHARED_LIBS)
    set_target_properties(jupitermag PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
    if(MSVC)
        target_sources(jupitermag PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/jupitermag.def)
    endif()
endif()

target_include_directories(jupitermag
    PUBLIC
        $<BUILD_INTERFACE:${LIBJUPITERMAG_INCLUDE_DIR}>
        $<BUILD_INTERFACE:${LIBCON2020_ROOT}/include>
        $<BUILD_INTERFACE:${LIBINTERNALFIELD_ROOT}/include>
        $<BUILD_INTERFACE:${LIBSPLINE_ROOT}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
        ${LIBCON2020_ROOT}/src
        ${LIBSPLINE_ROOT}/src
        ${LIBJUPITERMAG_SRC_DIR}
)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(jupitermag PRIVATE -Wextra)
endif()

target_link_libraries(jupitermag
    PUBLIC
        con2020::con2020
        internalfield
        spline::spline
)

if (UNIX AND NOT APPLE)
    target_link_libraries(jupitermag PUBLIC m)
endif()

install(TARGETS jupitermag
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
