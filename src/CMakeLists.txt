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

target_include_directories(jupitermag
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

target_compile_options(jupitermag PRIVATE -Wextra)

target_link_libraries(jupitermag
    PUBLIC
        con2020::con2020
        internalfield
        spline::spline
        m
)

install(TARGETS jupitermag
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
