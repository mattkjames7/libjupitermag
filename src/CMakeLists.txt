set(LIBJUPITERMAG_SOURCES
    model.cc
    trace.cc
    interptraceclosestpos.cc
    footprint.cc
    coordconv.cc
    libjupitermag.cc
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

install(TARGETS jupitermag_shared jupitermag_static
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
