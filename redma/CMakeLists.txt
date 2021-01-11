file(GLOB allfiles "${CMAKE_CURRENT_SOURCE_DIR}/*")

foreach(filename ${allfiles})
    if(IS_DIRECTORY ${filename})
        add_subdirectory(${filename})
    endif()
endforeach()

# note: setting this to shared may cause problems at runtime when loading meshes
add_library(core STATIC
            ${HEADERS}
            ${SOURCES})

target_link_libraries(core ${LINKED_LIBRARIES})
