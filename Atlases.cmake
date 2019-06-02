cmake_minimum_required(VERSION 2.7)

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

set(ATLAS_URL "https://biomedic.doc.ic.ac.uk/brain-development/downloads/dHCP/atlases-dhcp-structural-pipeline-v1.zip")
set(ATLAS_MD5 77e924bc17a4906f5814874009f5eca6)

if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/atlases)
    ExternalProject_Add(atlases
      URL ${ATLAS_URL}
      URL_MD5 ${ATLAS_MD5}
      PREFIX atlases
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ""
    )

    add_custom_target(atlases_move ALL
      ${CMAKE_COMMAND} -E rename ${CMAKE_CURRENT_BINARY_DIR}/atlases/src/atlases ${CMAKE_CURRENT_SOURCE_DIR}/atlases
    )

    ADD_DEPENDENCIES(atlases_move atlases)
endif()
