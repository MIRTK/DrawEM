set(ATLASES_URL "https://biomedic.doc.ic.ac.uk/brain-development/downloads/dHCP/atlases-dhcp-structural-pipeline-v1.zip")
set(ATLASES_MD5 "77e924bc17a4906f5814874009f5eca6")
set(ATLASES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/atlases")

set(_default ON)
if (IS_DIRECTORY "${ATLASES_DIR}")
  set(_default OFF)
endif ()
option(DrawEM_ATLASES "Whether to download DrawEM atlases during build step" ${_default})
unset(_default)

if (DrawEM_ATLASES)

  include(ExternalProject)

  ExternalProject_Add(atlases
    URL ${ATLASES_URL}
    URL_MD5 ${ATLASES_MD5}
    PREFIX atlases
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )

  ExternalProject_Add_Step(atlases rename
    COMMAND ${CMAKE_COMMAND} -E rename "${CMAKE_CURRENT_BINARY_DIR}/atlases/src/atlases" "${ATLASES_DIR}"
    COMMENT "Move extracted atlases to source directory"
    DEPENDEES download
  )

endif ()
