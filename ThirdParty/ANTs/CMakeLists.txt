cmake_minimum_required(VERSION 2.7)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0046 NEW)
  cmake_policy(SET CMP0048 NEW)
endif(COMMAND cmake_policy)

PROJECT( N4 )

IF(NOT ITK_DIR)
  INCLUDE(External_ITK.cmake)
  SET(External_ITK 1)
  MESSAGE(STATUS "Build ITK as external project for build of DrawEM::N4")
ELSE()
  MESSAGE(STATUS "Use system-wide ITK installation for build of DrawEM::N4")
ENDIF()

# Set up ITK
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
    "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

INCLUDE_DIRECTORIES(${ITK_INCLUDE_DIRS})

SET(SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/antsCommandLineParser" "${CMAKE_CURRENT_SOURCE_DIR}/antsCommandLineOption" )

ADD_EXECUTABLE( N4 N4.cxx ${SOURCES} )
TARGET_LINK_LIBRARIES(N4 ${ITK_LIBRARIES})
set_target_properties(N4 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${MIRTK_TOOLS_DIR}")

IF(External_ITK)
  ADD_DEPENDENCIES(N4 ITK)
ENDIF()
