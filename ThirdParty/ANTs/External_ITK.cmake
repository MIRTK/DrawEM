cmake_minimum_required(VERSION 2.7)

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

ExternalProject_Add(
  ITK
  GIT_REPOSITORY "https://github.com/InsightSoftwareConsortium/ITK.git"
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/ThirdParty/ITK"
  CMAKE_ARGS -DBUILD_EXAMPLES=OFF -DBUILD_SHARED_LIBS=ON -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=${ITK_DIR}
  UPDATE_COMMAND ""
  PATCH_COMMAND ""
  INSTALL_COMMAND ""
)

ExternalProject_Get_Property(ITK SOURCE_DIR)
ExternalProject_Get_Property(ITK BINARY_DIR)
set(ITK_SOURCE_DIR ${SOURCE_DIR})
set(ITK_DIR ${BINARY_DIR})
