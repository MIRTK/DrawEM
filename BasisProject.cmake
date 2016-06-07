################################################################################
# @file  BasisProject.cmake
# @brief Sets basic information about the MIRTK module and calls basis_project().
#
# This file defines basic information about a project by calling 
# the basis_project() function. This basic information, also known as metadata, 
# is used by CMake BASIS to setup the project. The dependencies to other modules
# have to be specified here such that the top-level IRTK project can analyze the
# inter-module dependencies, as well as dependencies on third-party libraries.
#
# @sa http://opensource.andreasschuh.com/cmake-basis/standard/modules.html
#
# @ingroup BasisSettings
################################################################################

# Note: The #<*> dependency patterns are required by the basisproject tool and
#       should be kept on a separate line as last commented argument of the
#       corresponding options of the basis_project() command. The TEMPLATE
#       option and set argument are also required by this tool and should not
#       be changed manually. The argument is updated by basisproject --update.

basis_project (

  # ----------------------------------------------------------------------------
  # meta-data
  NAME        "DrawEM"
  PACKAGE     "MIRTK"
  AUTHORS     "Antonios Makropoulos"
  DESCRIPTION "Tools for the segmentation of the developing brain"
  COPYRIGHT   "2016 Imperial College London, Antonios Makropoulos"
  LICENSE     "Apache License Version 2.0"
  CONTACT     "Antonios Makropoulos <a.makropoulos11@imperial.ac.uk>"
  TEMPLATE    "mirtk-module/1.0"

  # ----------------------------------------------------------------------------
  # dependencies
  DEPENDS
    MIRTK{Common,Numerics,Image}
    #<dependency>
  OPTIONAL_DEPENDS
    TBB{tbb}
    #<optional-dependency>
  TOOLS_DEPENDS
    MIRTK{IO}
  TEST_DEPENDS
    #<test-dependency>
  OPTIONAL_TEST_DEPENDS
    #<optional-test-dependency>

)
