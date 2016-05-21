# ============================================================================
# Developing brain Region Annotation With Expectation-Maximization (Draw-EM)
#
# Copyright 2013-2016 Imperial College London
# Copyright 2013-2016 Christian Ledig
# Copyright 2013-2016 Antonios Makropoulos
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ============================================================================

##############################################################################
# @file  CMakeLists.txt
# @brief Build configuration of MIRTK Draw-EM tools.
##############################################################################

macro(add_image_command cmd)
  mirtk_add_executable(${cmd})
  mirtk_target_dependencies(${cmd} LibCommon LibNumerics LibImage LibIO ${ARGN})
endmacro()

add_image_command(calculate-filtering)
add_image_command(calculate-gradients)
add_image_command(change-label)
add_image_command(measure-volume)
add_image_command(padding)


macro(add_drawem_command cmd)
  add_image_command(${cmd} LibDrawEM)
endmacro()

add_drawem_command(fill-holes)
add_drawem_command(fill-holes-nn-based)
add_drawem_command(em)
add_drawem_command(draw-em)
add_drawem_command(em-hard-segmentation)
add_drawem_command(kmeans)
add_drawem_command(normalize)
add_drawem_command(split-labels)
