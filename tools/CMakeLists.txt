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


mirtk_add_executable(calculate-filtering  DEPENDS LibCommon LibNumerics LibImage LibIO)
mirtk_add_executable(calculate-gradients  DEPENDS LibCommon LibNumerics LibImage LibIO)
mirtk_add_executable(change-label  DEPENDS LibCommon LibNumerics LibImage LibIO)
mirtk_add_executable(fill-holes  DEPENDS LibCommon LibNumerics LibImage LibIO LibDrawEM)
mirtk_add_executable(fill-holes-nn-based  DEPENDS LibCommon LibNumerics LibImage LibIO LibDrawEM)
mirtk_add_executable(measure-volume  DEPENDS LibCommon LibNumerics LibImage LibIO)
mirtk_add_executable(padding  DEPENDS LibCommon LibNumerics LibImage LibIO)

mirtk_add_executable(em  DEPENDS LibCommon LibNumerics LibImage LibIO LibDrawEM)
mirtk_add_executable(draw-em  DEPENDS LibCommon LibNumerics LibImage LibIO LibDrawEM)
mirtk_add_executable(em-hard-segmentation  DEPENDS LibCommon LibNumerics LibImage LibIO LibDrawEM)
mirtk_add_executable(kmeans  DEPENDS LibCommon LibNumerics LibImage LibIO LibDrawEM)
mirtk_add_executable(normalize  DEPENDS LibCommon LibNumerics LibImage LibIO LibDrawEM)
mirtk_add_executable(split-labels  DEPENDS LibCommon LibNumerics LibImage LibIO LibDrawEM)
