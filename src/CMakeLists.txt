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

set(HEADERS
  ${BINARY_INCLUDE_DIR}/mirtk/${PROJECT_NAME}Export.h
  BiasCorrection.h
  BiasField.h
  BSplineBiasField.h
  DrawEM.h
  EMBase.h
  HashProbabilisticAtlas.h
  ImageHistogram1D.h
  Gaussian.h
  KMeans.h
  MeanShift.h
  NormalizeNyul.h
  PolynomialBiasField.h
  ProbabilisticAtlas.h
)

set(SOURCES
  BiasCorrection.cc
  BiasField.cc
  BSplineBiasField.cc
  DrawEM.cc
  EMBase.cc
  HashProbabilisticAtlas.cc
  ImageHistogram1D.cc
  Gaussian.cc
  KMeans.cc
  MeanShift.cc
  NormalizeNyul.cc
  PolynomialBiasField.cc
  ProbabilisticAtlas.cc
)

set(DEPENDS
  LibCommon
  LibNumerics
  LibImage
  LibIO
)

mirtk_add_library(AUTO_REGISTER)
