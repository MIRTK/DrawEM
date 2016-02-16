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
  mirtkBiasCorrection.h
  mirtkBiasField.h
  mirtkBSplineBiasField.h
  mirtkDrawEM.h
  mirtkEMBase.h
  mirtkHashProbabilisticAtlas.h
  mirtkImageHistogram1D.h
  mirtkGaussian.h
  mirtkKMeans.h
  mirtkMeanShift.h
  mirtkNormalizeNyul.h
  mirtkPolynomialBiasField.h
  mirtkProbabilisticAtlas.h
)

set(SOURCES
  mirtkBiasCorrection.cc
  mirtkBiasField.cc
  mirtkBSplineBiasField.cc
  mirtkDrawEM.cc
  mirtkEMBase.cc
  mirtkHashProbabilisticAtlas.cc
  mirtkImageHistogram1D.cc
  mirtkGaussian.cc
  mirtkKMeans.cc
  mirtkMeanShift.cc
  mirtkNormalizeNyul.cc
  mirtkPolynomialBiasField.cc
  mirtkProbabilisticAtlas.cc
)

set(DEPENDS
  LibCommon
  LibNumerics
  LibImage
  LibPointSet
)

mirtk_add_library()


