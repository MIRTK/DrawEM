#!/bin/bash
# ============================================================================
# Developing brain Region Annotation With Expectation-Maximization (Draw-EM)
#
# Copyright 2013-2020 Imperial College London
# Copyright 2013-2020 Antonios Makropoulos
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

export TISSUE_ATLAS_CONNECTIVITIES=$DRAWEMDIR/parameters/fetal/connectivities.mrf
export TISSUE_ATLAS_T2_DIR=$DRAWEMDIR/atlases/fetal/T2
export TISSUE_ATLAS_TISSUES_DIR=$DRAWEMDIR/atlases/fetal/atlas-9
export TISSUE_ATLAS_TISSUES=`ls $TISSUE_ATLAS_TISSUES_DIR -v |grep structure`
export TISSUE_ATLAS_OUTLIER_TISSUES=structure4
export TISSUE_ATLAS_CSF_TISSUES="structure1"
export TISSUE_ATLAS_GM_TISSUES=structure2
export TISSUE_ATLAS_WM_TISSUES="structure3"
export TISSUE_ATLAS_DGM_TISSUES=structure7
export TISSUE_ATLAS_CAVUM=structure9
export TISSUE_ATLAS_VENTRICLES=structure5
export TISSUE_ATLAS_NONCORTICAL="$TISSUE_ATLAS_VENTRICLES structure6 structure7 structure8 $TISSUE_ATLAS_CAVUM"
export TISSUE_ATLAS_MIN_AGE=23
export TISSUE_ATLAS_MAX_AGE=37
export BET_THRESHOLD=0.1
