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

export MULTICHANNEL_REGISTRATION=1
export HIGH_WM_VENTRICLE_CORRECTION=1
export HEMISPHERE_HOLE_CORRECTION=1

export CONNECTIVITIES=$DRAWEMDIR/parameters/ALBERT/connectivities.mrf
export LOOKUP_TABLE=$DRAWEMDIR/parameters/ALBERT/LUT.txt
export ATLAS_T2_DIR=$DRAWEMDIR/atlases/ALBERTs/T2
export ATLAS_TISSUES_DIR=$DRAWEMDIR/atlases/ALBERTs/tissues-v4
export ATLAS_SEGMENTATIONS_DIR=$DRAWEMDIR/atlases/ALBERTs/segmentations-v4
export ATLAS_GM_POSTERIORS_DIR=$DRAWEMDIR/atlases/ALBERTs/gm-posteriors-v3
export ATLASES=`ls $ATLAS_T2_DIR | sed -e 's:.nii.gz::g'`
export ATLAS_TISSUES="csf gm wm outlier hwm lwm"
export OUTLIER_TISSUES=outlier
export CSF_TISSUES=csf
export GM_TISSUES=gm
export WM_TISSUES="wm hwm lwm"
export HIGH_WM_TISSUE=hwm

####### The following are needed for the dHCP pipeline surface reconstruction #######
export CSF_TISSUE_LABEL=1
export GM_TISSUE_LABEL=2
export WM_TISSUE_LABEL=3
export BG_TISSUE_LABEL=4

export HIPPOCAMPI="1 2"
export AMYGDALA="3 4"
export DEEP_SUBCORTICAL_GM="40 41 42 43 44 45 46 47 85 86 87"
export DEEP_GM="$HIPPOCAMPI $AMYGDALA $DEEP_SUBCORTICAL_GM"
export CORPUS_CALLOSUM="48"
export BRAINSTEM="19"
export CEREBELLUM="17 18"
#####################################################################################

export CSF_LABEL=83
export OUTLIER_LABEL=84
export CORTICAL_WM="51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82"
export CORTICAL_GM="5 6 7 8 9 10 11 12 13 14 15 16 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39"
export CORTICAL="$CORTICAL_GM $CORTICAL_WM"
export VENTRICLES="49 50"
export NONCORTICAL="$DEEP_GM $VENTRICLES $BRAINSTEM $CEREBELLUM $CORPUS_CALLOSUM"
export NONCORTICAL="$NONCORTICAL 88 89 90"
export ALL_LABELS="$OUTLIER_LABEL $CSF_LABEL $CORTICAL $NONCORTICAL"
export LEFT_HEMI_LABELS="1 3 5 7 9 11 13 15 17 21 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 64 66 68 70 72 74 76 78 80 82 87"
export RIGHT_HEMI_LABELS="2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 50 52 54 56 58 60 62 63 65 67 69 71 73 75 77 79 81 86"
export ALL_LABELS_TO_TISSUE_LABELS="1 83 1 32 5 6 7 8 9 10 11 12 13 14 15 16 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 2 33 48 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 3 1 84 4 2 49 50 5 2 17 18 6 11 40 41 42 43 44 45 46 47 85 86 87 7 1 19 8 4 1 2 3 4 9"

if [ $TISSUE_ATLAS_NAME == "fetal" ];then
    # we will add cavum as a correction using the fetal atlas
    export SEPARATE_CAVUM_FROM_CSF_CORRECTION=1
    export CAVUM=91
    export ALL_LABELS="$ALL_LABELS $CAVUM"
    export ALL_LABELS_TO_LABELS="1 83 1 36 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 2 33 48 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 3 3 44 45 85 3 1 49 4 1 50 5 1 $CAVUM 6 1 19 7 1 17 8 1 18 9 1 90 10 2 41 47 11 2 40 46 12 2 43 87 13 2 42 86 14 1 88 15 1 89 16 1 48 17 1 84 0"
    export ALL_LABELS_TO_TISSUE_LABELS="$ALL_LABELS_TO_TISSUE_LABELS 1 $CAVUM 1"
else
    export ALL_LABELS_TO_LABELS="1 51 5 1 52 6 1 53 7 1 54 8 1 55 9 1 56 10 1 57 11 1 58 12 1 59 13 1 60 14 1 61 15 1 62 16 1 63 20 1 64 21 1 65 22 1 66 23 1 67 24 1 68 25 1 69 26 1 70 27 1 71 28 1 72 29 1 73 30 1 74 31 1 75 32 1 76 33 1 77 34 1 78 35 1 79 36 1 80 37 1 81 38 1 82 39 3 83 84 85 0 1 86 42 1 87 43"
fi

export WHITE_SURFACE_TISSUE_LABELS="$WM_TISSUE_LABEL 5 7 9"
export PIAL_SURFACE_TISSUE_LABELS="$WHITE_SURFACE_TISSUE_LABELS $GM_TISSUE_LABEL"
