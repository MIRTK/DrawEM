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


[ $# -ge 1 ] || { echo "usage: $BASH_SOURCE <subject>" 1>&2; exit 1; }
subj=$1

scriptdir=$(dirname "$BASH_SOURCE")
sdir=segmentations-data
mkdir -p $sdir/corrections || exit 1

run mirtk calculate $sdir/posteriors/gm/$subj.nii.gz -div 100 -mul $sdir/tissue-posteriors/$TISSUE_ATLAS_GMAT/$subj.nii.gz -out $sdir/corrections/$subj-gm-in-gmat-prob.nii.gz
run mirtk calculate $sdir/posteriors/gm/$subj.nii.gz -sub $sdir/corrections/$subj-gm-in-gmat-prob.nii.gz -out $sdir/posteriors/gm/$subj.nii.gz
run mirtk calculate $sdir/posteriors/wm/$subj.nii.gz -add $sdir/corrections/$subj-gm-in-gmat-prob.nii.gz -out $sdir/posteriors/wm/$subj.nii.gz

run mirtk calculate $sdir/posteriors/wm/$subj.nii.gz -sub $sdir/posteriors/gm/$subj.nii.gz -binarize 0 -mul segmentations/$subj-initial.nii.gz -out $sdir/corrections/$subj-wm-gt-gm-labels.nii.gz
run mirtk padding segmentations/$subj-initial.nii.gz $sdir/corrections/$subj-wm-gt-gm-labels.nii.gz segmentations/$subj-initial.nii.gz 1 $SUPER_GM_LABEL $SUPER_WM_LABEL
