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


[ $# -eq 1 ] || { echo "usage: $(basename "$0") <subject>" 1>&2; exit 1; }

subj=$1
sdir=segmentations-data

run(){
  echo "$@"
  "$@" || exit 1
}

for tissue in gm wm;do
    cortical_var=CORTICAL_${tissue^^}
    cortical_labels=${!cortical_var}

    mkdir -p $sdir/cortical-$tissue || exit 1
    for label in $cortical_labels;do 
        mkdir -p $sdir/labels/seg$label-extended || exit 1
    done

    #max prob of cortical structures + mrf regularization
    segnum=0; labels=""; structs="";
    for label in $cortical_labels;do 
        let segnum++
        labels="$labels $sdir/labels/seg$label/$subj.nii.gz "
        structs="$structs $sdir/labels/seg$label-extended/$subj.nii.gz "
        segnumbers="$segnumbers 1 $segnum $label" 
    done

    run mirtk em-hard-segmentation $segnum $labels $sdir/cortical-$tissue/$subj.nii.gz -mrftimes 1 -posteriors $structs
    run mirtk padding $sdir/cortical-$tissue/$subj.nii.gz $sdir/cortical-$tissue/$subj.nii.gz $sdir/cortical-$tissue/$subj.nii.gz $segnumbers
done
