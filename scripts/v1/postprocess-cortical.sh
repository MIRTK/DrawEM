#!/bin/bash
# ============================================================================
# Developing brain Region Annotation With Expectation-Maximization (Draw-EM)
#
# Copyright 2013-2016 Imperial College London
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


[ $# -eq 1 ] || { echo "usage: $(basename "$0") <subject>" 1>&2; exit 1; }

subj=$1
sdir=segmentations-data

run(){
  echo "$@"
  "$@" || exit 1
}

corts=(`cat $DRAWEMDIR/parameters/cortical.csv`)
wmcorts=(`cat $DRAWEMDIR/parameters/cortical-wm.csv`)
numcorts=${#corts[*]}

mkdir -p $sdir/cortical || exit 1
for ((n=0;n<$numcorts;n++));do 
  r=${corts[$n]}; 
  mkdir -p $sdir/labels/seg$r-extended || exit 1
done


#max prob of cortical structures + mrf regularization
segnum=0; labels=""; structs="";
for ((n=0;n<$numcorts;n++));do 
  let segnum++
  r=${corts[$n]};
  w=${wmcorts[$n]};
  #merge wm,gm probs - we'll split them later using the segmentation posteriors 
  run mirtk calculate $sdir/labels/seg$w/$subj.nii.gz -add $sdir/labels/seg$r/$subj.nii.gz -out $sdir/labels/seg$r-extended/$subj.nii.gz;
  labels="$labels $sdir/labels/seg$r-extended/$subj.nii.gz "; 
  structs="$structs $sdir/labels/seg$r-extended/$subj.nii.gz "; 
  segnumbers="$segnumbers 1 $segnum $r" 
done

run mirtk ems-hard-segmentation $segnum $labels $sdir/cortical/$subj.nii.gz -mrftimes 1 -posteriors $structs || exit 1
run mirtk padding $sdir/cortical/$subj.nii.gz $sdir/cortical/$subj.nii.gz $sdir/cortical/$subj.nii.gz $segnumbers || exit 1




