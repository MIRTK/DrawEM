#!/bin/bash
# ============================================================================
# Developing brain Region Annotation With Expectation-Maximization (Draw-EM)
#
# Copyright 2013-2016 Imperial College London
# Copyright 2013-2016 Andreas Schuh
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


[ $# -ge 1 -a $# -le 2 ] || { echo "usage: $(basename "$0") <subject> [<suffix>]" 1>&2; exit 1; }

subj=$1
suffix=""
[ $# -eq 1 ] || suffix=$2

sdir=segmentations-data
scriptdir=$(dirname "$BASH_SOURCE")

f=segmentations/$subj-initial.nii.gz


run(){
  echo "$@"
  "$@" || exit 1
}

if [ ! -f segmentations/"$subj"_all_labels$suffix.nii.gz ];then
# creating the all labels file (initial segmentation + cortical division)
$scriptdir/postprocess-cortical.sh $subj
run mirtk padding $sdir/cortical/$subj.nii.gz $f segmentations/$subj-cortical-wm.nii.gz 2000 0 -invert 
run mirtk padding $sdir/cortical/$subj.nii.gz $f segmentations/$subj-cortical-gm.nii.gz 1000 0 -invert 

corts=(`cat $DRAWEMDIR/parameters/cortical.csv`)
wmcorts=(`cat $DRAWEMDIR/parameters/cortical-wm.csv`)
gmtowmnumbers=""
for ((r=0;r<${#corts[*]};r++));do
gmtowmnumbers="$gmtowmnumbers 1 ${corts[$r]} ${wmcorts[$r]}"
done
run mirtk padding segmentations/$subj-cortical-wm.nii.gz segmentations/$subj-cortical-wm.nii.gz segmentations/$subj-cortical-wm.nii.gz $gmtowmnumbers 
run mirtk padding $f $f segmentations/"$subj"_all_labels_ini$suffix.nii.gz 2 1000 2000 0
run mirtk calculate segmentations/"$subj"_all_labels_ini$suffix.nii.gz -add segmentations/$subj-cortical-gm.nii.gz -add segmentations/$subj-cortical-wm.nii.gz -out segmentations/"$subj"_all_labels$suffix.nii.gz
# cleanup
rm segmentations/$subj-cortical-gm.nii.gz segmentations/$subj-cortical-wm.nii.gz segmentations/"$subj"_all_labels_ini$suffix.nii.gz
fi

if [ ! -f segmentations/"$subj"_labels$suffix.nii.gz ];then
# creating the labels file
run mirtk padding segmentations/"$subj"_all_labels$suffix.nii.gz segmentations/"$subj"_all_labels$suffix.nii.gz segmentations/"$subj"_labels$suffix.nii.gz $DRAWEMDIR/parameters/all-labels-to-labels.csv  
fi

if [ ! -f segmentations/"$subj"_tissue_labels$suffix.nii.gz ];then
# creating the tissue labels labels
run mirtk padding segmentations/"$subj"_all_labels$suffix.nii.gz segmentations/"$subj"_all_labels$suffix.nii.gz segmentations/"$subj"_tissue_labels$suffix.nii.gz $DRAWEMDIR/parameters/all-labels-to-tissue-labels.csv  
fi

