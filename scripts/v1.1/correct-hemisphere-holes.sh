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


[ $# -eq 1 ] || { echo "usage: $BASH_SOURCE <subject>" 1>&2; exit 1; }
subj=$1

sdir=segmentations-data
mkdir -p $sdir/corrections || exit 1

subcorts=`cat $DRAWEMDIR/parameters/subcortical-all.csv`

run(){
  echo "$@"
  "$@" || exit 1
}

if [ -f segmentations/${subj}_L_white.nii.gz -a -f segmentations/${subj}_R_white.nii.gz ];then

# if small gm components exist inside the white surface..
run mirtk calculate segmentations/$subj-initial.nii.gz -mul 0 -out $sdir/corrections/$subj-gmtochange.nii.gz
for h in L R;do
run mirtk calculate segmentations/${subj}_${h}_white.nii.gz -mul segmentations/$subj-initial.nii.gz  -out $sdir/corrections/$subj-gmmask.nii.gz
run mirtk padding $sdir/corrections/$subj-gmtochange.nii.gz $sdir/corrections/$subj-gmmask.nii.gz $sdir/corrections/$subj-gmtochange.nii.gz 1000 1
done
rm $sdir/corrections/$subj-gmmask.nii.gz

volcorr=`mirtk measure-volume $sdir/corrections/$subj-gmtochange.nii.gz`
if [ "$volcorr" != "" ];then 

# ..redistribute small components' probability into wm/dgm
num=1;
inprobs="$sdir/posteriors/wm/$subj.nii.gz 2000"
outprobs="$sdir/posteriors/wm/$subj.nii.gz"
for r in ${subcorts}; do
inprobs="$inprobs $sdir/posteriors/seg$r/$subj.nii.gz $r"
outprobs="$outprobs $sdir/posteriors/seg$r/$subj.nii.gz"
let num=num+1
done

run mirtk change-label segmentations/$subj-initial.nii.gz $sdir/corrections/$subj-gmtochange.nii.gz $num $inprobs segmentations/$subj-initial.nii.gz $sdir/posteriors/gm/$subj.nii.gz $outprobs
run mirtk padding $sdir/posteriors/gm/$subj.nii.gz $sdir/corrections/$subj-gmtochange.nii.gz $sdir/posteriors/gm/$subj.nii.gz 1 0

fi

fi
