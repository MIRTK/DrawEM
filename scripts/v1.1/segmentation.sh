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


[ $# -eq 1 ] || { echo "usage: $(basename "$0") <subject>"; exit 1; }
subj=$1

run()
{
  echo "$@"
  "$@" || exit 1
}

if [ ! -f segmentations/$subj-initial.nii.gz ];then
sdir=segmentations-data

subcorts=`cat $DRAWEMDIR/parameters/subcortical-all.csv`
tissues="outlier csf gm wm hwm lwm"


mkdir -p segmentations $sdir/posteriors logs || exit 1;
for r in ${subcorts};do mkdir -p $sdir/posteriors/seg$r || exit 1; done
for str in ${tissues};do mkdir -p $sdir/posteriors/$str || exit 1; done


structs=""; saveposts=""; posts=""; num=0;
# subcortical
for r in ${subcorts};do 
structs="$structs $sdir/labels/seg$r/$subj.nii.gz";
post=$sdir/posteriors/seg$r/$subj.nii.gz
posts="$posts $post "
saveposts="$saveposts -saveprob $num $post "; 
num=$(($num+1)); 
done
# tissues
for str in ${tissues};do
structs="$structs $sdir/labels/$str/$subj.nii.gz";
post=$sdir/posteriors/$str/$subj.nii.gz
posts="$posts $post "
saveposts="$saveposts -saveprob $num $post "; 
num=$(($num+1));
done



# segmentation
run mirtk draw-em N4/$subj.nii.gz 27 $structs segmentations/$subj-em.nii.gz -padding 0 -mrf $DRAWEMDIR/parameters/connectivities.mrf -tissues 1 21 1 22 1 23 3 24 25 26  -hui -postpenalty $sdir/MADs/$subj-subspace.nii.gz $saveposts 1>logs/$subj-em 2>logs/$subj-em-err

# add hwm and lwm posterior probability to wm
run mirtk calculate $sdir/posteriors/hwm/$subj.nii.gz -add $sdir/posteriors/lwm/$subj.nii.gz -add $sdir/posteriors/wm/$subj.nii.gz -out $sdir/posteriors/wm/$subj.nii.gz

# posteriors [0,1] -> [0,100]
for post in ${posts};do
run mirtk calculate $post -mul 100 -out $post 
done

# set whole wm to 2000 and whole gm to 1000 and final label numbers for the rest (subcortical, csf, out)
run mirtk padding segmentations/$subj-em.nii.gz segmentations/$subj-em.nii.gz segmentations/$subj-initial.nii.gz $DRAWEMDIR/parameters/seg-numbers.csv
fi
