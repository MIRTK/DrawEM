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


if [ ! -f segmentations/$subj-initial.nii.gz ];then
sdir=segmentations-data

mkdir -p segmentations $sdir/posteriors logs || exit 1;
for r in ${NONCORTICAL};do mkdir -p $sdir/posteriors/seg$r || exit 1; done
for str in ${ATLAS_TISSUES};do mkdir -p $sdir/posteriors/$str || exit 1; done


structs=""; saveposts=""; posts=""; num=0;
# subcortical
for r in ${NONCORTICAL};do 
structs="$structs $sdir/labels/seg$r/$subj.nii.gz";
post=$sdir/posteriors/seg$r/$subj.nii.gz
posts="$posts $post "
saveposts="$saveposts -saveprob $num $post "; 
num=$(($num+1)); 
done
# tissues
tissue_num=$num
tissues_parameter=""
for group in ${ATLAS_TISSUES_GROUPING};do
    tissues_parameter="$tissues_parameter $group"
    for ((group_i=0; group_i<$group; group_i++));do
        tissues_parameter="$tissues_parameter $tissue_num"
        tissue_num=$(($tissue_num+1))
    done
done

for str in ${ATLAS_TISSUES};do
structs="$structs $sdir/labels/$str/$subj.nii.gz";
post=$sdir/posteriors/$str/$subj.nii.gz
posts="$posts $post "
saveposts="$saveposts -saveprob $num $post ";
num=$(($num+1));
done


# segmentation
run mirtk draw-em N4/$subj.nii.gz $num $structs segmentations/$subj-em.nii.gz -padding 0 -mrf $CONNECTIVITIES -tissues $tissues_parameter -hui -postpenalty $sdir/MADs/$subj-subspace.nii.gz $saveposts 1>logs/$subj-em 2>logs/$subj-em-err

# add posterior probability of sub-tissues to tissues
atlas_tissues_arr=($ATLAS_TISSUES)
tissue_num=0
for group in ${ATLAS_TISSUES_GROUPING};do
    addem=""
    out=$sdir/posteriors/${atlas_tissues_arr[$tissue_num]}/$subj.nii.gz
    for ((group_i=0; group_i<$group; group_i++));do
        addem=$addem"-add $sdir/posteriors/${atlas_tissues_arr[$tissue_num]}/$subj.nii.gz ";
        tissue_num=$(($tissue_num+1))
    done
    if [ $group -gt 1 ];then
        addem=`echo $addem|sed -e 's:^-add::g'`
        run mirtk calculate $addem -out $out
    fi
done

# posteriors [0,1] -> [0,100]
for post in ${posts};do
run mirtk calculate $post -mul 100 -out $post 
done

run mirtk padding segmentations/$subj-em.nii.gz segmentations/$subj-em.nii.gz segmentations/$subj-initial.nii.gz $EM_LABELS_TO_ALL_LABELS
fi
