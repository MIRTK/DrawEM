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


[ $# -eq 1 ] || { echo "usage: $(basename "$0") <subject>" 1>&2; exit 1; }
subj=$1

rdir=posteriors
sdir=segmentations-data


subcorts=(`cat $DRAWEMDIR/parameters/subcortical-all.csv`)
all=(`cat $DRAWEMDIR/parameters/all-labels.csv `)
numsubcorts=${#subcorts[*]}


mkdir -p $sdir/posteriors/temp || exit 1
for ((n=0;n<${#all[*]};n++));do r=${all[$n]}; mkdir -p $rdir/seg$r || exit 1;done


# out, csf
cp $sdir/posteriors/csf/$subj.nii.gz $rdir/seg24/$subj.nii.gz || exit 1

for tissue in csf gm wm;do 
  mkdir -p $rdir/$tissue
  cp $sdir/posteriors/$tissue/$subj.nii.gz $rdir/$tissue/$subj.nii.gz || exit 1
done

for tissue in gm wm;do 
    # cortical wm, gm
    corts=(`cat $DRAWEMDIR/parameters/cortical-$tissue.csv`)
    numcorts=${#corts[*]}

    addem=""
    for ((n=0;n<$numcorts;n++));do
        r=${corts[$n]};
        addem=$addem"-add $sdir/labels/seg$r-extended/$subj.nii.gz ";
    done
    addem=`echo $addem|sed -e 's:^-add::g'`
    run mirtk calculate $addem -out $sdir/posteriors/temp/$subj-$tissue.nii.gz

    for ((n=0;n<$numcorts;n++));do 
        r=${corts[$n]};
        run mirtk calculate $sdir/labels/seg$r-extended/$subj.nii.gz -div-with-zero $sdir/posteriors/temp/$subj-$tissue.nii.gz -mul $sdir/posteriors/$tissue/$subj.nii.gz -out $rdir/seg$r/$subj.nii.gz   
    done 
    rm $sdir/posteriors/temp/$subj-$tissue.nii.gz
done

# subcortical
for ((n=0;n<$numsubcorts;n++));do 
    r=${subcorts[$n]};
    cp $sdir/posteriors/seg$r/$subj.nii.gz $rdir/seg$r/$subj.nii.gz || exit 1
done 
