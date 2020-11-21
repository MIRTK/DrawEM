#!/bin/bash
# ============================================================================
# Developing brain Region Annotation With Expectation-Maximization (Draw-EM)
#
# Copyright 2013-2020 Imperial College London
# Copyright 2013-2020 Andreas Schuh
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

rdir=posteriors
sdir=segmentations-data

mkdir -p $sdir/posteriors/temp|| exit 1
for label in $ALL_LABELS;do mkdir -p $rdir/seg$label $sdir/posteriors/seg$label || exit 1;done

cp $sdir/posteriors/csf/$subj.nii.gz $sdir/posteriors/seg$CSF_LABEL/$subj.nii.gz || exit 1
cp $sdir/posteriors/outlier/$subj.nii.gz $sdir/posteriors/seg$OUTLIER_LABEL/$subj.nii.gz || exit 1


for tissue in gm wm;do 
    # cortical wm, gm
    cortical_var=CORTICAL_${tissue^^}
    cortical_labels=${!cortical_var}

    addem=""
    for label in $cortical_labels;do 
        addem=$addem"-add $sdir/labels/seg$label-extended/$subj.nii.gz ";
    done
    addem=`echo $addem|sed -e 's:^-add::g'`
    run mirtk calculate $addem -out $sdir/posteriors/temp/$subj-$tissue-sum.nii.gz

    for label in $cortical_labels;do
        run mirtk calculate $sdir/labels/seg$label-extended/$subj.nii.gz -div-with-zero $sdir/posteriors/temp/$subj-$tissue-sum.nii.gz -mul $sdir/posteriors/$tissue/$subj.nii.gz -out $sdir/posteriors/seg$label/$subj.nii.gz   
    done 
    rm $sdir/posteriors/temp/$subj-$tissue-sum.nii.gz
done

# make sure all posteriors sum up to 100
addem=""
for label in $ALL_LABELS;do
    addem=$addem"-add $sdir/posteriors/seg$label/$subj.nii.gz ";
done
addem=`echo $addem|sed -e 's:^-add::g'`
run mirtk calculate $addem -out $sdir/posteriors/temp/$subj-sum.nii.gz

for label in $ALL_LABELS;do
    mirtk calculate $sdir/posteriors/seg$label/$subj.nii.gz -div-with-zero $sdir/posteriors/temp/$subj-sum.nii.gz -mul 100 -out $rdir/seg$label/$subj.nii.gz
done 
rm $sdir/posteriors/temp/$subj-sum.nii.gz

# copy tissues
for tissue in outlier csf gm wm;do
    mkdir -p $rdir/$tissue
    cp $sdir/posteriors/$tissue/$subj.nii.gz $rdir/$tissue/$subj.nii.gz || exit 1
done