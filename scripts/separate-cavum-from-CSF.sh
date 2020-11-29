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
mkdir -p $sdir/corrections $sdir/posteriors/seg$CAVUM || exit 1

addem=""
for subtissue in $TISSUE_ATLAS_CAVUM $TISSUE_ATLAS_CSF_TISSUES $TISSUE_ATLAS_VENTRICLES;do
    addem=$addem"-add $sdir/tissue-posteriors/$subtissue/$subj.nii.gz "
done
addem=`echo $addem|sed -e 's:^-add::g'`
run mirtk calculate $addem -out $sdir/corrections/$subj-ventricles-csf-tissue-posterior.nii.gz
run mirtk calculate $sdir/tissue-posteriors/$TISSUE_ATLAS_CAVUM/$subj.nii.gz -div $sdir/corrections/$subj-ventricles-csf-tissue-posterior.nii.gz -out $sdir/corrections/$subj-cavum-part.nii.gz

structures="csf"
for ven in $VENTRICLES;do
    structures="$structures seg$ven"
done

addem=""
for structure in $structures;do
    addem=$addem"-add $sdir/posteriors/$structure/$subj.nii.gz "
done
addem=`echo $addem|sed -e 's:^-add::g'`
run mirtk calculate $addem -out $sdir/corrections/$subj-csf-posterior.nii.gz

run mirtk calculate $sdir/corrections/$subj-csf-posterior.nii.gz -mul $sdir/corrections/$subj-cavum-part.nii.gz -out $sdir/posteriors/seg$CAVUM/$subj.nii.gz
for structure in $structures;do
    run mirtk calculate $sdir/corrections/$subj-cavum-part.nii.gz -sub 1 -mul -1 -mul $sdir/posteriors/$structure/$subj.nii.gz -out $sdir/posteriors/$structure/$subj.nii.gz
done


num_no_cavum=`echo $structures|wc -w`
structures="$structures seg$CAVUM"
posteriors=""
for structure in $structures;do
    posteriors=$posteriors" $sdir/posteriors/$structure/$subj.nii.gz";
done
num=`echo $structures|wc -w`
run mirtk em-hard-segmentation $num $posteriors $sdir/corrections/$subj-ventricles-cavum-hard.nii.gz

run mirtk padding $sdir/corrections/$subj-ventricles-cavum-hard.nii.gz segmentations/$subj-initial.nii.gz $sdir/corrections/$subj-ventricles-cavum-hard.nii.gz $num_no_cavum $CSF_LABEL $VENTRICLES 0 -invert
run mirtk padding segmentations/$subj-initial.nii.gz $sdir/corrections/$subj-ventricles-cavum-hard.nii.gz segmentations/$subj-initial.nii.gz 1 $num $CAVUM

rm $sdir/corrections/$subj-ventricles-csf-tissue-posterior.nii.gz $sdir/corrections/$subj-cavum-part.nii.gz $sdir/corrections/$subj-csf-posterior.nii.gz $sdir/corrections/$subj-ventricles-cavum-hard.nii.gz
