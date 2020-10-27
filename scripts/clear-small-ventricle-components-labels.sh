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

# cleaning up small ventricle components
high_wm_em_label=0
for label in $NONCORTICAL $OUTLIER_TISSUES $CSF_TISSUES $GM_TISSUES $WM_TISSUES;do
    if [ $label == $HIGH_WM_TISSUE ];then break; fi
    let high_wm_em_label++
done

num_ventricles=`echo $VENTRICLES|wc -w`
run mirtk padding segmentations/$subj-em.nii.gz segmentations/$subj-em.nii.gz $sdir/corrections/$subj-hwm-init.nii.gz $high_wm_em_label 0 -invert
run mirtk padding segmentations/$subj-initial.nii.gz segmentations/$subj-initial.nii.gz $sdir/corrections/$subj-ven-init.nii.gz $num_ventricles $VENTRICLES 0 -invert $num_ventricles $VENTRICLES 1
$scriptdir/clear-small-components.sh $sdir/corrections/$subj-ven-init.nii.gz $sdir/corrections/$subj-ven.nii.gz
run mirtk calculate $sdir/corrections/$subj-ven-init.nii.gz -sub $sdir/corrections/$subj-ven.nii.gz -out $sdir/corrections/$subj-ven-diff.nii.gz

# if small ventricle components are surrounded by hwm..
run mirtk fill-holes-nn-based $sdir/corrections/$subj-hwm-init.nii.gz $sdir/corrections/$subj-ven-diff.nii.gz $sdir/corrections/$subj-hwm-fillh.nii.gz
run mirtk calculate $sdir/corrections/$subj-ven-diff.nii.gz -mul $sdir/corrections/$subj-hwm-fillh.nii.gz -out $sdir/corrections/$subj-ventohwm.nii.gz

volcorr=`mirtk measure-volume $sdir/corrections/$subj-ventohwm.nii.gz`
if [ "$volcorr" != "" ];then 
    # ..redistribute small ventricle components' probability into hwm/wm
    addem=""
    for label in $VENTRICLES;do
        addem=$addem"-add $sdir/posteriors/seg$label/$subj.nii.gz ";
    done
    addem=`echo $addem|sed -e 's:^-add::g'`
    run mirtk calculate $addem -mul $sdir/corrections/$subj-ventohwm.nii.gz -out $sdir/corrections/$subj-ventohwm-prob.nii.gz
    for tissue in hwm wm;do
        run mirtk calculate $sdir/corrections/$subj-ventohwm-prob.nii.gz -add $sdir/posteriors/$tissue/$subj.nii.gz -out $sdir/posteriors/$tissue/$subj.nii.gz
    done
    for label in $VENTRICLES;do
        run mirtk padding $sdir/posteriors/seg$label/$subj.nii.gz $sdir/corrections/$subj-ventohwm.nii.gz $sdir/posteriors/seg$label/$subj.nii.gz 1 0
    done

    # ..fix the segmentation too
    run mirtk padding segmentations/$subj-initial.nii.gz $sdir/corrections/$subj-ventohwm.nii.gz segmentations/$subj-initial.nii.gz 1 $SUPER_WM_LABEL

    # clean up
    rm $sdir/corrections/$subj-ventohwm-prob.nii.gz
fi

# clean up
rm $sdir/corrections/$subj-ven.nii.gz $sdir/corrections/$subj-ven-diff.nii.gz $sdir/corrections/$subj-hwm-init.nii.gz $sdir/corrections/$subj-hwm-fillh.nii.gz $sdir/corrections/$subj-ven-init.nii.gz
