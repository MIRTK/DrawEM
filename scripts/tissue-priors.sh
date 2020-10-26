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


[ $# -ge 2 ] || { echo "usage: $(basename "$0") <subject> <age> [<#jobs>]" 1>&2; exit 1; }
subj=$1
age=$2
njobs=1
if [ $# -gt 2 ];then njobs=$3;fi

sdir=segmentations-data

structures=""; save_posteriors=""; posteriors=""; num=-1;
tissues_parameter=""

add_structure()
{
    structure=$1
    mkdir -p $sdir/tissue-posteriors/$structure || exit 1
    num=$(($num+1))
    structures="$structures $sdir/template/$structure/$subj.nii.gz"
    posterior=$sdir/tissue-posteriors/$structure/$subj.nii.gz
    posteriors="$posteriors $posterior "
    save_posteriors="$save_posteriors -saveprob $num $posterior "
}

add_tissue(){
    tissue_var=$1
    tissue_labels=${!tissue_var}
    num_subtissues=`echo $tissue_labels |wc -w`
    tissues_parameter="$tissues_parameter $num_subtissues"
    for subtissue in $tissue_labels;do
        add_structure $subtissue;
        tissues_parameter="$tissues_parameter $num"
    done
}

if [ ! -f $sdir/gm-posteriors/$subj.nii.gz  ];then
    echo "creating $subj tissue priors"

    [ $age -lt $TEMPLATE_MAX_AGE ] || { age=$TEMPLATE_MAX_AGE; }
    [ $age -gt $TEMPLATE_MIN_AGE ] || { age=$TEMPLATE_MIN_AGE; }

    # registration of atlas template
    template_dof=dofs/$subj-tissue-atlas-$age-n.dof.gz
    if [ ! -f $template_dof ];then
        run mirtk register N4/$subj.nii.gz $TISSUE_ATLAS_T2_DIR/template-$age.nii.gz -dofout $template_dof -parin $DRAWEMDIR/parameters/ireg.cfg -threads $njobs -v 0
    fi

    mkdir -p $sdir/tissue-initial-segmentations $sdir/gm-posteriors || exit 1
    for tissue in ${TISSUE_ATLAS_TISSUES};do
        mkdir -p $sdir/template/$tissue || exit 1
        run mirtk transform-image $TISSUE_ATLAS_TISSUES_DIR/$tissue/$age.nii.gz $sdir/template/$tissue/$subj.nii.gz -dofin $template_dof -target N4/$subj.nii.gz -interp Linear
    done

    # subcortical
    for r in ${TISSUE_ATLAS_NONCORTICAL};do
        add_structure $r
    done
    # tissues
    add_tissue TISSUE_ATLAS_OUTLIER_TISSUES
    add_tissue TISSUE_ATLAS_CSF_TISSUES
    add_tissue TISSUE_ATLAS_GM_TISSUES
    add_tissue TISSUE_ATLAS_WM_TISSUES

    num_structures=$(($num+1))
    run mirtk draw-em N4/$subj.nii.gz $num_structures $structures $sdir/tissue-initial-segmentations/$subj.nii.gz -padding 0 -mrf $TISSUE_ATLAS_CONNECTIVITIES  -tissues $tissues_parameter -hui -relaxtimes 2 $save_posteriors  1>logs/$subj-tissue-em 2>logs/$subj-tissue-em-err

    addem=""
    for subtissue in ${TISSUE_ATLAS_GM_TISSUES};do
        addem=$addem"-add $sdir/tissue-posteriors/$subtissue/$subj.nii.gz ";
    done
    addem=`echo $addem|sed -e 's:^-add::g'`
    run mirtk calculate $addem -mul 100 -out $sdir/gm-posteriors/$subj.nii.gz
fi
