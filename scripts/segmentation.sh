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

sdir=segmentations-data

structures=""; save_posteriors=""; posteriors=""; num=-1;
tissues_parameter=""
em_labels_to_all_labels=""

add_structure()
{
    structure=$1
    mkdir -p $sdir/posteriors/$structure || exit 1;
    num=$(($num+1)); 
    structures="$structures $sdir/labels/$structure/$subj.nii.gz";
    posterior=$sdir/posteriors/$structure/$subj.nii.gz
    posteriors="$posteriors $posterior "
    save_posteriors="$save_posteriors -saveprob $num $posterior "; 
}

add_noncortical()
{
    structure_nr=$1
    add_structure seg$structure_nr;
    em_labels_to_all_labels="$em_labels_to_all_labels 1 $(($num+1)) $structure_nr"
}

add_tissue(){
    tissue_var=$1
    output_all_label=$2
    tissue_labels=${!tissue_var}
    num_subtissues=`echo $tissue_labels |wc -w`
    tissues_parameter="$tissues_parameter $num_subtissues"
    em_labels_to_all_labels="$em_labels_to_all_labels $num_subtissues"
    for subtissue in $tissue_labels;do
        add_structure $subtissue;
        tissues_parameter="$tissues_parameter $num"
        em_labels_to_all_labels="$em_labels_to_all_labels $(($num+1))"
    done
    em_labels_to_all_labels="$em_labels_to_all_labels $output_all_label"
}

add_tissue_posteriors(){
    # add posterior probability of sub-tissues to tissues
    tissue_var=$1
    output_tissue=$2
    tissue_labels=${!tissue_var}
    addem=""
    for subtissue in $tissue_labels;do
        addem=$addem"-add $sdir/posteriors/$subtissue/$subj.nii.gz ";
    done
    addem=`echo $addem|sed -e 's:^-add::g'`
    run mirtk calculate $addem -out $sdir/posteriors/$output_tissue/$subj.nii.gz
}


if [ ! -f segmentations/$subj-initial.nii.gz ];then
    # subcortical
    for r in ${NONCORTICAL};do
        add_noncortical $r 
    done
    # tissues
    add_tissue OUTLIER_TISSUES $OUTLIER_LABEL
    add_tissue CSF_TISSUES $CSF_LABEL
    add_tissue GM_TISSUES $SUPER_GM_LABEL
    add_tissue WM_TISSUES $SUPER_WM_LABEL

    # segmentation
    mkdir -p segmentations || exit 1;
    num_structures=$(($num+1))
    run mirtk draw-em N4/$subj.nii.gz $num_structures $structures segmentations/$subj-em.nii.gz -padding 0 -mrf $CONNECTIVITIES -tissues $tissues_parameter -hui -postpenalty $sdir/MADs/$subj-subspace.nii.gz $save_posteriors 1>logs/$subj-em 2>logs/$subj-em-err

    # posteriors [0,1] -> [0,100]
    for post in ${posteriors};do
        run mirtk calculate $post -mul 100 -out $post 
    done

    # add posterior probability of sub-tissues to tissues
    add_tissue_posteriors OUTLIER_TISSUES outlier
    add_tissue_posteriors CSF_TISSUES csf
    add_tissue_posteriors GM_TISSUES gm
    add_tissue_posteriors WM_TISSUES wm

    run mirtk padding segmentations/$subj-em.nii.gz segmentations/$subj-em.nii.gz segmentations/$subj-initial.nii.gz $em_labels_to_all_labels
fi
