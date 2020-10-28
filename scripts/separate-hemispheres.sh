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


[ $# -eq 1 ] || { echo "usage: $(basename "$0") <subject>"; exit 1; }
subj=$1

scriptdir=$(dirname "$BASH_SOURCE")


if [ ! -f segmentations/${subj}_L_white.nii.gz -o ! -f segmentations/${subj}_R_white.nii.gz -o ! -f segmentations/${subj}_L_pial.nii.gz -o ! -f segmentations/${subj}_R_pial.nii.gz ];then
    # create initial labels files
    suffix=-sephemi
    $scriptdir/postprocess.sh $subj $suffix

    # left, right, both hemispheres
    run mirtk calculate segmentations/"$subj"_labels$suffix.nii.gz -mul 0 -out segmentations/$subj-empty.nii.gz
    run mirtk padding segmentations/$subj-empty.nii.gz segmentations/"$subj"_labels$suffix.nii.gz segmentations/$subj-L-hemisphere.nii.gz `echo $LEFT_HEMI_LABELS|wc -w` `echo $LEFT_HEMI_LABELS` 1
    run mirtk padding segmentations/$subj-empty.nii.gz segmentations/"$subj"_labels$suffix.nii.gz segmentations/$subj-R-hemisphere.nii.gz `echo $RIGHT_HEMI_LABELS|wc -w` `echo $RIGHT_HEMI_LABELS` 1

    # compute a cutting plane based on dmap 
    for h in L R;do
        run mirtk calculate-distance-map segmentations/$subj-$h-hemisphere.nii.gz segmentations/$subj-$h-hemisphere-dmap.nii.gz
        run mirtk smooth-image segmentations/$subj-$h-hemisphere-dmap.nii.gz segmentations/$subj-$h-hemisphere-dmap.nii.gz 1 -float
    done
    run mirtk calculate segmentations/$subj-R-hemisphere-dmap.nii.gz -sub segmentations/$subj-L-hemisphere-dmap.nii.gz -mask-below 0 -inside 1 -outside 0 -out segmentations/$subj-hemisphere-cut.nii.gz

    for surf in white pial;do
        # white / pial surfaces
        surf_tissue_labels=${surf^^}_SURFACE_TISSUE_LABELS
        num_labels=`echo ${!surf_tissue_labels}|wc -w`
        run mirtk padding segmentations/"$subj"_tissue_labels$suffix.nii.gz segmentations/"$subj"_tissue_labels$suffix.nii.gz segmentations/${subj}_${surf}_init.nii.gz $num_labels ${!surf_tissue_labels} 0 -invert $num_labels ${!surf_tissue_labels} 1

        h_num=0
        for h in L R;do
            # left, right white / pial surfaces
            run mirtk padding segmentations/${subj}_${surf}_init.nii.gz segmentations/$subj-hemisphere-cut.nii.gz segmentations/${subj}_${h}_${surf}_init.nii.gz $h_num 0
            # keep only connected components with volume > 5% volume of first component 
            $scriptdir/clear-small-components.sh segmentations/${subj}_${h}_${surf}_init.nii.gz segmentations/${subj}_${h}_${surf}_unfilled.nii.gz
            # fill holes
            run mirtk fill-holes segmentations/${subj}_${h}_${surf}_unfilled.nii.gz segmentations/${subj}_${h}_${surf}.nii.gz
            let h_num++
        done
    done

    # clean up
    for h in L R;do
        for surf in white pial;do
            rm segmentations/${subj}_${h}_${surf}_init.nii.gz segmentations/${subj}_${h}_${surf}_unfilled.nii.gz
        done
        rm segmentations/$subj-$h-hemisphere.nii.gz segmentations/$subj-$h-hemisphere-dmap.nii.gz
    done
    rm segmentations/$subj-hemisphere-cut.nii.gz segmentations/$subj-empty.nii.gz segmentations/${subj}_white_init.nii.gz segmentations/${subj}_pial_init.nii.gz
    rm segmentations/"$subj"_all_labels$suffix.nii.gz segmentations/"$subj"_labels$suffix.nii.gz segmentations/"$subj"_tissue_labels$suffix.nii.gz
fi
