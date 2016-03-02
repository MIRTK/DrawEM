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


[ $# -eq 1 ] || { echo "usage: $(basename "$0") <subject>"; exit 1; }
subj=$1

run(){
  echo "$@"
  "$@" || exit 1
}

scriptdir=$(dirname "$BASH_SOURCE")


if [ ! -f segmentations/${subj}_L_white.nii.gz -o ! -f segmentations/${subj}_R_white.nii.gz -o ! -f segmentations/${subj}_L_pial.nii.gz -o ! -f segmentations/${subj}_R_pial.nii.gz ];then

# create initial labels files
suffix=-sephemi
$scriptdir/postprocess.sh $subj $suffix


# structures for the left, right hemisphere and those that extend to both
left=`cat $DRAWEMDIR/parameters/structures-left.csv`
right=`cat $DRAWEMDIR/parameters/structures-right.csv`
both=`cat $DRAWEMDIR/parameters/structures-both.csv`
numleft=`echo $left|wc -w`
numright=`echo $right|wc -w`
numboth=`echo $both|wc -w`

# left, right, both hemispheres
run mirtk calculate segmentations/"$subj"_labels$suffix.nii.gz -mul 0 -out segmentations/$subj-empty.nii.gz
run mirtk padding segmentations/$subj-empty.nii.gz segmentations/"$subj"_labels$suffix.nii.gz segmentations/$subj-L-hemisphere.nii.gz $numleft `echo $left` 1
run mirtk padding segmentations/$subj-empty.nii.gz segmentations/"$subj"_labels$suffix.nii.gz segmentations/$subj-R-hemisphere.nii.gz $numright `echo $right` 1
run mirtk padding segmentations/$subj-empty.nii.gz segmentations/"$subj"_labels$suffix.nii.gz segmentations/$subj-both-hemisphere.nii.gz $numboth `echo $both` 1
run mirtk padding segmentations/$subj-both-hemisphere.nii.gz segmentations/"$subj"_all_labels$suffix.nii.gz segmentations/$subj-both-hemisphere.nii.gz 85 1

# compute a cutting plane based on dmap 
run mirtk calculate-distance-map segmentations/$subj-L-hemisphere.nii.gz segmentations/$subj-L-hemisphere-dmap.nii.gz
run mirtk calculate-distance-map segmentations/$subj-R-hemisphere.nii.gz segmentations/$subj-R-hemisphere-dmap.nii.gz
run mirtk smooth-image segmentations/$subj-L-hemisphere-dmap.nii.gz segmentations/$subj-L-hemisphere-dmap.nii.gz 1 -float
run mirtk smooth-image segmentations/$subj-R-hemisphere-dmap.nii.gz segmentations/$subj-R-hemisphere-dmap.nii.gz 1 -float
run mirtk calculate segmentations/$subj-R-hemisphere-dmap.nii.gz -sub segmentations/$subj-L-hemisphere-dmap.nii.gz -mask-below 0 -inside 1 -outside 0 -out segmentations/$subj-hemisphere-cut.nii.gz

# split structures that are in both hemispheres into left and right according to the cutting plane
run mirtk calculate segmentations/$subj-hemisphere-cut.nii.gz -mul segmentations/$subj-both-hemisphere.nii.gz -add segmentations/$subj-L-hemisphere.nii.gz -out segmentations/$subj-L-hemisphere.nii.gz
run mirtk calculate segmentations/$subj-hemisphere-cut.nii.gz -sub 1 -mul -1 -mul segmentations/$subj-both-hemisphere.nii.gz -add segmentations/$subj-R-hemisphere.nii.gz -out segmentations/$subj-R-hemisphere.nii.gz
rm segmentations/$subj-both-hemisphere.nii.gz segmentations/$subj-hemisphere-cut.nii.gz segmentations/$subj-L-hemisphere-dmap.nii.gz  segmentations/$subj-R-hemisphere-dmap.nii.gz segmentations/$subj-empty.nii.gz

rmlabelswhite="5 1 2 4 6 8"
rmlabelspial="4 1 4 6 8"

for h in L R;do
    for surf in white pial;do
	# left, right wm+dgm / gm+wm+cgm
	rmlabels=rmlabels$surf
        run mirtk padding segmentations/$subj-$h-hemisphere.nii.gz segmentations/"$subj"_tissue_labels$suffix.nii.gz segmentations/${subj}_${h}_${surf}_init.nii.gz ${!rmlabels} 0

        # keep only connected components with volume > 5% volume of first component	
	$scriptdir/clear-small-components.sh segmentations/${subj}_${h}_${surf}_init.nii.gz segmentations/${subj}_${h}_${surf}_unfilled.nii.gz

	# fill holes
	run mirtk fill-holes segmentations/${subj}_${h}_${surf}_unfilled.nii.gz segmentations/${subj}_${h}_${surf}.nii.gz

        # clean up
        rm segmentations/${subj}_${h}_${surf}_init.nii.gz segmentations/${subj}_${h}_${surf}_unfilled.nii.gz
    done
    # clean up
    rm segmentations/$subj-$h-hemisphere.nii.gz
done

# clean up
rm segmentations/"$subj"_all_labels$suffix.nii.gz segmentations/"$subj"_labels$suffix.nii.gz segmentations/"$subj"_tissue_labels$suffix.nii.gz

fi


