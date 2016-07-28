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


if [ -n "$FSLDIR" ]; then
  [ -f "$FSLDIR/bin/bet" ] || { echo "FSLDIR environment variable invalid!" 1>&2; exit 1; }
else
  FSLDIR="$(cd "$(dirname "$(which bet)")"/.. && pwd)"
  [ -f "$FSLDIR/bin/bet" ] || { echo "FSLDIR environment variable not set!" 1>&2; exit 1; }
  export PATH="$FSLDIR/bin:$PATH"
fi

run(){
  echo "$@"
  "$@" || exit 1
}


sdir=segmentations-data

mkdir -p segmentations N4 dofs bias || exit 1

if [ ! -f N4/$subj.nii.gz ];then 
  #convert image and rescale
  run mirtk convert-image T2/$subj.nii.gz N4/${subj}_rescaled.nii.gz -rescale 0 1000 -double 

  if [ ! -f segmentations/${subj}_brain_mask.nii.gz ];then
    #brain extract
    run bet N4/${subj}_rescaled.nii.gz segmentations/${subj}_brain.nii.gz -R -f 0.1 -m 
    rm segmentations/${subj}_brain.nii.gz
  fi

  #bias correct
  run $DRAWEMDIR/ThirdParty/ITK/N4 3 -i N4/${subj}_rescaled.nii.gz -x segmentations/${subj}_brain_mask.nii.gz -o "[N4/${subj}_corr.nii.gz,bias/$subj.nii.gz]" -c "[50x50x50,0.001]" -s 2 -b "[100,3]" -t "[0.15,0.01,200]"
  run mirtk calculate N4/${subj}_corr.nii.gz -mul segmentations/${subj}_brain_mask.nii.gz -out N4/${subj}_corr.nii.gz 
  
  #rescale image
  run mirtk convert-image N4/${subj}_corr.nii.gz N4/$subj.nii.gz -rescale 0 1000 -double 

  rm N4/${subj}_rescaled.nii.gz N4/${subj}_corr.nii.gz
fi
