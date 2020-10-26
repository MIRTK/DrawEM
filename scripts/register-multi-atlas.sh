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


[ $# -ge 2 -a $# -le 3 ] || { echo "usage: $(basename "$0") <subject> <age> <#jobs>" 1>&2; exit 1; }
subj=$1
age=$2
njobs=1
if [ $# -gt 2 ];then njobs=$3;fi

sdir=segmentations-data
mkdir -p dofs

# initial tissue segmentation if we do multi-channel registration
if [ $MULTICHANNEL_REGISTRATION -eq 1 ];then
    need_tissue_seg=0
    for atlas in ${ATLASES};do
      if [ ! -f dofs/$subj-$atlas-n.dof.gz ];then
        need_tissue_seg=1
        break
      fi
    done

    if [ $need_tissue_seg -eq 1 ];then
        script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
        $script_dir/tissue-priors.sh $subj $age $njobs
    fi
fi

# registration of multiple atlases
for atlas in ${ATLASES};do
  dof=dofs/$subj-$atlas-n.dof.gz
  if [ ! -f $dof ];then
    if [ $MULTICHANNEL_REGISTRATION -eq 1 ];then
        run mirtk register N4/$subj.nii.gz $ATLAS_T2_DIR/$atlas.nii.gz $sdir/gm-posteriors/$subj.nii.gz $ATLAS_GM_POSTERIORS_DIR/$atlas.nii.gz -parin $DRAWEMDIR/parameters/ireg-multichannel-structural.cfg  -dofout $dof -threads $njobs -v 0
    else
        run mirtk register N4/$subj.nii.gz $ATLAS_T2_DIR/$atlas.nii.gz -parin $DRAWEMDIR/parameters/ireg.cfg  -dofout $dof -threads $njobs -v 0
    fi
  fi
done

