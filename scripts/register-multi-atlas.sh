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


[ $# -ge 2 -a $# -le 3 ] || { echo "usage: $(basename "$0") <subject> <age> <#jobs>" 1>&2; exit 1; }
subj=$1
age=$2
njobs=1
if [ $# -gt 2 ];then njobs=$3;fi

sdir=segmentations-data
mkdir -p dofs

# initial tissue segmentation if we do multi-channel registration
if [ "$ATLAS_REGISTRATION_CHANNELS" != "" ];then
    script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    $script_dir/tissue-priors.sh $subj $age $njobs

    energy=`cat $DRAWEMDIR/parameters/ireg-structural.cfg | sed 's/\[/@/g' | grep "Energy function ="`
    energy_channels=""
    img1=1
    img2=2
    for channel in $ATLAS_REGISTRATION_CHANNELS;do
        img1=$(($img1+2))
        img2=$(($img2+2))
        energy_channels="$energy_channels + 3.0 SSD[GM Sim](I($img1), I($img2) o T)"
    done
    new_energy="$energy $energy_channels"
    cat $DRAWEMDIR/parameters/ireg-structural.cfg | sed 's/\[/@/g'|sed -e "s:$energy:$new_energy:g" | sed 's/@/\[/g'> $sdir/ireg.cfg
else
    cp $DRAWEMDIR/parameters/ireg-structural.cfg $sdir/ireg.cfg
fi

# registration of multiple atlases
for atlas in ${ATLASES};do
    dof=dofs/$subj-$atlas-n.dof.gz
    if [ ! -f $dof ];then
        pairs="N4/$subj.nii.gz $ATLAS_T2_DIR/$atlas.nii.gz"
        if [ "$ATLAS_REGISTRATION_CHANNELS" != "" ];then
            for channel in $ATLAS_REGISTRATION_CHANNELS;do
                pairs="$pairs $sdir/tissue-posteriors/$channel/$subj.nii.gz $ATLAS_PROBABILITIES_DIR/$channel/$atlas.nii.gz"
            done
        fi
        run mirtk register $pairs -parin $sdir/ireg.cfg -dofout $dof -threads $njobs -v 0
    fi
done
