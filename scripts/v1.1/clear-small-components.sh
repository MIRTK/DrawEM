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

[ $# -ge 2 ] || { echo "usage: $BASH_SOURCE <input_mask> <output_mask> [<keep_percent_ratio>]" 1>&2; exit 1; }
f=$1
outf=$2
keepratio=0.05
if [ $# -gt 2 ];then keepratio=$3;fi

run(){
  echo "$@"
  "$@" || exit 1
}

if [ ! -f $outf ];then

sdir=segmentations-data
mkdir -p $sdir/corrections || exit 1
base=`basename $f| sed -e 's:.nii.gz::g'`

# connected components
run mirtk extract-connected-components $f $sdir/corrections/$base-lccs.nii.gz -all -output-component-labels -connectivity 6
mirtk measure-volume $sdir/corrections/$base-lccs.nii.gz > $sdir/corrections/$base-lccs-comps
comp1=`cat $sdir/corrections/$base-lccs-comps |grep "^1 "|cut -d' ' -f2`
comps=`cat $sdir/corrections/$base-lccs-comps |wc -l`
retainvol=`echo "$comp1*$keepratio" | /usr/bin/bc`
retainvol=`printf "%.*f\n" 0 $retainvol` #round
retain="1"

# keep only connected components with volume > 5% volume of first component
for ((r=2;r<$comps;r++));do 
    comp=`cat $sdir/corrections/$base-lccs-comps |grep "^$r "|cut -d' ' -f2`
    if [ `echo "$comp<$retainvol" | /usr/bin/bc` -eq 1 ];then break;fi
    retain="$retain $r"
done

numretain=`echo $retain|wc -w`
run mirtk padding $f $sdir/corrections/$base-lccs.nii.gz $outf $numretain $retain 0 -invert

rm $sdir/corrections/$base-lccs.nii.gz $sdir/corrections/$base-lccs-comps
fi
