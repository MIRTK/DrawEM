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


[ $# -eq 1 ] || { echo "usage: $(basename "$0") <image>" 1>&2; exit 1; }

image=$1

image_base=`echo $image | sed -e "s:.nii.gz$::g" | sed -e "s:.nii$::g"`
for thr in 0.1 0.2 0.3 0.4 0.5 0.6 0.7;do
    bet $image ${image_base}_thr_${thr}.nii.gz -R -f $thr -n -m
    rm ${image_base}_thr_${thr}.nii.gz
done
