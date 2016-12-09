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

sdir=segmentations-data
[ -n "$sdir" ] || { echo "$BASH_SOURCE: sdir must not be empty!" 1>&2; exit 1; }

rm -f $sdir/MADs/$subj.nii.gz $sdir/MADs/$subj-grad.nii.gz $sdir/MADs/$subj-subspace.nii.gz
rm -f $sdir/atlas-weights/$subj-ALBERT_*.nii.gz
rm -f $sdir/corrections/$subj-gmtochange.nii.gz $sdir/corrections/$subj-ventohwm.nii.gz
rm -f $sdir/cortical/$subj.nii.gz
rm -f $sdir/gm-posteriors/$subj.nii.gz
rm -f $sdir/labels/*/$subj.nii.gz
rm -f $sdir/posteriors/*/$subj.nii.gz
rm -f $sdir/template/*/$subj.nii.gz
rm -f $sdir/tissue-initial-segmentations/$subj.nii.gz
rm -f $sdir/tissue-posteriors/*/$subj.nii.gz
rm -f $sdir/transformations/T2-$subj-ALBERT_*.nii.gz $sdir/transformations/$subj-ALBERT_*.nii.gz $sdir/transformations/tissues-$subj-ALBERT_*.nii.gz
rm -f segmentations/$subj-em.nii.gz segmentations/$subj-initial.nii.gz
rm -f logs/$subj logs/$subj-err logs/$subj-em logs/$subj-em-err logs/$subj-tissue-em logs/$subj-tissue-em-err

#may fail if other files present, which is ok
# rmdir $sdir/atlas-weights                    2> /dev/null
# rmdir $sdir/corrections                      2> /dev/null
# rmdir $sdir/cortical                         2> /dev/null
# rmdir $sdir/gm-posteriors                    2> /dev/null
# rmdir $sdir/labels/*                         2> /dev/null
# rmdir $sdir/labels                           2> /dev/null
# rmdir $sdir/MADs                             2> /dev/null
# rmdir $sdir/posteriors/*                     2> /dev/null
# rmdir $sdir/posteriors                       2> /dev/null
# rmdir $sdir/template/*                       2> /dev/null
# rmdir $sdir/template                         2> /dev/null
# rmdir $sdir/tissue-initial-segmentations     2> /dev/null
# rmdir $sdir/tissue-posteriors/*              2> /dev/null
# rmdir $sdir/tissue-posteriors                2> /dev/null
# rmdir $sdir/transformations                  2> /dev/null
# rmdir $sdir                                  2> /dev/null
# rmdir logs                                   2> /dev/null

exit 0











