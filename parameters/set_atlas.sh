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

tissue_atlas=$1
atlas=$2

. $DRAWEMDIR/parameters/configuration.sh

atlas_exists=`echo " $AVAILABLE_ATLASES "|grep " $atlas "`
if [ "$atlas_exists" == "" ];then
    echo "Unknown atlas: $atlas" >&2;
    exit;
fi

atlas_exists=`echo " $AVAILABLE_TISSUE_ATLASES "|grep " $tissue_atlas "`
if [ "$atlas_exists" == "" ];then
    echo "Unknown atlas: $tissue_atlas" >&2;
    exit;
fi
# load configuration
export ATLAS_NAME=$atlas
export TISSUE_ATLAS_NAME=$tissue_atlas
. $DRAWEMDIR/parameters/$atlas/configuration.sh
. $DRAWEMDIR/parameters/$tissue_atlas/configuration.sh
