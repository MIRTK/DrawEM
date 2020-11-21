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

# Available atlases from Draw-EM
export AVAILABLE_ATLASES="ALBERT MCRIB"
export AVAILABLE_TISSUE_ATLASES="neonatal"

export SUPER_GM_LABEL=10000
export SUPER_WM_LABEL=20000

# log function
run()
{
  echo "$@"
  "$@"
  if [ ! $? -eq 0 ]; then
    echo "$@ : failed"
    exit 1
  fi
}

# make run function global
typeset -fx run
