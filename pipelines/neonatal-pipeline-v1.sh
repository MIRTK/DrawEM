#! /bin/bash
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

usage()
{
  base=$(basename "$0")
  echo "usage: $base subjectID [options]
This script runs the neonatal segmentation pipeline of Draw-EM.
Required Options:
  -a / -age  <number>           Subject age in weeks. This is used to select the appropriate template for the initial registration. 
			        If the age is <28w or >44w, it will be set to 28w or 44w respectively.

  Note: The age of the subject can be alternatively provided by a file named 'ages.csv' inside the data directory.
        The format of the file must be a space-delimited file with lines: 'subjectID age'

Additional options:
  -c / -cleanup  <0/1>          Whether cleanup of temporary files is required (default: 1)
  -d / -data-dir  <directory>   The directory used to run the script. 
			        It must contain a T2 folder and the T2/subjectID.nii.gz image of the subject.
  -p / -save-posteriors  <0/1>  Whether the structures' posteriors are required (default: 0)
  -t / -threads  <number>       Number of threads (CPU cores) allowed for the registration to run in parallel (default: 1)
  -v / -verbose  <0/1>          Whether the script progress is reported (default: 1)
  -h / -help / --help           Print usage.
"
  exit;
}



if [ -n "$DRAWEMDIR" ]; then
  [ -d "$DRAWEMDIR" ] || { echo "DRAWEMDIR environment variable invalid!" 1>&2; exit 1; }
else
  export DRAWEMDIR="$(cd "$(dirname "$BASH_SOURCE")"/.. && pwd)"
fi

[ $# -ge 1 ] || { usage; }
subj=$1

case "$subj" in
    -h|-help|--help) usage; ;;
esac

age=""
cleanup=1 # whether to delete temporary files once done
datadir=`pwd`
posteriors=0   # whether to output posterior probability maps
threads=1
verbose=1

atlasname=non-rigid-v2
version="v1"

while [ $# -gt 0 ]; do
  case "$2" in
    -a|-age)  shift; age=$2; ;;
    -c|-cleanup)  shift; cleanup=$2; ;;
    -d|-data-dir)  shift; datadir=$2; ;;
    -p|-save-posteriors) shift; posteriors=$2; ;;
    -t|-threads)  shift; threads=$2; ;; 
    -v|-verbose)  shift; verbose=$2; ;; 
    -h|-help|--help) usage; ;;
    -*) echo "$0: Unrecognized option $1" >&2; usage; ;;
     *) break ;;
  esac
  shift
done


cd "$datadir" || { echo "$0: Directory $datadir does not exist" >&2;exit 1; }

[ -n "$age" -o ! -f ages.csv ] || { age=`cat ages.csv | tr '\t' ' ' | grep "^$subj " | cut -d' ' -f2`; }
[ -n "$age" ] || { echo "Subject age must be provided as argument or from the $datadir/ages.csv file" >&2; exit 1; }
age=`printf "%.*f\n" 0 $age` #round
[ $age -gt 28 ] || { age=28; }
[ $age -lt 44 ] || { age=44; }


[ $verbose -le 0 ] || { echo "DrawEM multi atlas $version
Subject:    $subj 
Age:        $age
Directory:  $datadir 
Posteriors: $posteriors 
Cleanup:    $cleanup 
Threads:    $threads

$BASH_SOURCE $@
----------------------------"; }

mkdir -p logs || exit 1

run()
{
  [ $verbose -le 0 ] || echo -n "$@..."
  "$DRAWEMDIR/scripts/$version/$@" 1>>logs/$subj 2>>logs/$subj-err
  if [ $? -eq 0 ]; then
    [ $verbose -le 0 ] || echo " done"
  else
    [ $verbose -le 0 ] || echo " failed: see log file logs/$subj-err for details"
    exit 1
  fi
}

rm -f logs/$subj logs/$subj-err
run preprocess.sh        $subj $age
# phase 1 tissue segmentation
run tissue-priors.sh     $subj $age $atlasname $threads
# registration using gm posterior + image
run register-multi-atlas-using-gm-posteriors.sh $subj $age $threads
# structural segmentation
run labels-multi-atlas.sh   $subj
run segmentation.sh      $subj
# post-processing
run separate-hemispheres.sh  $subj
run correct-segmentation.sh  $subj
run postprocess.sh       $subj

# if probability maps are required
[ "$posteriors" == "0" -o "$posteriors" == "no" -o "$posteriors" == "false" ] || run postprocess-pmaps.sh $subj

# cleanup
if [ "$cleanup" == "1" -o "$cleanup" == "yes" -o "$cleanup" == "true" ] && [ -f "segmentations/${subj}_labels.nii.gz" ];then
  run clear-data.sh $subj
  rm -f logs/$subj logs/$subj-err
  rmdir logs 2> /dev/null # may fail if other log files exist
fi

exit 0
