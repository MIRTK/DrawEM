#! /bin/bash

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
