#! /bin/bash

atlasname=$1

. $DRAWEMDIR/parameters/configuration.sh

atlas_exists=`echo " $AVAILABLE_ATLASES "|grep " $atlasname "`
if [ "$atlas_exists" == "" ];then
	echo "Unknown atlas: $atlasname" >&2;
    exit;
fi
# load configuration
export ATLAS_NAME=$atlasname
. $DRAWEMDIR/parameters/$ATLAS_NAME/configuration.sh
