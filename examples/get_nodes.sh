#!/bin/bash

# Get number of nodes of a given network
# Author: Fabricio Murai
# Creation date: 09/11/2015
# Last modified: 09/11/2015 11:58 AM EST

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"


if [[ $# -lt 1 ]]; then
    echo "Dataset path not provided. Terminating."
    exit 2
fi

INFOFILE=$DIR/basic_stats.txt
DATASET=`basename $1`
LINE=`grep $DATASET $INFOFILE`

if [[ $? -ne 0 ]]; then
    echo "Dataset $DATASET not found in $INFOFILE."
    exit 1
else
    NODES=`echo $LINE | awk '{print $2}'`
    echo $NODES
    exit 0
fi
