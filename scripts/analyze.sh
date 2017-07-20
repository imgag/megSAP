#!/bin/bash

#help
if [ $# -lt 1 ]
then
  echo "Usage:"
  echo "1) go to project folder e.g. /mnt/share/projects/gs/IND/:"
  echo "2) call with the following arguments:"
  echo "     sample name, e.g. GS090020"
  echo "     processing system, e.g. SeqCapEZv2"
  echo "   Notice: all other parameters are passed on to the analyze php script, e.g. "
  echo "     -backup"
  exit
fi

if [ ! -d "Sample_$1" ]; then
	echo "ERROR: Could not find folder Sample_$1"
	exit
fi
cd Sample_$1

#perform analysis
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
php $DIR/../src/Pipelines/analyze.php -folder . -name $1 ${@:2} --log analyze_$(date +%Y%m%d%H%M%S).log