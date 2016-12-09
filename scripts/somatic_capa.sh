#!/bin/bash

#help
if [ $# -lt 1 ]
then
  echo "Usage:"
  echo "1) go to project folder e.g. /mnt/share/projects/gs/IND/:"
  echo "2) call with the following arguments:"
  echo "     tumor sample e.g. GS090020"
  echo "     optional: reference sample e.g. GS090021"
  echo "     optional: system"
  exit
fi


#create output folder
if [[ "$2" != "" && "$2" != "-" ]]; then
	OUT=Somatic_$1-$2
else
	OUT=Somatic_$1
fi
if [ ! -d "$OUT" ]; then
  mkdir $OUT
fi

#add reference sample
if [[ "$2" != "" && "$2" != "-" ]]; then
    EXTRA="$EXTRA -nn $2"
fi

#perform analysis
if [ ! -d "$OUT" ]; then
	echo "Exit analyze.sh: Could not find folder $OUT."
	exit
fi
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "$DIR/../src/Pipelines/somatic_capa.php -pf . -tn $1 -o_folder $OUT $EXTRA ${@:3} --log $OUT/capa_$(date +%Y%m%d%H%M%S).log"
php $DIR/../src/Pipelines/somatic_capa.php -pf . -tn $1 -o_folder $OUT $EXTRA ${@:3} --log $OUT/capa_$(date +%Y%m%d%H%M%S).log