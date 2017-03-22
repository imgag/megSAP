#!/bin/bash

#help
if [ $# -lt 1 ]
then
  echo "Usage:"
  echo "1) go to project folder e.g. /mnt/share/projects/gs/IND/:"
  echo "2) call with the following arguments:"
  echo "     normal DNA ID, e.g. GS090022_01 followed by"
  echo "     multiple tumor DNA sample, e.g. GS090020_01"
  echo "   Notice: All samples have to be within the NGSD "
  exit
fi

#paths
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SCRIPTPATH=`pwd -P`
NO=$1

#project path
PROJECTPATH=`pwd -P`

#bash output
COMMAND_STATUS=/mnt/users/all/http/SampleStatus/data/

#map reference
php $DIR/../src/Pipelines/analyze.php -folder Sample_$NO -name $NO -steps ma --log Sample_$NO/analyze_$(date +%Y%m%d%H%M%S).log
status=$?
if [ $status -ne 0 ]; then
	echo "Error calculating reference sample."
	exit 1
fi

for TU in "${@:2}"; do
	#create output folder
	OUT=Somatic_$TU-$NO
	if [ ! -d "$OUT" ]; then
	  mkdir $OUT
	fi
	chmod 775 $OUT

	#perform analysis
	COMMAND="php $DIR/../src/Pipelines/somatic_dna.php -p_folder . -t_id $TU -n_id $NO -o_folder $OUT -smn --log $OUT/somatic_multitum_$(date +%Y%m%d%H%M%S).log"
#	only DNA, no need to restrict to NGSlong
#	JOB=$(qsub -V -b y -wd $SCRIPTPATH -m n -M christopher.schroeder@med.uni-tuebingen.de -q default_srv016,default_srv017,default_srv018 -e $COMMAND_STATUS -o $COMMAND_STATUS $COMMAND)
#	echo -e "$(date +%Y%m%d%H%M%S) \t $COMMAND \n\t -> $JOB" >> $COMMAND_STATUS/commands.txt
	$COMMAND
done