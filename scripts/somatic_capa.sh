#!/bin/bash

#help
if [ $# -lt 1 ]
then
  echo "Usage:"
  echo "1) go to project folder e.g. /mnt/share/projects/gs/IND/:"
  echo "2) call with the following arguments:"
  echo "     tumor sample e.g. GS090020"
  echo "     optional: reference sample e.g. GS090021"
  echo "     optional if samples should not be queued: NOQUEUE"
  exit
fi

#extract all arguments
args=( "$@" )

#check if NOQUEUE is available
noqueue=false
for (( i=0; i < ${#args[@]}; ++i))
do
	if [[ ${args[$i]} == "NOQUEUE" ]]
	then
		noqueue=true
		unset args[$i]
	fi
done
args=( "${args[@]}" )

#project path
PROJECTPATH=`pwd -P`

#bash output
COMMAND_STATUS=/mnt/users/all/http/SampleStatus/data/

#create output folder
OUT=Sample_${args[0]}
if [[ "${args[1]}" != "" && "${args[1]}" != "na" ]]
then
	OUT=Somatic_${args[0]}-${args[1]}
fi
if [ ! -d "$OUT" ]; then
  mkdir -m 775 $OUT
fi

#perform analysis
if [[ "${args[1]}" != "" && "${args[1]}" != "-" ]]; then
    EXTRA="$EXTRA -n_id ${args[1]}"
fi
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
COMMAND="php $DIR/../src/Pipelines/somatic_capa.php -p_folder . -t_id ${args[0]} -o_folder $OUT $EXTRA ${args[@]:2} --log $OUT/somatic_capa_$(date +%Y%m%d%H%M%S).log"
if [[ "$noqueue" == true ]]
then
	echo $COMMAND
	$COMMAND
else
	#only DNA no need to restrict to NGSlong
	JOB=$(qsub -V -b y -wd $PROJECTPATH -m n -M christopher.schroeder@med.uni-tuebingen.de -q default_srv016,default_srv017,default_srv018 -e $COMMAND_STATUS -o $COMMAND_STATUS $COMMAND)
	echo -e "$(date +%Y%m%d%H%M%S) \t $COMMAND \n\t -> $JOB" >> $COMMAND_STATUS/commands.txt
	echo $JOB
fi
