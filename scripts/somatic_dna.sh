#!/bin/bash

#help
if [ $# -lt 1 ]
then
  echo "Usage:"
  echo "1) go to project folder e.g. /mnt/share/projects/gs/IND/:"
  echo "2) call with the following arguments:"
  echo "     tumor DNA sample name, e.g. GS090020_01"
  echo "     normal DNA sample name, e.g. GS090022_01 (if none available use na or skip this and the following parts)"
  echo "   optional arguments are:"
  echo "   Notice: all other parameters are passed on to the analyze php script, e.g. "
  echo "     -backup"
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
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
COMMAND="php $DIR/../src/Pipelines/somatic_dna.php -p_folder . -t_id ${args[0]} -n_id ${args[1]} -o_folder $OUT --log $OUT/somatic_dna_$(date +%Y%m%d%H%M%S).log ${args[@]:2}"
if [[ "$noqueue" == true ]]
then
	echo $COMMAND
	$COMMAND
else
	#only DNA no need to restrict to NGSlong
	JOB=$(qsub -V -b y -wd $PROJECTPATH -m n -M christopher.schroeder@med.uni-tuebingen.de -q NGSlong,srv016_long,srv018long -e $COMMAND_STATUS -o $COMMAND_STATUS $COMMAND)
	echo -e "$(date +%Y%m%d%H%M%S) \t $COMMAND \n\t -> $JOB" >> $COMMAND_STATUS/commands.txt
	echo $COMMAND
fi