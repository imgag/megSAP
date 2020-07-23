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

if [ ! -d "Sample_${args[0]}" ]; then
	echo "ERROR: Could not find folder Sample_${args[0]}"
	exit
fi
cd Sample_${args[0]}

#bash output
COMMAND_STATUS=/mnt/users/all/http/SampleStatus/data/

#perform analysis
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
COMMAND="php $DIR/../src/Pipelines/analyze.php -folder . -name ${args[0]} ${args[@]:1} --log analyze_$(date +%Y%m%d%H%M%S).log"

#project path
PROJECTPATH=`pwd -P`

if [[ "$noqueue" == true ]]
then
	echo $COMMAND
	$COMMAND
else
	JOB=$(qsub -V -b y -wd $PROJECTPATH -m n -M christopher.schroeder@med.uni-tuebingen.de -q default_srv016,default_srv017,default_srv018 -e $COMMAND_STATUS -o $COMMAND_STATUS $COMMAND)
	echo -e "$(date +%Y%m%d%H%M%S) \t $COMMAND \n\t -> $JOB" >> $COMMAND_STATUS/commands.txt
	echo $JOB
fi
