#!/bin/bash

#help
if [ $# -lt 1 ]
then
  echo "Usage:"
  echo "1) go to project folder e.g. /mnt/share/projects/gs/IND/:"
  echo "2) call with the following arguments:"
  echo "     tumor DNA sample name, e.g. GS090020_01"
  echo "     relapse DNA sample name, e.g. GS090022_01"
  echo "     normal DNA sample name, e.g. GS090022_01"
  echo "   Notice: all other parameters are passed on to the analyze php script, e.g. "
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
OUT=Somatic_${args[0]}-${args[1]}-${args[2]}
if [ ! -d "$OUT" ]; then
  mkdir $OUT
fi
chmod 775 $OUT

#perform analysis
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
COMMAND="php $DIR/../src/Pipelines/somatic_triplett.php -p_folder . -t_dna_id ${args[0]} -r_dna_id ${args[1]} -n_dna_id ${args[2]} -o_folder $OUT --log $OUT/somatic_trio_$(date +%Y%m%d%H%M%S).log  ${args[@]:3}"
if [[ "$noqueue" == true ]]
then
	$COMMAND
else
	#no RNA data therefore no restrict to NGSlong
	JOB=$(qsub -V -b y -wd $PROJECTPATH -m n -M christopher.schroeder@med.uni-tuebingen.de -q default_srv016,default_srv017,default_srv018 -e $COMMAND_STATUS -o $COMMAND_STATUS $COMMAND)	
	echo -e "$(date +%Y%m%d%H%M%S) \t $COMMAND \n\t -> $JOB" >> $COMMAND_STATUS/commands.txt
	echo $COMMAND
fi