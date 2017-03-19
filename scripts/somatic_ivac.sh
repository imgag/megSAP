#!/bin/bash
#help
if [ $# -lt 1 ]
then
  echo "Usage:"
  echo "if job should not be submitted to queue use 'NOQUEUE' as argument"
  echo "1) go to project folder e.g. /mnt/share/projects/gs/IND/:"
  echo "2) call with the following arguments:"
  echo "     tumor DNA sample name, e.g. GS090020_01"
  echo "     reference DNA sample name, e.g. GS090022_01"
  echo "     tumor RNA sample name, e.g. GS090020_01"
  echo "     reference RNA sample name, e.g. GS090022_01"
  echo "   Notice: all other parameters are passed on to the analyze php script, e.g. "
  echo "     systems"
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

#create output folder
OUT=Somatic_${args[0]}-${args[1]}

#project path
PROJECTPATH=`pwd -P`

#bash output
COMMAND_STATUS=/mnt/users/all/http/SampleStatus/data/

#look up RNA data, checke folders (if * in paths cannot be resolved strange error messages appear)
RNA="";
if [[ "${args[2]}" != "" && "${args[2]}" != "-"  && "${args[2]}" != "na" ]]
then
	MESSAGE=""
	T_RNA_FOLDER="$PROJECTPATH/Sample_${args[2]}"
	if [ ! -d "$T_RNA_FOLDER" ]; then
		MESSAGE="$MESSAGE Could not find tumor RNA folder $T_RNA_FOLDER. "
	else
		MESSAGE=""
	fi
	RNA="$RNA -t_rna_fo $T_RNA_FOLDER -t_rna_id ${args[2]}"

	#normal RNA
	if [[ "${args[3]}" != "" && "${args[3]}" != "-" && "${args[3]}" != "na" ]]
	then
		N_RNA_FOLDER="$PROJECTPATH/Sample_${args[3]}"
		if [ ! -d "$N_RNA_FOLDER" ]; then
			MESSAGE="$MESSAGE Could not find normal RNA folder $N_RNA_FOLDER. "
		else
			MESSAGE=""
		fi
		RNA="$RNA -n_rna_fo $N_RNA_FOLDER -n_rna_id ${args[3]}"
	fi

	if [[ "$MESSAGE" != "" ]]
	then
		echo "$MESSAGE"
		exit
	fi
fi

if [ ! -d "$OUT" ]; then
  mkdir $OUT
fi

#perform analysis
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
COMMAND="php $DIR/../src/Pipelines/somatic_ivac.php -p_folder . -t_dna_id ${args[0]} -n_dna_id ${args[1]} -o_folder $OUT $RNA --log $OUT/somatic_ivac_$(date +%Y%m%d%H%M%S).log ${args[@]:4}"
QUEUE="-q default_srv016,default_srv017,default_srv018";

if [[ "$noqueue" == true ]]
then
	$COMMAND
else
	JOB=$(qsub -V -b y -wd $PROJECTPATH -m n -M christopher.schroeder@med.uni-tuebingen.de $QUEUE -e $COMMAND_STATUS -o $COMMAND_STATUS $COMMAND)	
	echo -e "$(date +%Y%m%d%H%M%S) \t $COMMAND \n\t -> $JOB" >> $COMMAND_STATUS/commands.txt
	echo $COMMAND
fi
