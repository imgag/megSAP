#!/bin/bash

set -e -u -o pipefail

# project specific parameters for somatic_dna_rna
project_args="-germline_suffix _adme -germline_target /mnt/projects/research/eMed-HCC/+IKP/ADME_Genes_hg19_eMED_20161206.bed"

# SGE parameters
_status='/mnt/users/all/http/SampleStatus/data'
_queue='default_srv017,default_srv018'

# resolve script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# usage information
usage() {
    cat <<EOT
Usage: $0 [NOQUEUE] tumor-dna-id normal-dna-id [tumor-rna-id] [normal-rna-id] [additional somatic_dna_rna arguments]

Sample folders need to be present in the current working directory.
If NOQUEUE is specified, the pipeline call is not submitted to SGE.
EOT
}

#extract all arguments
args=( "$@" )

# check for NOQUEUE argument
NOQUEUE=false
for (( i=0; i < ${#args[@]}; ++i))
do
	if [[ ${args[$i]} == "NOQUEUE" ]]
	then
		NOQUEUE=true
		unset args[$i]
	fi
done
args=( "${args[@]}" )

# check for at least 2 remaining parameters
if [[ $# -lt 2 ]]; then
    usage
    exit 2
fi

# output directory
OUT_DIR="Somatic_${args[0]}-${args[1]}"
mkdir -p $OUT_DIR

# log file path
LOG="${OUT_DIR}/somatic_dna_rna_$(date +%Y%m%d%H%M%S).log"

# pipeline command
COMMAND="php ${SCRIPT_DIR}/../src/Pipelines/somatic_dna_rna.php -p_folder . -t_dna_id ${args[0]} -n_dna_id ${args[1]} -t_rna_id ${args[2]:-na} -n_rna_id ${args[3]:-na} ${project_args} ${args[@]:4} --log $LOG"

# call
if $NOQUEUE; then
    $COMMAND
else
    JOB=$(qsub -pe smp 1 -V -b y -wd $(realpath $PWD) -m n -M christopher.schroeder@med.uni-tuebingen.de -e $_status -o $_status -q $_queue $COMMAND)
    echo -e "$(date +%Y%m%d%H%M%S) \t $COMMAND \n\t -> $JOB" >> "${_status}/commands.txt"
    echo $COMMAND
fi