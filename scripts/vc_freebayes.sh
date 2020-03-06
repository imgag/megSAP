#help
if [ $# -lt 2 ]
then
  echo "Usage:"
  echo " vc_freebayes.sh [chr:start-end] [bam] <options>"
  echo "Default options for germline calling:"
  echo " --min-alternate-fraction 0.1 --min-mapping-quality 1 --min-base-quality 15 --min-alternate-qsum 90"
  exit
fi

/mnt/share/opt/freebayes-1.2.0/bin/freebayes -b $2 -r $1 -f /mnt/share/data/genomes/GRCh37.fa ${@:3}


