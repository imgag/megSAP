#help
if [ $# -lt 2 ]
then
  echo "Usage:"
  echo " vc_freebayes.sh [chr:start-end] [bam] <options>"
  echo "Options:"
  echo " --min-alternate-fraction [0.2]"
  echo " --min-mapping-quality [1]"
  exit
fi

/mnt/share/opt/freebayes-1.2.0/bin/freebayes -b $2 -r $1 -f /mnt/share/data/genomes/GRCh37.fa --min-base-quality 20 --min-alternate-qsum 90 ${@:3}
