#help
if [ $# -lt 2 ]
then
  echo "Usage:"
  echo " vc_freebayes.sh [chr:start-end] [bam] <options>"
  echo "Options:"
  echo " --min-alternate-fraction [0.1]"
  echo " --min-mapping-quality [1]"
  echo " --min-base-quality [10]"
  echo " --min-alternate-qsum [90]"
  exit
fi

/mnt/share/opt/freebayes-1.1.0/bin/freebayes -f /mnt/share/data/genomes/GRCh37.fa -r $1 -b $2  ${@:3}
