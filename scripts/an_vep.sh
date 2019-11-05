#help
if [ $# -lt 2 ]
then
  echo "Usage:"
  echo " an_vep.sh [input] [output.vcf] <options>"
  exit
fi

vep_path=/mnt/share/opt/ensembl-vep-release-98.3/
vep_data_path=/mnt/share/data/dbs/ensembl-vep-98/

PERL5LIB=$vep_path/Bio/:$vep_path/cpan/lib/perl5/:$PERL5LIB
$vep_path/vep -i $1 -o $2 -vcf --no_stats --force_overwrite --offline --cache --dir_cache $vep_data_path/cache/ --fasta /tmp/local_ngs_data/GRCh37.fa --species homo_sapiens --assembly GRCh37
