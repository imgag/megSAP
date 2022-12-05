
library("scarHRD")
library("optparse")
 
# segmentation file format:
# SampleID Chromosome Start_position End_position total_cn A_cn B_cn ploidy
 
option_list = list(
  make_option(c("-s", "--seg"), type="character", default=NULL, 
              help="Segmentation file name", metavar="character"),
  make_option(c("-o", "--outDir"), type="character", default=NULL, 
              help="Output file name", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

scar_score(opt$seg ,reference = "grch38", seqz=FALSE, outputdir = opt$outDir)

