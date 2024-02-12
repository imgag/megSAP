

import argparse
import sys
import os
import tempfile
import shutil

from distutils.dir_util import copy_tree
from distutils.file_util import copy_file

# set max threads to be used to remove warning 
os.environ["NUMEXPR_MAX_THREADS"] = "10"

from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerExtractor import sigpro as sig


def main(args):
    
    if args.installGenome:
        genInstall.install(args.ref) # only once per genome
        return
    
    test_input(args)

    # create temp folder that will be deleted after use
    if args.mode == "snv":
        calculate_snv_signatures(args)

    if args.mode == "cnv":
        calculate_cnv_signatures(args)


def calculate_cnv_signatures(args):
    with tempfile.TemporaryDirectory() as tmpdir:
    # tmpdir = "/mnt/storage3/users/ahott1a1/tmp_data/megsap_test/"
    # if True:
    
        os.chdir(tmpdir)  # in cwd the result folder will be created so change to a temp folder
        result_folder = "results_cnv/"

        tmpfile = cnv_prepare_clincnv(tmpdir, args.input)
        
        if tmpfile == "":
            print("Warning: No CNVs called! - Skipping CNV signature calculation.")
            return
        
        sig.sigProfilerExtractor("seg:FACETS", result_folder, tmpfile, reference_genome=args.ref,
                                 minimum_signatures=args.minSig, maximum_signatures=args.maxSig,
                                 nmf_replicates=args.nmfRep, cpu=args.threads, seeds=args.seeds)


        copy_cnv_result_files(args, result_folder)

def cnv_prepare_clincnv(tmpdir, cnvFile):
    lines = []
    with open(cnvFile) as file:
        for line in file:
            if line.startswith("#"):
                continue

            line = line.strip()
            if len(line) == 0:
                continue
            parts = line.split("\t")
            new_line_parts = [parts[0], parts[1], parts[2], parts[5], parts[8], os.path.basename(cnvFile).replace("_clincnv.tsv", ""), "0"]
            lines.append("\t".join(new_line_parts))
    
    if len(lines) == 0:
        return ""

    tmpfile = tmpdir + "/cnv_input.seg"
    with open(tmpfile, "w") as file:
        # "chr start end tcn.em lcn.em sample id
        file.write("chr\tstart\tend\ttcn.em\tlcn.em\tsample\tid")
        for line in lines:
            file.write(line + "\n")

    return tmpfile


def copy_cnv_result_files(args, result_folder):
    if args.fullOutput:
        copy_tree(result_folder, args.outFolder)

    else:
        file_names = [
                      "/JOB_METADATA.txt",
                      "/Seeds.txt",
                      "/CNV48/Suggested_Solution/COSMIC_CNV48_Decomposed_Solution/CNV_Decomposition_Plots.pdf",
                      "/CNV48/Suggested_Solution/COSMIC_CNV48_Decomposed_Solution/De_Novo_map_to_COSMIC_CNV48.csv"
                      ]
        for file in file_names:
            copy_file(result_folder + file, args.outFolder)


def calculate_snv_signatures(args):
    with tempfile.TemporaryDirectory() as tmpdir:
    # tmpdir = "/mnt/storage3/users/ahott1a1/tmp_data/megsap_test/"
    # if True:
        os.chdir(tmpdir)  # in cwd the result folder will be created so change to a temp folder
        print("tmp working dir: " +  tmpdir)
           
        prepare_input_files(args, tmpdir)
        result_folder = "results_snv/"

        sig.sigProfilerExtractor("vcf", result_folder, args.input, reference_genome=args.ref,
                                 minimum_signatures=args.minSig, maximum_signatures=args.maxSig,
                                 nmf_replicates=args.nmfRep, cpu=args.threads, seeds=args.seeds)


        copy_snv_result_files(args, result_folder)


def copy_snv_result_files(args, result_folder):
    if args.fullOutput:
        copy_tree(result_folder, args.outFolder)

    else:
        file_names = [
                      "/JOB_METADATA.txt",
                      "/Seeds.txt",
                      "/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/SBS96_Decomposition_Plots.pdf",
                      "/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/De_Novo_map_to_COSMIC_SBS96.csv",
                      "/ID83/Suggested_Solution/COSMIC_ID83_Decomposed_Solution/ID83_Decomposition_Plots.pdf",
                      "/ID83/Suggested_Solution/COSMIC_ID83_Decomposed_Solution/De_Novo_map_to_COSMIC_ID83.csv",
                      "/DBS78/Suggested_Solution/COSMIC_DBS78_Decomposed_Solution/De_Novo_map_to_COSMIC_DBS78.csv",
                      "/DBS78/Suggested_Solution/COSMIC_DBS78_Decomposed_Solution/DBS78_Decomposition_Plots.pdf"
                      ]
        for file in file_names:
            if os.path.exists(result_folder + file):
                copy_file(result_folder + file, args.outFolder)

def prepare_input_files(args, tmpdir):
    if os.path.isdir(args.input):
        return

    if os.path.isfile(args.input):
        input_folder = tmpdir + "/vcf_files/"
        os.mkdir(input_folder)
        shutil.copy(args.input, input_folder)
        args.input = input_folder

    return


def test_input(args):
    if not os.path.isfile(args.input) and not os.path.isdir(args.input):
        raise ValueError("Input: '" + args.input + "' is neither a file nor a directory.")

    if os.path.isdir(args.input) and args.mode == "cnv":
        raise ValueError("The input cannot be a directory for the cnv signature calculation.")

    args.input = os.path.abspath(args.input)

    if not os.path.isdir(args.outFolder):
        os.mkdir(args.outFolder)

    args.outFolder = os.path.abspath(args.outFolder)

    if args.minSig < 1 or args.maxSig < 1:
        raise ValueError("Neither the minimum nor the maximum number of signatures ia allowed to be smaller than one.")

    if args.nmfRep < 1:
        raise ValueError("Number of nmf replicates cannot be lower than 1.")

    if args.maxSig < args.minSig:
        raise ValueError("The minimum number of signatures has to be less or equal to the maximum number of signatures.")
    
    if args.seeds != "random" and not os.path.isfile(args.seeds):
        raise ValueError("The given file for argument --seeds does not exist.")
    
    if args.seeds != "random":
        args.seeds = os.path.abspath(args.seeds)
    
    
    

if __name__ == '__main__':
    # Read parameters
    parser = argparse.ArgumentParser(description='Extracting genome mutational signatures from a vcf file.')
    parser.add_argument('--in', required=True, dest='input',
                        help='VCF file or folder with multiple VCF files for SNV Signatures'
                             ' or somatic ClinCNV segmentation file for CNV signatures.')
    parser.add_argument('--outFolder', required=True, dest='outFolder', help='Output folder.')

    # optional
    parser.add_argument('--minSignatures', required=False, dest='minSig', type=int, default=1,
                        help='Minimum number of Signatures tested.')
    parser.add_argument('--maxSignatures', required=False, dest='maxSig', type=int, default=10,
                        help='Maximum number of Signatures tested.')
    parser.add_argument('--nmfReplicates', required=False, dest='nmfRep', type=int, default=50,
                        help='Number of nmf replicates.')
    parser.add_argument('--reference', required=False, dest='ref', type=str, default="GRCh38",
                        choices=["GRCh37", "GRCh38", "mm9", "mm10"], help="Reference genome used.")
    parser.add_argument('--mode', required=False, dest='mode', type=str, default="snv",
                        choices=["snv", "cnv"], help="The kind of signatures to be calculated.")
    parser.add_argument('--fullOutput', required=False, dest='fullOutput', action='store_true', help="copies the full result folder to the given outFolder.")
    parser.add_argument('--installGenome', required=False, dest='installGenome', action='store_true', help="Installs the given reference and stops. Only needs to be done once - throws error if already done.")
    parser.add_argument('--threads', required=False, dest='threads', type=int, default=4,
                        help='Number of threads used.')
    parser.add_argument('--seeds', required=False, dest='seeds', default="random",
                        help='VCF file or folder with multiple VCF files for SNV Signatures'
                             ' or somatic ClinCNV segmentation file for CNV signatures.')
    args = ''
    try:
        args = parser.parse_args()
    except IOError as io:
        print(io)
        sys.exit('Error reading parameters.')

    main(args)
