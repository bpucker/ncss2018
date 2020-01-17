### Katharina Frey ###
### katharina.frey@uni-bielefeld.de ###
### v0.1 ###

import os
from argparse import ArgumentParser
import glob
import sys

__usage__ = """
python masterscript_ncss.py
--script_dir <FULL_PATH_TO_FOLDER_WITH_SCRIPTS>
--input_dir <FULL_PATH_TO_FOLDER_WITH_GENOME_AND_ANNOTATION_FILES>
	WARNING: file names need to match species names (example: Arabidopsis_thaliana.fna)
--python_dir <FULL_PATH_TO_PYTHON_INSTALLATION>
	default = python
--output_dir <FULL_PATH_TO_OUPUT_FOLDER>
--species_file <FULL_PATH_TO_FILE_CONTAINING_A_LIST_OF_SPECIES_NAMES>
	example: Arabisopsis_thaliana (separated with \\n from next species)
--bam_dir <FULL_PATH_TO_FOLDER_WITH_BAM_FILES>
	WARNING: file names need to match genome file names (example: Arabidopsis_thaliana.bam)
					"""

input_parameters = ArgumentParser()
input_parameters.add_argument("--script_dir", dest="script_directory", required=True)
input_parameters.add_argument("--input_dir", dest="input_directory", required=True)
input_parameters.add_argument("--python_dir", dest="python_directory", default="python")
input_parameters.add_argument("--output_dir", dest="output_directory", required=True)
input_parameters.add_argument("--species_file", dest="species_f", required=True)
input_parameters.add_argument("--bam_dir", dest="bam_directory", required=True)

if "--help" in sys.argv or "-h" in sys.argv:
    print(__usage__)
    sys.exit(1)

args = input_parameters.parse_args()

if args.script_directory is None:
	print("\n'--script_dir' was not set'")
	print(__usage__)
elif args.input_directory is None:
	print("\n'--input_dir' was not set'")
	print(__usage__)
elif args.python_directory is None:
	print("\n'--python_dir' was not set'")
	print(__usage__)	
elif args.output_directory is None:
	print("\n'--output_dir' was not set'")
	print(__usage__)
elif args.species_f is None:
	print("\n'--species_file' was not set'")
	print(__usage__)	
elif args.bam_directory is None:
	print("\n'--bam_dir' was not set'")
	print(__usage__)
else:
	#NCBI_splice_site_check.py
	files = glob.glob(args.input_directory + "*.fna")
	for genomefile in files:
		species = genomefile[(len(args.input_directory)):-4]
		if os.path.isfile(args.input_directory + species + ".g2t") == False:
			os.popen(args.python_directory + " " + args.script_directory + "NCBI_splice_site_check.py --in " + args.input_directory + " > " + args.output_directory + "NCBI_splice_site_check.log")
		else:
			print(args.input_directory + species + ".g2t already exists. Continue with further steps...")
	
	#splice_site_divergence_check.py
	if os.path.isdir(args.output_directory + "ss_diverg_check/") == False:
		os.popen(args.python_directory + " " + args.script_directory + "splice_site_divergence_check.py --data_dir " + args.input_directory + " --output_dir " + args.output_directory + "ss_diverg_check/ --species_file " + args.species_f  + " > " + args.output_directory + "splice_site_divergence_check.log")
	else:
		print(args.output_directory + "ss_diverg_check/ already exists. Continue with further steps...")
	
	#construct_RNA_seq_coverage_file.py & check_ncss.py
	if os.path.isdir(args.output_directory + "coverage_files/") == False:
		os.makedirs(args.output_directory + "coverage_files/")
	files = glob.glob(args.bam_directory + "*.bam")
	for bamfile in files:
		species = bamfile[(len(args.bam_directory)):-4]
		if os.path.isfile(args.output_directory + "coverage_files/" + species + ".cov") == False:
			os.popen(args.python_directory + " " + args.script_directory + "construct_RNA_seq_coverage_file.py --in " + bamfile + " --out " + args.output_directory + "coverage_files/" + species + ".cov --bam_is_sorted" + " > " + args.output_directory + "construct_RNA_seq_coverage_file.log")
		else:
			print(args.output_directory + species + ".cov already exists. Continue with further steps...")
		if os.path.isdir(args.output_directory + "check_ncss/" +species + "/") == False:
			os.popen(args.python_directory + " " + args.script_directory + "check_ncss.py --cov " + args.output_directory + "coverage_files/" + species + ".cov --doc " + args.input_directory + species + ".txt --out " + args.output_directory + "check_ncss/" +species + "/" + " > " + args.output_directory + "check_ncss.log")
		else:
			print(args.output_directory + "check_ncss/" + species + "/ already exists. Continue with further steps...")
	
	#get_RNA_seq_amount.py
	if os.path.isfile(args.output_directory + "RNA_seq_amount.txt") == False:
		os.popen(args.python_directory + " " + args.script_directory + "get_RNA_seq_amount.py --cov_dir " + args.output_directory + "coverage_files/ --out " + args.output_directory + "RNA_seq_amount.txt" + " > " + args.output_directory + "get_RNA_seq_amount.log")
	else:
		print(args.output_directory + "RNA_seq_amount.txt already exists. Continue with further steps...")
	
	#analyze_supported_splice_sites.py
	if os.path.isdir(args.output_directory + "supported_splice_sites/") == False:
		os.popen(args.python_directory + " " + args.script_directory + "analyze_supported_splice_sites.py --data " + args.input_directory + " --sssd " + args.output_directory + "check_ncss/ --cov_rep " + args.output_directory + "RNA_seq_amount.txt --out " + args.output_directory + "supported_splice_sites/" + " > " + args.output_directory + "analyze_supported_splice_sites.log")
	else:
		print(args.output_directory + "supported_splice_sites/ already exists. Continue with further steps...")
	
	#intron_size_analysis.py	
	if os.path.isdir(args.output_directory + "intron_size_analysis/") == False:
		os.popen(args.python_directory + " " + args.script_directory + "intron_size_analysis.py --in " + args.input_directory + " --out " + args.output_directory + "intron_size_analysis/" + " > " + args.output_directory + "intron_size_analysis.log")
	else:
		print(args.output_directory + "intron_size_analysis/ already exists. Continue with further steps...")

print("---finished---")
