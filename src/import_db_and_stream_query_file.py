# use the metagenome conda env
# stream query fasta or fastq files to the AA k-mer object

### import modules
import marisa_trie as mt
from hydra import WritingBloomFilter, ReadingBloomFilter
import numpy as np
import pandas as pd
import os
import time
from argparse import ArgumentTypeError
import argparse
import pickle
import sys


def local_tests():
	os.chdir("output")
	src_folder = '/Users/shaopeng/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Github/CMash_Rust_implementation/src'
	ref_file = 'fmh_scale_10_KO_record.pkl'
	bf_file = 'prefilter_bf_KO_record.db'
	out_label = 'test'
	query_file = '/Users/shaopeng/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Github/CMash_Rust_implementation/Amino_acid_index/full_path_3_faa_files.txt'
	query_type = 'fasta'

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Build AA FracMinHash sketches from a list of input ",
	                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-r', '--ref_database', type=str, help="K-mer database to scan for")
	parser.add_argument('-b', '--bf_prefilter', type=str, help="A Bloomfilter for prefilter")
	parser.add_argument('-l', '--label', type=str, help="Keyword for out file", default='KO_record')
	parser.add_argument('-i', '--input_file', type=str, help="A file containing paths of input fasta/fastq files")
	parser.add_argument('-p', '--input_type', type=str, help="File type: fasta/fastq")
	
	### load object functions
	src_folder = os.path.realpath(os.path.dirname(__file__))
	sys.path.insert(1, src_folder)
	print("Loading object functions form %s" % src_folder)
	from build_merged_TST_from_ref_AA_w_KOs import *
	
	
	### read parameters
	args = parser.parse_args()
	# kmer database
	ref_file = args.ref_database
	
	temp_start = time.time()
	if not os.path.exists(ref_file):
		raise ArgumentTypeError("K-mer database %s does not exist" % (ref_file))
	else:
		with open(ref_file, "rb") as fp:
			TST_obj = pickle.load(fp)
			TST_obj.print_stats()
	
	# bf prefilter
	bf_file = args.bf_prefilter
	if not os.path.exists(bf_file):
		raise ArgumentTypeError("Bloomfilter %s does not exist" % (bf_file))
	else:
		TST_obj.prefilter_bf = ReadingBloomFilter(bf_file)
	
	temp_end = time.time()
	print("Total loading time is %s seconds\n\n\n\n\n" % round((temp_end - temp_start), 3))
	
	
	
	# stream input files
	query_file = args.input_file
	query_type = args.input_type
	out_label = args.label
	print("The input file is: %s\n" % query_file)
	print("The input type is: %s\n" % query_type)
	print("Output label is: %s\n\n" % out_label)
	
	
	input_list = check_files(query_file)
	if len(input_list) == 0:
		raise Exception("Didn't find any file from input")
	
	for seq_file in input_list:
		# obj contains time counting information
		basename_file = os.path.basename(seq_file)
		out_name = "CI_scan_" + out_label + "_" + basename_file + ".csv"
		TST_obj.stream_a_query_file(seq_file, type=query_type, out_name=out_name)

	
	
	print("Pipe done")
	
	
	




