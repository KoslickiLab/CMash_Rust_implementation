# use the metagenome conda env

### import modules
import marisa_trie as mt
from hydra import WritingBloomFilter, ReadingBloomFilter
import numpy as np
import pandas as pd
import os
import time
from Bio import SeqIO
import khmer
from argparse import ArgumentTypeError
import argparse
import pickle


### basic functions
def check_files(input_file):
	"""
	Read all file paths from an input file and check if they exist
	:param input_file:
	:return: a sorted list of file paths in the input file
	"""
	print("Reading paths from %s to a list.\n" %input_file)
	out_list = list()
	with open(input_file, 'r') as temp:
		for line in temp.readlines():
			line = line.strip()
			if not os.path.exists(line):
				raise Exception("Input file %s does not exist." % line)
			out_list.append(os.path.abspath(line))
		out_list = sorted(out_list, key=os.path.basename)
		print("There are %d files in the list.\n\n\n\n\n" %len(out_list))
		return (out_list)

def unit_test():
	"""
	Just check if it works as expected
	"""
	with open("ref_file.fasta", 'w') as f:
		f.write(">erl:AOC36_00005|chromosomal replication initiator DnaA|ko:K1\n")
		f.write("ACGTACG")
	
	# store kmers
	test_obj = KO_TST(ksize=4, fmh_scale=1, label="test")  # keep all 4-mers
	test_obj.load_kmer_from_fasta_into_dict("ref_file.fasta")
	assert list(test_obj.kmer_dict.keys()) == ["ACGT", "CGTA", "GTAC", "TACG"]
	
	# calculate CI
	with open("query_file.fasta", 'w') as f:
		f.write(">erl:AOC36adfaadfr DnaA|ko:K2\n")
		f.write("ACGTAMKJGNSJGK")
	# the first 2 4-mers are matches
	test_obj.stream_a_query_file("query_file.fasta", type="fasta", out_name="test_CI.csv")
	# this should give you 0.5
	temp_df = pd.read_csv("test_CI.csv")
	temp_ci = temp_df.loc[0, 'CI']
	assert temp_ci == 0.5
	for file in ["query_file.fasta", "ref_file.fasta", "prefilter_bf_test.db.desc", "prefilter_bf_test.db",
	             "test_CI.csv"]:
		os.remove(file)


### KO_TST object
class KO_TST(object):
	"""
	A TST object with auxiliary arrays for KO. It's tuned for our AA fasta files where KO is the last element in name.
	"""
	def __init__(self, max_prime=9999999999971., ksize=7, rev_comp=False, fmh_scale=10, threads=16, label="ko_aa_obj"):
		"""
        Parameters
        ----------
        max_prime : int
            The maximum prime number to use for AA kmers in KO.
        ksize : int
            The k mer length, default is 7.
        rev_comp : bool
            Use reverse complement, default is False
        fmh_scale : int
            Scale factor for FracMinHash, default is 10
        threads: int
            The number of threads to use. Default is 16.
        """
		
		self.input_fasta = []
		self.max_prime = max_prime
		self.ksize = ksize
		self.rev_comp = rev_comp
		self.fmh_scale = fmh_scale
		self.threads = threads
		self.label = label
		# attributes for kmer load
		self.KO_names = []
		self.KO_cardinality = []
		self.dict_KO_to_index = {}  # to find the index of a KO
		# 20^7 = 10^9 * 10% = 10^8
		self.prefilter_bf = WritingBloomFilter(num_elements=10 ** 8, max_fp_prob=0.01, ignore_case=True, filename="prefilter_bf_"+label+".db")
		# TST object:
		self.kmer_dict = {}
		self.TST = mt.Trie()

	
	def print_stats(self):
		"""
        Prints out the stats of the KO.
        """
		print("The stats of this KO_TST object are:")
		print('Input fasta: \n%s\n' % "\n".join(self.input_fasta))
		print('Max Prime: %s \n' % self.max_prime)
		print('ksize: %s \n' % self.ksize)
		print('rev_comp: %s \n' % self.rev_comp)
		print('fmh_scale: %s \n' % self.fmh_scale)
		print('Label: %s \n' % self.label)
		print('Threads: %s \n' % self.threads)


	def load_kmer_from_fasta_into_dict(self, fasta_filename):
		"""
        Loads a KO_TST object from the input fasta.
        """
		if not os.path.exists(fasta_filename):
			raise ValueError('input fasta file does not exist')
		
		if fasta_filename in self.input_fasta:
			raise Exception("The input file has already been added!")
		
		print("Loading k-mers from : " + fasta_filename)
		time_start = time.time()
		self.input_fasta.append(fasta_filename)
		
		### temp variables
		max_prime = self.max_prime
		frac_prime = max_prime / self.fmh_scale
		ko_list = self.KO_names
		card_list = self.KO_cardinality
		dict_ko_index = self.dict_KO_to_index
		bf = self.prefilter_bf
		k_value = self.ksize
		kmer_dict = self.kmer_dict
		rev_comp = self.rev_comp
		
		### load the input fasta
		for record in SeqIO.parse(fasta_filename, 'fasta'):
			# get ko from string like "erl:AOC36_00005|chromosomal replication initiator DnaA|ko:K02313"
			current_ko = record.description.split('|')[-1].split(':')[-1]
			# check if this KO has already been added
			if current_ko in dict_ko_index:
				index_ko = dict_ko_index[current_ko]
			else:
				# add new KO to the last of the list
				index_ko = len(ko_list)
				ko_list.append(current_ko)
				card_list.append(0)
				dict_ko_index[current_ko] = index_ko
			# process sequence
			seq = str(record.seq)
			for i in range(len(seq) - k_value + 1):
				kmer = seq[i:i + k_value]
				if rev_comp:
					kmer = min(kmer, khmer.reverse_complement(kmer))
				
				# discard kmer above the frac ratio
				if khmer.hash_no_rc_murmur3(kmer) % max_prime > frac_prime:
					continue
					
				bf.add(kmer)
				# add kmer to the kmer_dict
				if kmer not in kmer_dict:
					kmer_dict[kmer] = {index_ko}
					card_list[index_ko] += 1
				else: # kmer has been found, check if it's from a known KO (dup)
					if index_ko not in kmer_dict[kmer]:
						kmer_dict[kmer].add(index_ko)
						card_list[index_ko] += 1
				
		### timeing
		time_end = time.time()
		print("Total time is: %s seconds \n" % round((time_end - time_start), 3))
		
		
	def find_a_kmer(self, kmer):
		"""
		Find if a kmer is in the database
		"""
		if self.rev_comp:
			kmer = min(kmer, khmer.reverse_complement(kmer))
			
		if kmer not in self.prefilter_bf:
			return False
		elif kmer not in self.kmer_dict:
			return False
		else:
			return [self.KO_names[x] for x in self.kmer_dict[kmer]]
		
		
	def stream_a_query_file(self, query_file, type, out_name="CI_out.csv"):
		"""
        Stream a fasta or fastq query file for CI in all KOs
        Type: fasta or fastq
        """
		if not os.path.exists(query_file):
			raise ValueError('input file does not exist')
		
		print("Start streaming query file: " + query_file)
		time_start = time.time()
		
		k_value = self.ksize
		prefilter = self.prefilter_bf
		kmer_dict = self.kmer_dict
		rev_comp = self.rev_comp
		
		temp_bf = WritingBloomFilter(num_elements=10 ** 8, max_fp_prob=0.01, ignore_case=True)
		hit_array = [0] * len(self.KO_names)
		
		for record in SeqIO.parse(query_file, type):
			seq = str(record.seq)
			for i in range(len(seq) - k_value + 1):
				kmer = seq[i:i + k_value]
				if rev_comp:
					kmer = min(kmer, khmer.reverse_complement(kmer))
				if kmer in prefilter:  # prefilter to accelerate unmatches
					if kmer not in temp_bf: # this is not a kmer that has already been added
						if kmer in kmer_dict:
							temp_bf.add(kmer)
							for index in kmer_dict[kmer]:
								hit_array[index] += 1
							
		ci_array = [round(i / (j+0.0000001), 3) for i, j in zip(hit_array, self.KO_cardinality)]
		temp_df = pd.DataFrame({'KO_name':self.KO_names, 'CI':ci_array})
		out_df = temp_df.sort_values('CI', ascending=False)
		out_df.to_csv(out_name, header=True, index=False)
		#del temp_bf, hit_array, temp_df, out_df
		
		### timeing
		time_end = time.time()
		print("Total time is: %s seconds \n\n" % round((time_end - time_start), 3))
	
		
	def no_prefilter_stream_a_query_file(self, query_file, type, out_name="CI_out.csv"):
			"""
	        Stream a fasta or fastq query file for CI in all KOs
	        Type: fasta or fastq
	        """
			if not os.path.exists(query_file):
				raise ValueError('input file does not exist')
			
			print("Start streaming query file: " + query_file)
			time_start = time.time()
			
			k_value = self.ksize
			prefilter = self.prefilter_bf
			kmer_dict = self.kmer_dict
			rev_comp = self.rev_comp
			
			temp_bf = WritingBloomFilter(num_elements=10 ** 8, max_fp_prob=0.01, ignore_case=True)
			hit_array = [0] * len(self.KO_names)
			
			for record in SeqIO.parse(query_file, type):
				seq = str(record.seq)
				for i in range(len(seq) - k_value + 1):
					kmer = seq[i:i + k_value]
					if rev_comp:
						kmer = min(kmer, khmer.reverse_complement(kmer))
					if kmer not in temp_bf:  # this is not a kmer that has already been added
						if kmer in kmer_dict:
							temp_bf.add(kmer)
							for index in kmer_dict[kmer]:
								hit_array[index] += 1
			
			ci_array = [round(i / (j+0.0000001), 3) for i, j in zip(hit_array, self.KO_cardinality)]
			temp_df = pd.DataFrame({'KO_name': self.KO_names, 'CI': ci_array})
			out_df = temp_df.sort_values('CI', ascending=False)
			out_df.to_csv(out_name, header=True, index=False)
			# del temp_bf, hit_array, temp_df, out_df
			
			### timeing
			time_end = time.time()
			print("Total time is: %s seconds \n\n" % round((time_end - time_start), 3))
		
		
	def export_to_pkl(self, export_file_name):
		"""
        Explort to a hdf5 file.
        """
		with open(export_file_name, 'wb') as outp:
			pickle.dump(self, outp, pickle.HIGHEST_PROTOCOL)
		

def local_tests():
	# two test files
	aa_fasta = "/Users/shaopeng/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Github/CMash_Rust_implementation/Amino_acid_index/top200.faa"
	small_fasta = "/Users/shaopeng/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Github/CMash_Rust_implementation/Amino_acid_index/top3.faa"
	nt_fastq = "/Users/shaopeng/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Github/CMash_Rust_implementation/Amino_acid_index/test_nt.fastq"
	fivek = "/Users/shaopeng/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Github/CMash_Rust_implementation/Amino_acid_index/fivek_record.faa"

	# local parameters
	ref_genome_files = [small_fasta, aa_fasta]
	ksize = 7
	max_prime = 9999999999971.
	fmh_scale = 10
	threads = 16
	rev_comp = False
	label = "KO_output"
	
	# test object
	fss = KO_TST()
	fss.print_stats()
	fss.load_kmer_from_fasta_into_dict(small_fasta)
	# check the KOs after loading
	print(fss.KO_names)
	print(fss.KO_cardinality)
	# load the large file with 200 records
	fss.load_kmer_from_fasta_into_dict(aa_fasta)  #0.14s for 200 records
	len(fss.KO_names)
	# find kmer
	fss.find_a_kmer('LSDGLDI')
	
	# stream a file
	fss.stream_a_query_file(aa_fasta, type='fasta', out_name="CI_out.csv")
	fss.no_prefilter_stream_a_query_file(aa_fasta, type='fasta', out_name="CI_out.csv")

	# store and load
	del fss.prefilter_bf  #pick can't handle BF obj
	fss.export_to_pkl(export_file_name="test_store.pkl")
	# load object
	with open("fmh_scale_10_KO_output.pkl", "rb") as fp:
		lsp = pickle.load(fp)
	# load bloom filter
	lsp.prefilter_bf = ReadingBloomFilter("prefilter_bf_KO_output.db")
	lsp.find_a_kmer('LSDGLDI')
	lsp.stream_a_query_file(fivek, type='fasta', out_name="CI_out.csv")
	lsp.no_prefilter_stream_a_query_file(fivek, type='fasta', out_name="CI_out.csv")
	
	# write and read BF filter
	aa = WritingBloomFilter(100, 0.01, ignore_case=True, filename="test_bf_100.bf")
	aa.add("AAA")
	aa.add("BBB")
	aa.add("CCC")
	bb = ReadingBloomFilter("test_bf_100.bf")
	

	
	
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Build AA FracMinHash sketches from a list of input ",
	                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-g', '--genome', type=str, help="File containing paths of ref sequences")
	parser.add_argument('-k', '--ksize', type=int, help="k-mer length", default=7)
	parser.add_argument('-p', '--max_prime', type=int, help="Prime value for hashing", default=9999999999971.)
	parser.add_argument('-s', '--fmh_scale', type=int, help="FracMinHash scale factor", default=10)
	parser.add_argument('-t', '--threads', type=int, help="thread number", default=16)
	parser.add_argument('-r', '--rev_comp', type=str, help="Use canonical kmer", default='False')
	parser.add_argument('-l', '--label', type=str, help="Keyword for out file", default='KO_record')
	
	### read parameters
	args = parser.parse_args()
	ref_genome_files = check_files(args.genome)
	if len(ref_genome_files) == 0:
		raise Exception("Didn't find genome!")
	
	ksize = args.ksize
	max_prime = args.max_prime
	fmh_scale = args.fmh_scale
	threads = args.threads
	label = args.label
	
	rev_comp = args.rev_comp
	rev_comp = rev_comp == 'True'
	if rev_comp:
		print("Using canonical kmers!")
		
		
		
	### build the TST object
	out_TST_obj = KO_TST(max_prime=max_prime, rev_comp=rev_comp, ksize=ksize, threads=threads, fmh_scale=fmh_scale, label=label)
	
	# load ref genomes into the TST object
	temp_start = time.time()
	for ref_genome in ref_genome_files:
		out_TST_obj.load_kmer_from_fasta_into_dict(fasta_filename=ref_genome)
	temp_end = time.time()
	print("Total loading time is %s seconds\n\n\n\n\n" % round((temp_end - temp_start), 3))
	del temp_start, temp_end
	
	# print status
	out_TST_obj.print_stats()
	
	# test run
	test_file = "/Users/shaopeng/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Github/CMash_Rust_implementation/Amino_acid_index/fivek_record.faa"
	out_TST_obj.stream_a_query_file(test_file, type="fasta", out_name="test_scan.csv")
	
	# save the object
	temp_start = time.time()
	print("Saving the object to pickle file")
	del out_TST_obj.prefilter_bf
	out_TST_obj.export_to_pkl(export_file_name="fmh_scale_"+str(out_TST_obj.fmh_scale)+"_"+out_TST_obj.label+'.pkl')
	temp_end = time.time()
	print("Total saving time is %s seconds\n\n\n\n\n" % round((temp_end - temp_start), 3))
	del temp_start, temp_end
	
	print("Job finished!")
	
		
	
		
	
	
	
