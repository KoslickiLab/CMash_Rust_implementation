# CMash production code via Rust
Develop production code for CMash, and apply it to metagenomic analyses.  
Contact: Shaopeng Liu (sml6467@psu.edu)



### Contents

- [Current production](#prod)
  1. [Python FracMinHash indexing for Amino acid by KO terms](#ko)

- [Ongoing tasks](#task)

  





<br>

### Current production <a name="prod"></a>

---

#### Python FracMinHash indexing for AA by KO <a name="ko"></a>

Part of our metagenome blueprint. We can build k-mer sketches for amino acids in all KO terms.

Please [check here](https://github.com/KoslickiLab/CMash_Rust_implementation/tree/main/production_py_AA_sketch/readme.md) for more details. 

<br>

A simple usage case in our GPU server (where I've already saved the ref-database):

1. build the "metagenome" conda env

2. use the following code:
   ```
   # fill your own variables:
   label=""  # keyword for output files
   query_list="" # a file containing abs paths of all query files, one per line
   input_type="" # fasta / fastq
   
   
   
   # here is an example
   label="test"
   query_list="/data/sml6467/github/CMash_Rust_implementation/output/test_data/query_list_3_small_fasta.txt"
   input_type="fasta"
   # a large file for speed testing
   # query_list2="/data/sml6467/github/CMash_Rust_implementation/output/test_data/query_list_1_file_1M.txt"
   # it takes ~156s to load the database, and ~360s to stream 1M records, and ~60s to sort and save csv
   
   
   
   # kmer sketch database
   py_script=/data/sml6467/github/CMash_Rust_implementation/src/import_db_and_stream_query_file.py
   ref_database=/data/sml6467/github/CMash_Rust_implementation/output/sketch_aa_7mer/fmh_scale_10_all_KO_AA.pkl
   bf_prefilter=/data/sml6467/github/CMash_Rust_implementation/output/sketch_aa_7mer/prefilter_bf_all_KO_AA.db
   
   
   
   # process data
   conda activate metagenome
   python ${py_script} -r ${ref_database} -b ${bf_prefilter} -l ${label} -i ${query_list} -p ${input_type}
   
   
   
   # each input file from query list will get a separate CSV output
   Sample output:
   Total loading time is 177.092 seconds
   
   Start streaming query file: /data/sml6467/github/CMash_Rust_implementation/output/test_data/record_1-10.faa
   Total time is: 0.176 seconds
   
   
   Start streaming query file: /data/sml6467/github/CMash_Rust_implementation/output/test_data/record_250-260.faa
   Total time is: 0.118 seconds
   
   
   Start streaming query file: /data/sml6467/github/CMash_Rust_implementation/output/test_data/record_50-60.faa
   Total time is: 0.113 seconds
   ```

   





<br>

### Ongoing tasks <a name="task"></a>	

---

- [ ] Complete Python functions for AA sketching (TST)
- [ ] Finish Rust Hashmap + list implementation of the Python code
- [ ] Add TST into Rust
- [ ] Add this as a feature of Sourmash and use their code base.

<br>




### Resources:
1. [CMash repo (Python)](https://github.com/dkoslicki/CMash)
2. [An example of CMash usage](https://github.com/KoslickiLab/CMASH-reproducibles)

