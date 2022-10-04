



## Documentation for AA FracMinHash sketching

### Python script:

1. [main object and wrapper to build reference database](https://github.com/KoslickiLab/CMash_Rust_implementation/blob/main/src/build_merged_TST_from_ref_AA_w_KOs.py)
2. [wrapper to stream query files](https://github.com/KoslickiLab/CMash_Rust_implementation/blob/main/src/import_db_and_stream_query_file.py)



### Algorithm overview:

---

Parameters:

| Param     | value          | Explaination                                           |
| --------- | -------------- | ------------------------------------------------------ |
| max_prime | 9999999999971. | A large prime number to mod for random permutation     |
| ksize     | 7              | k-mer length for amino acids                           |
| rev_comp  | False          | Use canonical k-mer                                    |
| fmh_scale | 10             | FracMinHash fraction: keep 1/fmh_scale of total k-mers |
| thread    | 16             | number of threads to use (not for here)                |
| label     | KO_term        | keyword for output files (ref_db)                      |



We will maintain a [KO_TST object](https://github.com/KoslickiLab/CMash_Rust_implementation/blob/a89244b9d3acd43aee327270ddaf0301e2d94bd2/src/build_merged_TST_from_ref_AA_w_KOs.py#L65):

1. key attributes:

```
1) ko_list: an array for all KO terms
2) card_list: an array for cardiority all KO terms, share the same index as 1)
3) dict_ko_index: a dict for KO:index pair that can quickly locate the index of a given KO
4) bf: a bloom filter as prefilter for all k-mers (size = 10^8, as 20^7 = 10^9)
5) kmer_dict: a dict for <k-mer> : <set of KO indices>, which will be used for data streaming
```



2. [Build ref-db](https://github.com/KoslickiLab/CMash_Rust_implementation/blob/a89244b9d3acd43aee327270ddaf0301e2d94bd2/src/build_merged_TST_from_ref_AA_w_KOs.py#L117)

```
For every record:
1. get the KO, locate in the ko_list OR append to the ko_list and card_list
2. get all k-mers in this record, for every k-mer
		a) use FracMinHash cutoff to select
		b) if it's a new kmer (not in kmer_dict), add to dict with ko, the card_list += 1
		c) if it exists, update the record with current ko
```



3. [Stream query file](https://github.com/KoslickiLab/CMash_Rust_implementation/blob/a89244b9d3acd43aee327270ddaf0301e2d94bd2/src/build_merged_TST_from_ref_AA_w_KOs.py#L196)

```
Maintain a temporary bloomfilter for dup count (temp_bf)
Maintain an array of hits for each KO (hit_array)

For every k-mer:
1. keep if in prefilter and NOT in temp_bf 
2. if it's in kmer_dict (double-check set identity, so the set size of pre-filter is less concerned), add to temp_bf and update hit_array (+1) based on indices attached to the k-mer 

Finally, the hit_array is the count of presence and can be used to get CI for all KOs.
```







