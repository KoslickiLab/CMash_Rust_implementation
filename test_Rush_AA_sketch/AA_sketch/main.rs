use std::env;
use bio::io::fasta;
use std::collections::HashMap;
use std::time::{Instant};



fn main() {
    println!("Running Rust with AA indexing");

    // read parameter
    let args: Vec<String> = env::args().collect();
    let input_file = &args[1];
    let ksize = &args[2];
    let _out_file = &args[3];
    let now = Instant::now();
    println!("The input file is {}", input_file);


    // record only
    let reader = fasta::Reader::from_file(input_file).unwrap();

    for result in reader.records() {

        let result_data = &result.unwrap();

        let sequences = result_data.seq();

        // create a hashmap
        let mut kmer_hm: HashMap<&[u8], Vec<usize>> = HashMap::new();

        // generate kmer and add to dict
        let iter = sequences.windows(ksize.parse().unwrap());
        for (i, kmer) in iter.enumerate() {
            if kmer_hm.contains_key(kmer) {
                // update value: append to array
                kmer_hm.entry(kmer).and_modify(|e| { e.push(i) });
            }
            else {
                kmer_hm.insert(kmer, vec![i]);
            }

        }


        // // print the hashmap for test
        // println!("Print hashmap!!!");
        // for (aa, loc) in kmer_hm.iter() {
        //     println!("aa seq: {:?}, loc: {:?}", aa, loc);
        // }

        // write to json?
        

        
    }

    println!("Time elapsed: {:?}", now.elapsed());



}
