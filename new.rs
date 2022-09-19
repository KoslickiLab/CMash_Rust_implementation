use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use std::collections::HashMap;

//use std::str::Split;

//use std::fs;

fn _hash_maps(){
    let mut test_hm = HashMap::new();
    test_hm.insert("key1",vec!["value1"]);
    test_hm.insert("key2",vec!["value2"]);
    test_hm.insert("key3",vec!["value3","value4"]);
    //println!("{}", test_hm["key1"]);
    println!("{:?}", test_hm["key1"]);

    if test_hm.contains_key("key3"){
        println!("{:?}", test_hm["key3"]);
        println!("{:?}", test_hm["key3"][1]); 
    }

    for (key,value) in &test_hm {
        println!("{} {:?}", key, value);
    }

}

fn _hashstr(data2:&str,k:usize)->Vec<u64>{

    //seperate into k-mers
    let seq_chars: Vec<char> = data2.chars().collect();
    let sequences: Vec<String> = seq_chars.windows(k)

    .map(|w| w.iter().collect::<String>()).collect();


    //find inverse of that k-mer
    let inverdata=data2.replace("A","t").replace("C","g").replace("T","a").replace("G","c")
    .replace("a","A").replace("c","C").replace("g","G").replace("t","T");

    let seq_chars2: Vec<char> = inverdata.chars().collect();
    let sequences2: Vec<String> = seq_chars2.windows(k)

    .map(|w| w.iter().collect::<String>()).collect();
 

    //hash the two we find
    let mut hashseq=Vec::new();
    for i in sequences {
        let mut hash = DefaultHasher::new();
        i.hash(&mut hash);

        hashseq.push(hash.finish());

    }

    let mut hashseqinver=Vec::new();
    for i in sequences2 {
        let mut hash = DefaultHasher::new();
        i.hash(&mut hash);
        
        hashseqinver.push(hash.finish());

    }


    //remain the smaller value (normal and its inverse) of a k-mer 
    let mut i=0;
    while i < hashseq.len()   {
        if hashseq[i] > hashseqinver[i] {
            hashseq[i] = hashseqinver[i];
        }
        i+=1;
    }

    //eliminate duplicate
    hashseq.sort();
    hashseq.dedup();

    return hashseq
}


fn main() {
    



    //let data = fs::read_to_string(r"C:\Users\Qijian\Desktop\rust\test2.txt").expect("failed to read file");




    let data1 = "ACTCGCAGATCAAAAAATGAAAAAA";
    let data2 = "TCGAACGTCAGATACTCGACCATCGA";
    let k=5;
    let hashseq2=_hashstr(data1,k);
    let hashseq1=_hashstr(data2,k);

    println!("{:?}\n", hashseq1);
    println!("{:?}\n", hashseq2);


    //This is to find how many similar object between 2 sequences
    let mut i2=0;
    let mut count=0;
    let len=hashseq2.len(); 
    'outer: for i in hashseq1{
        if i == hashseq2[i2]{
            count+=1;
            i2+=1;
            
        }
        else{
            if i2<len{
                while hashseq2[i2]<i {
                    i2+=1;
                    if i2>=len{break 'outer; }
                }
                if i == hashseq2[i2]{
                    println!("{:?}\n", hashseq2[i2]);
                    count+=1;
                    i2+=1;
                    if i2>=len{break 'outer; }
                }
            }
            else{
                break 'outer;
            }
        }
    }

    println!("A and B: {}", count);
    }