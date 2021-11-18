#[macro_use]
extern crate clap;
extern crate bio;
extern crate hashbrown;
extern crate rand;
extern crate rayon;

use flate2::read::MultiGzDecoder;
use bio::alignment::pairwise::banded;
use bio::io::fasta;
use bio::io::fasta::Record;
use bio::utils::TextSlice;
use std::path::Path;

use hashbrown::{HashMap, HashSet};

use clap::App;

fn main() {
    println!("Welcome to phasst phase!");
    let params = load_params();
}

#[derive(Clone)]
struct Params {
    linked_read_bam: Option<String>,
    hic_bam: Option<String>,
    long_read_bam: Option<String>,
    output: String,
    seed: u8,
    restarts: u32,
    threads: usize,
    fasta: String,
    ploidy: usize,
    vcf: String,
}


fn load_params() -> Params {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();

    let output = params.value_of("output").unwrap();
    let linked_read_bam = match params.value_of("linked_read_bam") {
        Some(x) => Some(x.to_string()),
        None => None
    };
    let hic_bam = match params.value_of("hic_bam") {
        Some(x) => Some(x.to_string()),
        None => None
    };
    let long_read_bam = match params.value_of("long_read_bam") {
        Some(x) => Some(x.to_string()),
        None => None
    };


    let threads = params.value_of("threads").unwrap_or("1");
    let threads = threads.to_string().parse::<usize>().unwrap();

    let seed = params.value_of("seed").unwrap_or("4"); // 4 is guarranteed random by dice roll https://xkcd.com/221/
    let seed = seed.to_string().parse::<u8>().unwrap();

    let restarts = params.value_of("restarts").unwrap_or("10");
    let restarts = restarts.to_string().parse::<u32>().unwrap();

    let fasta = params.value_of("fasta").unwrap();

    let ploidy = params.value_of("ploidy").unwrap_or("2");
    let ploidy = ploidy.to_string().parse::<usize>().unwrap();

    let vcf = params.value_of("vcf").unwrap();

    Params {
        output: output.to_string(),
        linked_read_bam: linked_read_bam,
        long_read_bam: long_read_bam,
        hic_bam: hic_bam,
        threads: threads,
        seed: seed,
        fasta: fasta.to_string(),
        ploidy: ploidy,
        vcf: vcf.to_string(),
        restarts: restarts,
    }
}