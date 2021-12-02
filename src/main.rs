#[macro_use]
extern crate clap;
extern crate bio;
extern crate hashbrown;
extern crate rand;
extern crate rayon;

use flate2::read::MultiGzDecoder;
use bio::alignment::pairwise::banded;
use bio::io::fasta;
use bio::utils::TextSlice;
use std::path::Path;
use failure::{Error, ResultExt};
use rust_htslib::bam::{self, Read, Record};
use rust_htslib::bcf::{self, Read as BcfRead};
use rayon::prelude::*;

use hashbrown::{HashMap, HashSet};

use clap::App;

fn main() {
   let res = _main();
   if let Err(v) = res {
        let fail = v.as_fail();
        println!("Phasstphase error. v{}.", env!("CARGO_PKG_VERSION"));
        println!("Error: {}", fail);

        for cause in fail.iter_causes() {
            println!("Info: caused by {}", cause);
        }

        println!("\n{}", v.backtrace());

        println!("If you think this is bug in Phasstphase, please file a bug at https://github.com/wheaton5/phasstphase, and include the information above and the command-line you used.");
        std::process::exit(1)
    }
}


fn _main() -> Result<(), Error> {
    println!("Welcome to phasst phase!");
    let params = load_params();
    if params.long_read_bam == None && params.linked_read_bam == None && params.hic_bam == None {
        eprintln!("Must supply at least one bam");
        std::process::exit(1);
    }
    let fai = params.fasta.to_string() + ".fai";
    let fa_index_iter = fasta::Index::from_file(&fai)
        .expect(&format!("error opening fasta index: {}", fai))
        .sequences();
    let mut chroms: Vec<String> = Vec::new();
    for chrom in fa_index_iter {
        chroms.push(chrom.name.to_string());
    }
    
    

    let mut chunks: Vec<ThreadData> = Vec::new();
    for (i, chrom) in chroms.iter().enumerate() {
        let data = ThreadData{
            index: i,
            long_read_bam: match &params.long_read_bam {
                Some(x) => Some(x.to_string()),
                None => None,
            },
            linked_read_bam: match &params.linked_read_bam {
                Some(x) => Some(x.to_string()),
                None => None,
            },
            hic_bam: match &params.hic_bam {
                Some(x) => Some(x.to_string()),
                None => None,
            },
            fasta: params.fasta.to_string(),
            chrom: chrom.to_string(),
            vcf: params.vcf.to_string(),
        };
        chunks.push(data);
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(params.threads)
        .build()
        .unwrap();
    let results: Vec<_> = pool.install(|| {
        chunks.par_iter()
            .map(|rec_chunk| {
                phase_chunk(&rec_chunk)
            })
            .collect()
    });
    Ok(())
}

fn phase_chunk<'a>(data: &ThreadData) -> Result<(), Error> {
    println!("thread {} chrom {}", data.index, data.chrom);
    for i in 1..1000000000 {
        //eprintln!("{}", i);
        let b = (i as f64).sqrt();
    }
    let long_read_bam_reader = match &data.long_read_bam {
        Some(x) => Some(bam::IndexedReader::from_path(x)?),
        None => None,
    };
    let linked_read_bam_reader = match &data.linked_read_bam {
        Some(x) => Some(bam::IndexedReader::from_path(x)?),
        None => None,
    };
    let hic_bam_reader = match &data.hic_bam {
        Some(x) => Some(bam::IndexedReader::from_path(x)?),
        None => None,
    };
    Ok(())
}

struct ThreadData {
    index: usize,
    long_read_bam: Option<String>,
    linked_read_bam: Option<String>,
    hic_bam: Option<String>,
    fasta: String,
    chrom: String,
    vcf: String,
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