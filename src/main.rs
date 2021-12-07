#[macro_use]
extern crate clap;
extern crate bio;
extern crate hashbrown;
extern crate rand;
extern crate rayon;

use bio::alignment::pairwise::banded;
use bio::io::fasta;
use bio::utils::TextSlice;
use failure::{Error, ResultExt};
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use rust_htslib::bam::{self, Read, Record};
use rust_htslib::bcf::{self, Read as BcfRead};
use std::path::Path;

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
        let data = ThreadData {
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
        chunks
            .par_iter()
            .map(|rec_chunk| phase_chunk(&rec_chunk))
            .collect()
    });
    Ok(())
}

fn phase_chunk(data: &ThreadData) -> Result<(), Error> {
    println!("thread {} chrom {}", data.index, data.chrom);
    let molecule_alleles = get_molecule_alleles_assignments(data);
    Ok(())
}

fn get_molecule_alleles_assignments(data: &ThreadData) -> Result<MoleculeAllelesWrapper, Error> {
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

   

    let mut vcf_reader = bcf::IndexedReader::from_path(data.vcf.to_string())?;
    let rid = vcf_reader.header().name2rid(data.chrom.as_bytes())?;
    vcf_reader.fetch(rid, 0, None)?; 
    let mut total = 0;
    let mut hets = 0;
    for (i, _rec) in vcf_reader.records().enumerate() {
        total += 1;
        let rec = _rec?;
        let pos = rec.pos();
        let rec_rid = rec.rid().expect("could not unwrap vcf record id");
        if rec_rid != rid {
            println!("vcf read off end of chromosome {} to {:?}", data.chrom, 
                std::str::from_utf8(vcf_reader.header().rid2name(rec_rid)?));
            break;
        }
        let genotypes = rec.genotypes()?;
        let genotype = genotypes.get(0); // assume only 1 and get the first one
        println!("{:?}", genotype);
        if is_heterozygous(genotype) {
            hets += 1;
            let molecule_assignments = get_molecule_allele_assignments(
                &long_read_bam_reader,
                &linked_read_bam_reader,
                &hic_bam_reader,
            );
        }
    }
    println!("done, saw {} records of which {} were hets in chrom {}", total, hets, data.chrom);
    Ok(MoleculeAllelesWrapper {
        hic_alleles: None,
        long_read_alleles: None,
        linked_read_alleles: None,
    })
}

fn get_molecule_allele_assignments(
    long_reads: &Option<bam::IndexedReader>,
    linked_reads: &Option<bam::IndexedReader>,
    hic_reads: &Option<bam::IndexedReader>,
) -> MoleculeAlleleWrapper {

    MoleculeAlleleWrapper{
        long_read_assignments: None,
        linked_read_assignments: None,
        hic_read_assignments: None,
    }
}

fn is_heterozygous(gt: bcf::record::Genotype) -> bool {
    if gt[0] == bcf::record::GenotypeAllele::Unphased(0)
        && gt[1] == bcf::record::GenotypeAllele::Unphased(1)
    {
        println!("0/1");
        return true;
    } else if gt[0] == bcf::record::GenotypeAllele::Unphased(1)
        && gt[1] == bcf::record::GenotypeAllele::Unphased(0)
    {
        println!("1/0");
        return true;
    } else if gt[0] == bcf::record::GenotypeAllele::Unphased(0)
        && gt[1] == bcf::record::GenotypeAllele::Phased(1)
    {
        println!("0|1");
        return true;
    } else if gt[0] == bcf::record::GenotypeAllele::Unphased(1)
        && gt[1] == bcf::record::GenotypeAllele::Phased(0)
    {
        println!("1|0");
        return true;
    }
    println!("not heterozygous");
    return false;
}

struct MoleculeAlleleWrapper {
    long_read_assignments: Option<ReadAssignments>,
    linked_read_assignments: Option<ReadAssignments>,
    hic_read_assignments: Option<ReadAssignments>,
}

struct ReadAssignments {
    assignments: Vec<ReadAssignment>,
}

struct ReadAssignment {
    id: String,
    allele: bool,
}


struct MoleculeAllelesWrapper {
    hic_alleles: Option<Vec<MoleculeAlleles>>,
    long_read_alleles: Option<Vec<MoleculeAlleles>>,
    linked_read_alleles: Option<Vec<MoleculeAlleles>>,
}

struct MoleculeAlleles {}

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
        None => None,
    };
    let hic_bam = match params.value_of("hic_bam") {
        Some(x) => Some(x.to_string()),
        None => None,
    };
    let long_read_bam = match params.value_of("long_read_bam") {
        Some(x) => Some(x.to_string()),
        None => None,
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
