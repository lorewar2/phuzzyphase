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
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bcf::{self, Read as BcfRead};
use std::path::Path;

use hashbrown::{HashMap, HashSet};

use clap::App;

const K: usize = 6; // kmer match length
const W: usize = 20; // Window size for creating the band
const MATCH: i32 = 1; // Match score
const MISMATCH: i32 = -5; // Mismatch score
const GAP_OPEN: i32 = -5; // Gap open score
const GAP_EXTEND: i32 = -1; // Gap extend score

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
    //if params.long_read_bam == None && params.linked_read_bam == None && params.hic_bam == None {
    //    eprintln!("Must supply at least one bam");
    //    std::process::exit(1);
    //}
    let fai = params.fasta.to_string() + ".fai";
    let fa_index_iter = fasta::Index::from_file(&fai)
        .expect(&format!("error opening fasta index: {}", fai))
        .sequences();
    let mut chroms: Vec<String> = Vec::new();
    let mut chrom_lengths: Vec<u64> = Vec::new();
    for chrom in fa_index_iter {
        chroms.push(chrom.name.to_string());
        chrom_lengths.push(chrom.len);
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
            min_mapq: params.min_mapq,
            min_base_qual: params.min_base_qual,
            chrom_length: chrom_lengths[i],
            window: params.window,
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
    let mut long_read_bam_reader = match &data.long_read_bam {
        Some(x) => Some(bam::IndexedReader::from_path(x)?),
        None => None,
    };
    let mut linked_read_bam_reader = match &data.linked_read_bam {
        Some(x) => Some(bam::IndexedReader::from_path(x)?),
        None => None,
    };
    let mut hic_bam_reader = match &data.hic_bam {
        Some(x) => Some(bam::IndexedReader::from_path(x)?),
        None => None,
    };
    let mut fasta = fasta::IndexedReader::from_file(&data.fasta).expect("cannot open fasta file");

    let mut molecule_alleles: MoleculeAlleleWrapper = MoleculeAlleleWrapper{
        long_read_assignments: match &data.long_read_bam {
            Some(_x) => Some(HashMap::new()),
            None => None,
        },
        linked_read_assignments: match &data.linked_read_bam {
            Some(_x) => Some(HashMap::new()),
            None => None,
        },
        
        hic_read_assignments: match &data.hic_bam {
            Some(_x) => Some(HashMap::new()),
            None => None,
        },
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
        let alleles = rec.alleles();
        if alleles.len() > 2 {
            continue; // ignore multi allelic sites
        }
        let reference = std::str::from_utf8(alleles[0])?;
        let alternative = std::str::from_utf8(alleles[1])?;
        let rec_rid = rec.rid().expect("could not unwrap vcf record id");

        let genotypes = rec.genotypes()?;
        let genotype = genotypes.get(0); // assume only 1 and get the first one
        if is_heterozygous(genotype) {
            hets += 1;
            get_molecule_allele_assignments(
                &data.chrom,
                pos as usize,
                reference.to_string(),
                alternative.to_string(),
                data.min_mapq,
                data.min_base_qual,
                &mut long_read_bam_reader,
                &mut linked_read_bam_reader,
                &mut hic_bam_reader,
                &mut fasta,
                data.window,
                data.chrom_length,
                &mut molecule_alleles,
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
    chrom: &String,
    pos: usize,
    ref_allele: String,
    alt_allele: String,
    min_mapq: u8,
    min_base_qual: u8,
    long_reads: &mut Option<bam::IndexedReader>,
    linked_reads: &mut Option<bam::IndexedReader>,
    hic_reads: &mut Option<bam::IndexedReader>,
    fasta: &mut fasta::IndexedReader<std::fs::File>,
    window: usize,
    chrom_length: u64,
    molecule_alleles: &mut MoleculeAlleleWrapper,
) {
    if (pos + window) as u64 > chrom_length {
        return;
    }
    match long_reads {
        Some(bam) => {
            let tid = bam.header().tid(chrom.as_bytes()).expect("cannot find chrom tid");
            bam.fetch((tid, pos as u32, (pos + 1) as u32)).expect("blah");
            let ref_start = (pos - window) as u64;
            let ref_end = (pos + window + ref_allele.len()) as u64;
            fasta.fetch(chrom, ref_start, ref_end);
            let mut ref_sequence: Vec<u8> = Vec::new();
            fasta.read(&mut ref_sequence).expect("failed to read fasta sequence");

            let mut alt_sequence: Vec<u8> = Vec::new();
            for i in 0..window as usize {
                alt_sequence.push(ref_sequence[i]);
            }
            for base in alt_allele.as_bytes() {
                alt_sequence.push(*base);
            }
            for i in (window as usize+ref_allele.len())..ref_sequence.len() {
                alt_sequence.push(ref_sequence[i]);
            }
            
            println!("double checking. ref allele {} alt {} chrom {} pos {}", ref_allele, alt_allele, chrom, pos);
            println!("ref_sequence {}", std::str::from_utf8(&ref_sequence).unwrap());
            println!("alt_sequence {}", std::str::from_utf8(&alt_sequence).unwrap());
            println!("");
            for _rec in bam.records() {
                let rec = _rec.expect("cannot read bam record");
                println!("\n{}", std::str::from_utf8(rec.qname()).unwrap());
                if rec.mapq() < min_mapq {
                    continue
                }
                if rec.is_secondary() || rec.is_supplementary() {
                    continue;
                }
                let mut read_start: usize = 10000000000; // I am lazy and for some reason dont know how to do things, so this is my bad solution
                let mut read_end: usize = rec.seq_len();
                let mut min_bq = 93;
                let qual = rec.qual();
                for pos_pair in rec.aligned_pairs() {
                    if (pos_pair[1] as u64) >= ref_start && (pos_pair[1] as u64) < ref_end {
                        if pos_pair[1] as usize >= pos && pos_pair[1] as usize <= pos + ref_allele.len().max(alt_allele.len()) {
                            min_bq = min_bq.min(qual[pos_pair[0] as usize]);
                        }
                        read_end = pos_pair[0] as usize;
                        if read_start == 10000000000 {
                            read_start = pos_pair[0] as usize;
                        }
                    } 
                }
                if min_bq < min_base_qual {
                    continue;
                }
                if read_start >= read_end {
                    eprintln!("what happened, read start {} read end {}",read_start, read_end);
                    continue;
                }
                let seq = rec.seq().encoded[read_start..read_end].to_vec();
                eprintln!("ref sequence {}", std::str::from_utf8(&ref_sequence).unwrap());
                eprintln!("alt sequence {}", std::str::from_utf8(&alt_sequence).unwrap());
                eprintln!("read segment {}\n", std::str::from_utf8(&seq).unwrap());
                let score = |a: u8, b: u8| if a == b { MATCH } else { MISMATCH };
                let mut aligner = banded::Aligner::new(GAP_OPEN, GAP_EXTEND, score, K, W);
                let ref_alignment = aligner.local(&seq, &ref_sequence);
                let alt_alignment = aligner.local(&seq, &alt_sequence);
                if ref_alignment.score > alt_alignment.score {
                    eprintln!("read supports ref allele");
                    match &mut molecule_alleles.long_read_assignments {
                        Some(long_read_assignment) => {
                            let readdata = long_read_assignment.entry(std::str::from_utf8(rec.qname()).expect("readname fail").to_string()).or_insert(Vec::new());
                            readdata.push(Allele{
                                locus: pos,
                                allele: true,
                            });
                        },
                        None => (),
                    }
                } else if alt_alignment.score > ref_alignment.score {
                    eprintln!("read supports alt allele");
                    match &mut molecule_alleles.long_read_assignments {
                        Some(long_read_assignment) => {
                            let readdata = long_read_assignment.entry(std::str::from_utf8(rec.qname()).expect("readname fail").to_string()).or_insert(Vec::new());
                            readdata.push(Allele{
                                locus: pos,
                                allele: false,
                            });
                        },
                        None => (),
                    }
                } else {
                    eprintln!("\nread had equal alignment scores ref {} alt {}", ref_allele, alt_allele);
                    
                }
                
            }
        }
        None => (),
    }
}

fn is_heterozygous(gt: bcf::record::Genotype) -> bool {
    if gt[0] == bcf::record::GenotypeAllele::Unphased(0)
        && gt[1] == bcf::record::GenotypeAllele::Unphased(1)
    {
        return true;
    } else if gt[0] == bcf::record::GenotypeAllele::Unphased(1)
        && gt[1] == bcf::record::GenotypeAllele::Unphased(0)
    {
        return true;
    } else if gt[0] == bcf::record::GenotypeAllele::Unphased(0)
        && gt[1] == bcf::record::GenotypeAllele::Phased(1)
    {
        return true;
    } else if gt[0] == bcf::record::GenotypeAllele::Unphased(1)
        && gt[1] == bcf::record::GenotypeAllele::Phased(0)
    {
        return true;
    }
    return false;
}



struct MoleculeAlleleWrapper {
    long_read_assignments: Option<HashMap<String, Vec<Allele>>>,
    linked_read_assignments: Option<HashMap<String, Vec<Allele>>>,
    hic_read_assignments: Option<HashMap<String, Vec<Allele>>>,
}

struct Allele {
    locus: usize,
    allele: bool, // ref is 0 and alt is 1
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
    chrom_length: u64,
    vcf: String,
    min_mapq: u8,
    min_base_qual: u8,
    window: usize,
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
    min_mapq: u8,
    min_base_qual: u8,
    window: usize,
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

    let min_mapq = params.value_of("min_mapq").unwrap_or("30");
    let min_mapq = min_mapq.to_string().parse::<u8>().unwrap();
    
    let min_base_qual = params.value_of("min_base_qual").unwrap_or("20");
    let min_base_qual = min_base_qual.to_string().parse::<u8>().unwrap();

    let restarts = params.value_of("restarts").unwrap_or("10");
    let restarts = restarts.to_string().parse::<u32>().unwrap();

    let fasta = params.value_of("fasta").unwrap();

    let ploidy = params.value_of("ploidy").unwrap_or("2");
    let ploidy = ploidy.to_string().parse::<usize>().unwrap();

    let allele_alignment_window = params.value_of("ploidy").unwrap_or("200");
    let allele_alignment_window = allele_alignment_window.to_string().parse::<usize>().unwrap();

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
        min_mapq: min_mapq,
        min_base_qual: min_base_qual,
        window: allele_alignment_window,
    }
}
