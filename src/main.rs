#[macro_use]
extern crate clap;
extern crate bio;
extern crate hashbrown;
extern crate rand;
extern crate rayon;

use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use std::process::Command;
use std::{thread, time};


use bio::alignment::pairwise::banded;
use bio::io::fasta;
use bio::utils::TextSlice;
use failure::{Error, ResultExt};
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use rust_htslib::bam::{self, Read, Record};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bcf::{self, Read as BcfRead};
use rust_htslib::bcf::{Format};
use std::path::Path;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::convert::TryInto;


use hashbrown::{HashMap, HashSet};

use clap::App;

const K: usize = 6; // kmer match length
const W: usize = 20; // Window size for creating the band
const MATCH: i32 = 1; // Match score
const MISMATCH: i32 = -5; // Mismatch score
const GAP_OPEN: i32 = -5; // Gap open score
const GAP_EXTEND: i32 = -1; // Gap extend score

fn main() {
    let result = _main();
    if let Err(v) = result {
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
    if Path::new(&params.output).is_dir() {
        eprintln!("restarting from partial output in {}", params.output);
    } else {
        fs::create_dir(params.output.to_string())?;
    }
    
    let fai = params.fasta.to_string() + ".fai";
    let fa_index_iter = fasta::Index::from_file(&fai)
        .expect(&format!("error opening fasta index: {}", fai))
        .sequences();
    let mut chroms: Vec<String> = Vec::new();
    chroms.push("20".to_string());
    let mut chrom_lengths: Vec<u64> = Vec::new();
    chrom_lengths.push(63025520);
    /*for chrom in fa_index_iter {
        chroms.push(chrom.name.to_string());
        chrom_lengths.push(chrom.len);
    }
    */
    let mut vcf_reader = bcf::IndexedReader::from_path(params.vcf.to_string())?;
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
            output: params.output.to_string(),
            vcf_out: format!("{}/chrom_{}.vcf", params.output, chrom),
            vcf_out_done: format!("{}/chrom_{}.vcf.done", params.output, chrom),
            phasing_window: params.phasing_window,
            seed: params.seed,
            ploidy: params.ploidy,
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
    
    if !Path::new(&data.vcf_out_done).exists() {
        get_all_variant_assignments(data);
    }
    
    let mut vcf_reader = bcf::IndexedReader::from_path(format!("{}.gz",data.vcf_out.to_string()))
        .expect("could not open indexed vcf reader on output vcf");
    let chrom = vcf_reader.header().name2rid(data.chrom.as_bytes()).expect("cant get chrom rid");
    let (num_variants, last_position) = inspect_vcf(&mut vcf_reader, &data);
    let mut cluster_centers = init_cluster_centers(num_variants, &data);
    let mut window_start: usize = 60000; //TODO dont do this lol, for testing
    let mut window_end: usize = 60000+data.phasing_window;
    let mut position_to_index: HashMap<usize, usize> = HashMap::new();
    let mut position_so_far: usize = 0;
    let mut cluster_center_delta: f32 = 10.0;
    while (window_start as u64) < data.chrom_length && window_start < last_position as usize {
        
        vcf_reader.fetch(chrom, window_start as u64, Some(window_end as u64)).expect("some actual error");
        println!("fetching region {}:{}-{}", data.chrom, window_start, window_end);
        let molecules = get_long_read_molecules(&mut vcf_reader, &mut position_to_index, &mut position_so_far);

        let mut iteration = 0;
        for haplotype in 0..cluster_centers.len() {
            print!("haplotype {}\t", haplotype);
            for variant in 0..10 {
                print!("var{}:{} | ", variant, cluster_centers[haplotype][variant]);
            }println!();
        }
        while cluster_center_delta > 0.01 {
            let posteriors = expectation(&molecules, &cluster_centers);
            cluster_center_delta = maximization(&molecules, posteriors, &mut cluster_centers);
            for haplotype in 0..cluster_centers.len() {
                print!("haplotype {}\t", haplotype);
                for variant in 0..10 {
                    print!("var{}:{} | ", variant, cluster_centers[haplotype][variant]);
                }println!();
            }
            println!("cluster center delta {}", cluster_center_delta);
            iteration += 1;
        }
        println!("converged in {} iterations", iteration);
        for haplotype in 0..cluster_centers.len() {
            print!("haplotype {}\t", haplotype);
            for variant in 0..10 {
                print!("var{}:{} | ", variant, cluster_centers[haplotype][variant]);
            }println!();
        }
        window_start += data.phasing_window/4;
        window_start = window_start.min(last_position as usize);
        window_end += data.phasing_window/4;
        window_end = window_end.min(last_position as usize);
        //break;
    }
    
    Ok(())
}

fn maximization(molecules: &Vec<Vec<Allele>>, posteriors: Vec<Vec<f32>>, cluster_centers: &mut Vec<Vec<f32>>) -> f32 {
    let mut updates: HashMap<usize, Vec<(f32, f32)>> = HashMap::new(); // variant index to vec across 
    let mut variant_molecule_count: HashMap<usize, usize> = HashMap::new();
    // haplotype clusters to a tuple of (numerator, denominator)
    for molecule_index in 0..molecules.len() {
        // if molecule does not support any haplotype over another, dont use it in maximization
        let mut different = false;
        for haplotype in 0..cluster_centers.len() {
            if posteriors[molecule_index][haplotype] != posteriors[molecule_index][0] { 
                different = true; 
            }
        }
        if !different { continue; }
        for allele in &molecules[molecule_index] {
            let variant_index = allele.index;
            let alt = allele.allele;
            let count = variant_molecule_count.entry(variant_index).or_insert(0);
            *count += 1;
            for haplotype in 0..cluster_centers.len() {
                let posterior = posteriors[molecule_index][haplotype];
                let numerators_denominators = updates.entry(variant_index).or_insert(vec![(0.0, 0.0); cluster_centers.len()]);
                if alt {
                    numerators_denominators[haplotype].0 += posterior;
                    numerators_denominators[haplotype].1 += posterior;
                } else {
                    numerators_denominators[haplotype].1 += posterior;
                }
            }
        }
    }
    let mut total_change = 0.0;
    for (variant_index, haplotypes) in updates.iter() {
        if variant_molecule_count.get(variant_index).unwrap() > &3 { //TODO Dont hard code stuff
            for haplotype in 0..haplotypes.len() {
                let numerators_denominators = haplotypes[haplotype];
                let allele_fraction = (numerators_denominators.0 / numerators_denominators.1).max(0.001).min(0.999);
                total_change += (allele_fraction - cluster_centers[haplotype][*variant_index]).abs();
                cluster_centers[haplotype][*variant_index] = allele_fraction;
            }
        }
    }
    total_change
}

// molecules is vec of molecules each a vec of alleles
// cluster centers is vec of haplotype clusters  by variant loci
// posteriors (return value) is molecules by clusters to posterior probability
fn expectation(molecules: &Vec<Vec<Allele>>, cluster_centers: &Vec<Vec<f32>>) -> Vec<Vec<f32>> {
    let mut posteriors: Vec<Vec<f32>> = Vec::new();

    for molecule in molecules.iter() {
        let mut log_probs: Vec<f32> = Vec::new(); // for each haplotype
        for haplotype in 0..cluster_centers.len() {
            let mut log_prob = 0.0; // log(0) = probability 1
            for allele in molecule.iter() {
                if allele.allele { // alt allele
                    log_prob += cluster_centers[haplotype][allele.index].ln();// adding in log space, multiplying in probability space
                } else {
                    log_prob += (1.0 - cluster_centers[haplotype][allele.index]).ln();
                }
            }
            log_probs.push(log_prob);
        }
        let bayes_log_denom = log_sum_exp(&log_probs);
        let mut mol_posteriors: Vec<f32> = Vec::new();

        for log_prob in log_probs {
            mol_posteriors.push((log_prob - bayes_log_denom).exp());
        }
        
        posteriors.push(mol_posteriors);
    }
    posteriors
}

fn log_sum_exp(p: &Vec<f32>) -> f32 {
    let max_p: f32 = p.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
    let sum_rst: f32 = p.iter().map(|x| (x - max_p).exp()).sum();
    max_p + sum_rst.ln()
}

fn get_long_read_molecules(vcf: &mut bcf::IndexedReader, 
        position_to_index: &mut HashMap<usize, usize>,
        position_so_far: &mut usize) -> Vec<Vec<Allele>> {
    let mut molecules: HashMap<String, Vec<Allele>> = HashMap::new();
    println!("check");
    for rec in vcf.records() {
        let rec = rec.expect("couldnt unwrap record");
        let pos = rec.pos();
        let var_index = position_to_index.entry(pos as usize).or_insert(*position_so_far);
        *position_so_far += 1;
        println!("pos {}", pos);
        match rec.format(b"RM").string() {
            Ok(rec_format) => {
                for rec in rec_format.iter() {
                    let ref_mols = std::str::from_utf8(rec).expect(":").to_string();
                    for read in ref_mols.split(";") {
                        let mol = molecules.entry(read.to_string()).or_insert(Vec::new());
                        mol.push(Allele{index: *var_index, allele: false});
                    }
                }
            },
            Err(_) => println!("fail"),
        }
        match rec.format(b"AM").string() {
            Ok(rec_format) => {
                for rec in rec_format.iter() {
                    let alt_mols = std::str::from_utf8(rec).expect(":").to_string();
                    for read in alt_mols.split(";") {
                        let mol = molecules.entry(read.to_string()).or_insert(Vec::new());
                        mol.push(Allele{index: *var_index, allele: true});
                    }
                }
            },
            Err(_) => println!("fail"),
        }
    }
    println!("check2");
    let mut to_return: Vec<Vec<Allele>> = Vec::new();
    for (_read_name, alleles) in molecules.iter() {
        if alleles.len() < 2 { continue; }
        let mut mol: Vec<Allele> = Vec::new();
        for allele in alleles {
            mol.push(*allele);
        }
        to_return.push(mol);
    }
    to_return
}

#[derive(Clone,Copy)]
struct Allele {
    index:  usize,
    allele: bool, // alt is true, ref is false
}


fn inspect_vcf(vcf: &mut bcf::IndexedReader, data: &ThreadData) -> (usize, i64) {
    let chrom = vcf.header().name2rid(data.chrom.as_bytes()).expect("cant get chrom rid");
    vcf.fetch(chrom, 0, None).expect("could not fetch in vcf");
    let mut num = 0;
    let mut last_pos = 0;
    for r in vcf.records() {
        num += 1;
        let rec = r.expect("could not unwrap vcf record");
        last_pos = rec.pos();
    }
    println!("vcf has {} variants for chrom {}", num, chrom);
    (num, last_pos)
}

fn init_cluster_centers(num: usize, data: &ThreadData) -> Vec<Vec<f32>> {
    let seed = [data.seed; 32];
    let mut rng: StdRng = SeedableRng::from_seed(seed);
    let mut cluster_centers: Vec<Vec<f32>> = Vec::new();
    for k in 0..data.ploidy {
        let mut cluster_center: Vec<f32> = Vec::new();
        for v in 0..num {
            cluster_center.push(0.5);
        }
        cluster_center[0] = rng.gen::<f32>().min(0.98).max(0.02);
        cluster_centers.push(cluster_center);
    }
    cluster_centers
}


fn get_all_variant_assignments(data: &ThreadData) -> Result<(), Error> {
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
    let header_view =  vcf_reader.header();
    let mut new_header = bcf::header::Header::from_template(header_view);
    new_header.push_record(br#"##fileformat=VCFv4.2"#);
    new_header.push_record(br#"##FORMAT=<ID=AM,Number=1,Type=String,Description="alt molecules">"#);
    new_header.push_record(br#"##FORMAT=<ID=RM,Number=1,Type=String,Description="ref molecules">"#);

 
    let mut vcf_writer = bcf::Writer::from_path(data.vcf_out.to_string(), 
        &new_header, true, Format::Vcf)?;
    let chrom = vcf_reader.header().name2rid(data.chrom.as_bytes())?;
    vcf_reader.fetch(chrom, 0, None)?;  // skip to chromosome for this thread
    let mut total = 0;
    let mut hets = 0;
    for (i, _rec) in vcf_reader.records().enumerate() {
        total += 1;
        let rec = _rec?;
        let pos = rec.pos();
        let alleles = rec.alleles();
        let mut new_rec = vcf_writer.empty_record();
        copy_vcf_record(&mut new_rec, &rec);
        if alleles.len() > 2 {
            continue; // ignore multi allelic sites
        }
        let reference = std::str::from_utf8(alleles[0])?;
        let alternative = std::str::from_utf8(alleles[1])?;
        let rec_chrom = rec.rid().expect("could not unwrap vcf record id");

        let genotypes = rec.genotypes()?;
        let genotype = genotypes.get(0); // assume only 1 and get the first one
        if is_heterozygous(genotype) {
            hets += 1;
            get_variant_assignments(
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
                &mut vcf_writer,
                &mut new_rec
            );
        }
        if hets > 200 { break; } //TODO remove, for small example
    }
    println!("done, saw {} records of which {} were hets in chrom {}", total, hets, data.chrom);
    //let output = Command::new("bgzip").arg(&data.vcf_out.to_string()).
    //    output().expect("lkjfd");
    let command = format!("#!/bin/bash\nbgzip -c {} > {}.gz && tabix -p vcf {}.gz", data.vcf_out, data.vcf_out, data.vcf_out);
    
    let index_script = format!("{}/index_{}.sh", data.output, data.chrom);
    println!("{}", &index_script);
    let mut f = File::create(&index_script).expect("Unable to create file");
    f.write_all(command.as_bytes()).expect("Unable to write data");
    println!("{}", format!("chmod 777 {}", &index_script));
    Command::new("chmod").args(&["777",&index_script]).output().expect("chmod?");
    Command::new(&format!("./{}", &index_script)).output().expect("seriously?");
    Command::new("rm").arg(&index_script).output().expect("cant delete");
    //let output = Command::new("bcftools").
    //    args(&["view", &data.vcf_out, "-Oz", "-o", &format!("{}.gz", data.vcf_out)]).
    //    output().expect("lkjf");
    //let result = cmd1.wait().expect("could not run bgzip");
    //println!("{}", output.status.success());

    //let output = Command::new("tabix").args(&["-p", "vcf", &format!("{}.gz", data.vcf_out)]).
    //   output().expect("why");
    //let output = Command::new("bcftools").args(&["index", &format!("{}.gz", data.vcf_out)]).
    //    output().expect("lkkll");
    //let result = cmd2.wait().expect("could not tabix index vcf");
    //println!("{}", output.status);
    fs::File::create(data.vcf_out_done.to_string())?;
    Ok(())
}

fn copy_vcf_record(new_rec: &mut bcf::record::Record, rec: &bcf::record::Record) {
    new_rec.set_rid(rec.rid());
    new_rec.set_pos(rec.pos());
    new_rec.set_id(&rec.id());
    for filter in rec.filters() {
        new_rec.push_filter(&filter);
    }
    new_rec.set_alleles(&rec.alleles()).expect("could not write alleles to new record???");
    new_rec.set_qual(rec.qual());
    let header = rec.header();
    for header_record in header.header_records() {
        match header_record {
            bcf::header::HeaderRecord::Filter{key, values} => {},
            bcf::header::HeaderRecord::Info{key, values} => {
                let mut format = FORMAT{Id: "blah".to_string(), Type: FORMAT_TYPE::Integer};
                for (x,y) in values {
                    match x.as_str() {
                        "ID" => format.Id = y,
                        "Type" => match y.as_str() {
                            "Integer" => format.Type = FORMAT_TYPE::Integer,
                            "Float" => format.Type = FORMAT_TYPE::Float,
                            "String" => format.Type = FORMAT_TYPE::String,
                            "Char" => format.Type = FORMAT_TYPE::Char,
                            &_ => (),
                        },
                        &_ => (),
                    }
                }
                match format.Type {
                    FORMAT_TYPE::Integer => {
                        match rec.info(&format.Id.as_bytes()).integer() {
                            Ok(rec_format) => {
                                for thingy in rec_format.iter() {
                                    new_rec.push_info_integer(&format.Id.as_bytes(), thingy).expect("fail1");
                                }
                            },
                            Err(_) => (),
                        }
                    },
                    FORMAT_TYPE::Float => {
                        match rec.info(&format.Id.as_bytes()).float() {
                            Ok(rec_format) => {
                                for thingy in rec_format.iter() {
                                    new_rec.push_info_float(&format.Id.as_bytes(), thingy).expect("fail1");
                                }
                            },
                            Err(_) => (),
                        }
                    },
                    FORMAT_TYPE::String => {
                        match rec.info(&format.Id.as_bytes()).string() {
                            Ok(rec_format) => {
                                new_rec.push_info_string(&format.Id.as_bytes(), &rec_format.expect("blerg")).expect("fail1");
                            },
                            Err(_) => (),
                        }
                    },
                    FORMAT_TYPE::Char => {
                        match rec.info(&format.Id.as_bytes()).string() {
                            Ok(rec_format) => {
                                new_rec.push_info_string(&format.Id.as_bytes(), &rec_format.expect("blerg2")).expect("fail1");
                            },
                            Err(_) => (),
                        }
                    },
                }
            },
            bcf::header::HeaderRecord::Format{key, values} => {
                let mut format = FORMAT {Id: "blah".to_string(), Type: FORMAT_TYPE::Integer};
                for (x,y) in values {
                    match x.as_str() {
                        "ID" => format.Id = y,
                        "Type" => match y.as_str() {
                            "Integer" => format.Type = FORMAT_TYPE::Integer,
                            "Float" => format.Type = FORMAT_TYPE::Float,
                            "String" => format.Type = FORMAT_TYPE::String,
                            "Char" => format.Type = FORMAT_TYPE::Char,
                            &_ => (),
                        },
                        &_ => (),
                    }
                }
                match format.Type {
                    FORMAT_TYPE::Integer => {
                        match rec.format(&format.Id.as_bytes()).integer() {
                            Ok(rec_format) => {
                                for thingy in rec_format.iter() {
                                    new_rec.push_format_integer(&format.Id.as_bytes(), thingy).expect("noooooooooo");
                                }
                            },
                            Err(_) => (),
                        }
                    },
                    FORMAT_TYPE::Float => {
                        match rec.format(&format.Id.as_bytes()).float() {
                            Ok(rec_format) => {
                                for thingy in rec_format.iter() {
                                    new_rec.push_format_float(&format.Id.as_bytes(), thingy).expect("fail1");
                                }
                            },
                            Err(_) => (),
                        }
                    },
                    FORMAT_TYPE::String => {
                        if format.Id == "GT".to_string() {
                            let sample_count = header.sample_count();
                            for i in 0..sample_count {
                                let gt = rec.genotypes().expect("lkjlkj").get(i as usize);
                                new_rec.push_genotypes(&gt);
                            }
                            
                        } else {
                            match rec.format(&format.Id.as_bytes()).string() {
                                Ok(rec_format) => {
                                    new_rec.push_format_string(&format.Id.as_bytes(), &rec_format).expect("fail1");
                                },
                                Err(_) => (),
                            }
                        }
                        
                    },
                    FORMAT_TYPE::Char => {
                        match rec.format(&format.Id.as_bytes()).string() {
                            Ok(rec_format) => {
                                new_rec.push_format_string(&format.Id.as_bytes(), &rec_format).expect("fail1");
                            },
                            Err(_) => (),
                        }
                    },
                }
            },
            bcf::header::HeaderRecord::Contig{key, values} => {},
            bcf::header::HeaderRecord::Structured{key, values} => {},
            bcf::header::HeaderRecord::Generic{key, value} => {},
        }
    }
}

struct FORMAT {
    Id: String,
    Type: FORMAT_TYPE,
}

#[derive(Debug)]
enum FORMAT_TYPE {
    Integer, Float, String, Char
}



fn get_variant_assignments<'a> (
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
    vcf_writer: &mut bcf::Writer,
    vcf_record: &mut bcf::record::Record,

) {
    if (pos + window) as u64 > chrom_length {
        return;
    }
    match long_reads {
        Some(bam) => {
            let tid = bam.header().tid(chrom.as_bytes()).expect("cannot find chrom tid");
            bam.fetch((tid, pos as u32, (pos + 1) as u32)).expect("blah"); // skip to region of bam of this variant position
            let ref_start = (pos - window) as u64;
            let ref_end = (pos + window + ref_allele.len()) as u64;
            fasta.fetch(chrom, ref_start, ref_end).expect("fasta fetch failed");
            let mut ref_sequence: Vec<u8> = Vec::new();
            fasta.read(&mut ref_sequence).expect("failed to read fasta sequence");
            //if pos == 69505 { println!("yes, this one"); }
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
            
            
            let mut read_names_ref: Vec<String> = Vec::new();
            let mut read_names_alt: Vec<String> = Vec::new();
            for _rec in bam.records() {
                let rec = _rec.expect("cannot read bam record");
                //if pos == 69505 {println!("testing molecule");}
                if rec.mapq() < min_mapq {
                    continue
                }
                if rec.is_secondary() || rec.is_supplementary() {
                    continue;
                }
                let mut read_start: Option<usize> = None; // I am lazy and for some reason dont know how to do things, so this is my bad solution
                let mut read_end: usize = rec.seq_len();
                let mut min_bq = 93;
                let qual = rec.qual();
                //println!("ref_start {} ref_end {}", ref_start, ref_end);
                for pos_pair in rec.aligned_pairs() {
                    //println!("\t{},{}",pos_pair[0],pos_pair[1]);
                    if (pos_pair[1] as u64) >= ref_start && (pos_pair[1] as u64) < ref_end {
                        if pos_pair[1] as usize >= pos && pos_pair[1] as usize <= pos + ref_allele.len().max(alt_allele.len()) {
                            min_bq = min_bq.min(qual[pos_pair[0] as usize]);
                        }
                        read_end = pos_pair[0] as usize;
                        if read_start == None {
                            read_start = Some(pos_pair[0] as usize); // assign this value only the first time in this loop that outter if statement is true
                            // getting the position in the read that corresponds to the reference region of pos-window
                        }
                    } 
                }
                
                if min_bq < min_base_qual {
                    continue;
                }
                if read_start == None {
                    println!("what happened, read start {:?} read end {}", read_start, read_end);
                    continue;
                }
                let read_start = read_start.expect("why read start is none");
                let seq = rec.seq().as_bytes()[read_start..read_end].to_vec();
                let score = |a: u8, b: u8| if a == b { MATCH } else { MISMATCH };
                let mut aligner = banded::Aligner::new(GAP_OPEN, GAP_EXTEND, score, K, W);
                let ref_alignment = aligner.local(&seq, &ref_sequence);
                let alt_alignment = aligner.local(&seq, &alt_sequence);
                if ref_alignment.score > alt_alignment.score {
                    read_names_ref.push(std::str::from_utf8(rec.qname()).expect("wtff").to_string());
                } else if alt_alignment.score > ref_alignment.score {
                    read_names_alt.push(std::str::from_utf8(rec.qname()).expect("wtf").to_string());
                } 
                //if pos == 69505 { println!("scores {} and {}",ref_alignment.score, alt_alignment.score);
                //println!("read_names_ref {}, read_names_alt {}", read_names_ref.len(), read_names_alt.len()); }
            }
            
            let concat_ref = read_names_ref.join(";");
            let concat_alt = read_names_alt.join(";");
            //if pos == 69505 { println!("{} {}", concat_ref, concat_alt); }
            vcf_record.push_format_string(b"AM", &[concat_alt.as_bytes()]).expect("blarg");
            vcf_record.push_format_string(b"RM", &[concat_ref.as_bytes()]).expect("gggg");
            vcf_writer.write(vcf_record).expect("nope");
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
    output: String,
    vcf_out: String,
    vcf_out_done: String,
    phasing_window: usize,
    seed: u8,
    ploidy: usize,
}

#[derive(Clone)]
struct Params {
    linked_read_bam: Option<String>,
    hic_bam: Option<String>,
    long_read_bam: Option<String>,
    output: String,
    seed: u8,
    threads: usize,
    fasta: String,
    ploidy: usize,
    vcf: String,
    min_mapq: u8,
    min_base_qual: u8,
    window: usize,
    phasing_window: usize,
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

    let threads = params.value_of("threads").unwrap();
    let threads = threads.to_string().parse::<usize>().unwrap();

    let seed = params.value_of("seed").unwrap(); 
    let seed = seed.to_string().parse::<u8>().unwrap();

    let min_mapq = params.value_of("min_mapq").unwrap();
    let min_mapq = min_mapq.to_string().parse::<u8>().unwrap();
    
    let min_base_qual = params.value_of("min_base_qual").unwrap();
    let min_base_qual = min_base_qual.to_string().parse::<u8>().unwrap();

    let fasta = params.value_of("fasta").unwrap();

    let ploidy = params.value_of("ploidy").unwrap();
    let ploidy = ploidy.to_string().parse::<usize>().unwrap();

    let allele_alignment_window = params.value_of("allele_alignment_window").unwrap();
    let allele_alignment_window = allele_alignment_window.to_string().parse::<usize>().unwrap();

    let phasing_window = params.value_of("phasing_window").unwrap();
    let phasing_window = phasing_window.to_string().parse::<usize>().unwrap();


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
        min_mapq: min_mapq,
        min_base_qual: min_base_qual,
        window: allele_alignment_window,
        phasing_window: phasing_window,
    }
}
