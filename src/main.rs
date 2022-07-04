use rust_htslib::{bam, faidx, bam::Read, bam::record::Aux};
use clap::{App, SubCommand, Arg};

fn main() {
    let matches = App::new("smrest")
        .version("0.1")
        .author("Jared Simpson <jared.simpson@oicr.on.ca>")
        .about("Toolkit for estimating somatic mutation rate from long reads")
        .subcommand(SubCommand::with_name("sim")
                .about("generate a bam file containing simulated mutations")
                .arg(Arg::with_name("mutation-rate")
                    .short('r')
                    .long("mutation-rate")
                    .takes_value(true)
                    .help("the rate of mutations to generate, expressed in number of mutations per megabase"))
                .arg(Arg::with_name("input-bam")
                    .required(true)
                    .index(1)
                    .help("the input bam file to process")))
        .subcommand(SubCommand::with_name("extract")
                .about("gather candidate somatic mutations from aligned reads")
                .arg(Arg::with_name("genome")
                    .short('g')
                    .long("genome")
                    .takes_value(true)
                    .help("the reference genome"))
                .arg(Arg::with_name("input-bam")
                    .required(true)
                    .index(1)
                    .help("the input bam file to process"))).get_matches();
    
    if let Some(matches) = matches.subcommand_matches("extract") {
        extract_mutations(matches.value_of("input-bam").unwrap(),
                          matches.value_of("genome").unwrap());
    }
}

fn base2index(base: char) -> i32 {
    return match base {
        'A' => 0,
        'C' => 1,
        'G' => 2,
        'T' => 3,
        _ => -1
    };
}

fn haplotype_index(record: &bam::Record) -> Option<i32> {
    match record.aux(b"HP") {
        Ok(value) => {
            if let Aux::I32(v) = value {
                return Some(v)
            } else {
                return None
            }
        }
        Err(_e) => return None
    }
}

fn extract_mutations(input_bam: &str, reference_genome: &str) {

    let mut bam = bam::IndexedReader::from_path(input_bam).unwrap();
    println!("chromosome\tposition\treference_base\tvariant_base\tcanonical_type\tcanonical_context\taligned_depth\t\
              hmajor_variant_count\thminor_variant_count\thmajor_vaf\thminor_vaf\th1_a\th1_c\th1_t\th1_t\th2_a\th2_c\th2_g\th2_t");

    let chromosome_name = "chr20";
    let min_variant_observations = 5;

    // set up header and faidx for pulling reference sequence
    let faidx = faidx::Reader::from_path(reference_genome).expect("Could not read reference genome:");

    let header = bam::Header::from_template(bam.header());
    let header_view = bam::HeaderView::from_header(&header);
    let tid = header_view.tid(chromosome_name.as_bytes()).unwrap();
    let chromosome_length = header_view.target_len(tid).unwrap() as usize;
    let mut chromosome_sequence = faidx.fetch_seq_string(chromosome_name, 0, chromosome_length).unwrap();
    chromosome_sequence.make_ascii_uppercase();
    let chromosome_bytes = chromosome_sequence.as_bytes();

    // go to chromosome of interest
    bam.fetch(chromosome_name).unwrap();
    for p in bam.pileup() {
        let pileup = p.unwrap();
        
        let mut aligned_depth = 0;

        // four bases, on two haplotypes
        let mut base_counts: [u32; 8] = [0; 8];
        let mut haplotype_depth: [u32; 2] = [0; 2];

        for alignment in pileup.alignments() {
            if let Some(qpos) = alignment.qpos() {
                if let Some(haplotype_index) = haplotype_index(&alignment.record()) {
                    let read_base = alignment.record().seq()[qpos] as char;
                    let bi = base2index(read_base);
                    let bci = (haplotype_index - 1) * 4 + bi;
                    base_counts[bci as usize] += 1;
                    aligned_depth += 1;
                    haplotype_depth[ (haplotype_index - 1) as usize ] += 1;
                }
            }
        }

        let reference_position = pileup.pos() as usize;
        let reference_base = chromosome_bytes[reference_position] as char;
        let reference_base_index = base2index(reference_base) as usize;
        let mut max_variant_count = 0;
        let mut max_variant_index = 0;
        let mut max_haplotype_index = 0;

        for i in 0usize..8usize {
            if i != reference_base_index && i != reference_base_index + 4 && base_counts[i] > max_variant_count {
                max_variant_count = base_counts[i];
                max_variant_index = i % 4;
                max_haplotype_index = i / 4;
            }
        }

        if max_variant_count >= min_variant_observations && reference_base != 'N' {

            // grab reference context
            //let mut reference_context = Vec::new();
            //reference_context.clone_from_slice( &chromosome_bytes[ (reference_position - 1)..(reference_position + 2)] );
            let mut reference_context = chromosome_bytes[ (reference_position - 1)..(reference_position + 2)].to_vec();

            let bases = "ACGT";
            let variant_base = bases.as_bytes()[max_variant_index as usize] as char;
            
            // grab mutation type
            let mut mutation_type: [char; 3] = [ reference_base, '>', variant_base ];

            // convert mutation type/context to canonical form C>x, T>x
            if mutation_type[0] != 'C' && mutation_type[0] != 'T' {
                mutation_type[0] = bio::alphabets::dna::complement(mutation_type[0] as u8) as char;
                mutation_type[2] = bio::alphabets::dna::complement(mutation_type[2] as u8) as char;
                reference_context = bio::alphabets::dna::revcomp(reference_context);
            }
            
            let mutation_type_str = String::from_iter(&mutation_type);
            let context_str = std::str::from_utf8(&reference_context).unwrap();

            // read support statistics
            let position = pileup.pos() + 1; // to match vcf

            
            let hmaj_vaf = max_variant_count as f32 / haplotype_depth[max_haplotype_index] as f32;
            
            let minor_haplotype_index = 1 - max_haplotype_index;
            let minor_variant_count = base_counts[ minor_haplotype_index * 4 + max_variant_index];
            let hmin_vaf = minor_variant_count as f32 / haplotype_depth[minor_haplotype_index] as f32;
            
            println!("{chromosome_name}\t{position}\t{reference_base}\t{variant_base}\t{mutation_type_str}\t\
                      {context_str}\t{aligned_depth}\t{max_variant_count}\t{minor_variant_count}\t\
                      {hmaj_vaf:.3}\t{hmin_vaf:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                      base_counts[0], base_counts[1], base_counts[2], base_counts[3],
                      base_counts[4], base_counts[5], base_counts[6], base_counts[7]);
            // depth {} {} {} {} {} {}", 
            //    pileup.pos(), reference_base, aligned_depth, max_variant_index, base_counts[0], base_counts[1], base_counts[2], base_counts[3]);
        }
    }
}

