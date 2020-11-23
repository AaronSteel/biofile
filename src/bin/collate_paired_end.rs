use biofile::bed::paired_end_collator::PairedEndCollator;
use clap::{clap_app, Arg};
use program_flow::{
    argparse::{extract_numeric_arg, extract_str_arg},
    eprint_named_vars, OrExit,
};

fn main() {
    let mut app = clap_app!(collate_paired_end =>
        (about: "Extracts paired-end reads and generate a single interval from \
        each pair of reads.")
    );
    app = app
        .arg(
            Arg::with_name("bed_path")
                .long("bed")
                .short("b")
                .takes_value(true)
                .required(true)
                .help("path to the bed file containing the paired-end reads"),
        )
        .arg(
            Arg::with_name("out_path")
                .long("out-path")
                .short("o")
                .takes_value(true)
                .required(true)
                .help("output file path."),
        )
        .arg(
            Arg::with_name("left_offset")
                .long("left")
                .short("l")
                .takes_value(true)
                .required(true)
                .long_help(
                    "Ater finding the mid-point M of each paired-end reads, \
                    the start coordinate for the single interval generated in \
                    the output will be M - right_offset. Note that in the BED \
                    format the start coordinate 0-based inclusive.",
                ),
        )
        .arg(
            Arg::with_name("right_offset")
                .long("right")
                .short("r")
                .takes_value(true)
                .required(true)
                .long_help(
                    "Ater finding the mid-point M of each paired-end reads, \
                    the end coordinate for the single interval generated in the \
                    output will be M + right_offset. Note that in the BED \
                    format the end coordinate is 0-based exclusive.",
                ),
        );
    let matches = app.get_matches();

    let bed_path = extract_str_arg(&matches, "bed_path");
    let out_path = extract_str_arg(&matches, "out_path");

    let left_offset = extract_numeric_arg(&matches, "left_offset")
        .unwrap_or_exit(Some("failed to extract the left offset argument"));

    let right_offset = extract_numeric_arg(&matches, "right_offset")
        .unwrap_or_exit(Some("failed to extract the right offset argument"));

    eprint_named_vars!(bed_path, out_path, left_offset, right_offset);

    let collator = PairedEndCollator::new(&bed_path);
    let result = collator
        .write(&out_path, left_offset, right_offset)
        .unwrap_or_exit(Some(
            "failed to collate the paired-end reads bed file",
        ));

    println!(
        "{}",
        format!("number of PCR duplicates: {}", result.num_pcr_duplicates)
    );
}
