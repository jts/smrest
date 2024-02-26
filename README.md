# smrest

smrest is a prototype somatic mutation caller for single molecule long reads. It uses haplotype phasing patterns for tumour samples that have a sigificant proportion of normal cells (purity > 0.3, < 0.8) to identify somatic mutations. For more details, see the preprint linked below.

## Citation

[Simpson, J.T., Detecting Somatic Mutations Without Matched Normal Samples Using Long Reads, BioRxiv](TODO)

## Compiling

This program is written in Rust and uses the [Cargo](https://doc.rust-lang.org/cargo/getting-started/installation.html) build system. After you have installed Cargo, you can compile this software from github as follows:

```
git clone https://github.com/jts/smrest.git
cd smrest
cargo build --release
```

## Usage

A Snakemake pipeline is provided in `workflow/Snakemake` to automate the process of calling heterozygous SNPs, phasing them with `whatshap`, then identifying somatic mutations.

TODO: demo data

## License

MIT

## Acknowledgements

This program reuses code originally developed by Edge et al for the [Longshot](https://github.com/pjedge/longshot) variant caller.
