# HumanReadRemover
If you want to submit your raw SARS-CoV-2 sequencing reads to public data repositories (e.g. to the NUM data hub [COGDatm](https://cogdat.dat)), you can use this small script to identify and filter out reads of potentially human origin.

## Usage
Both Nanopore and paired-end Illumina data are supported. The data input type is determined automatically based on the supplied input arguments (`--inputFASTQ` for Nanopore data, `--inputFASTQ_R1` and `--inputFASTQ_R2` for paired-end Illumina data).

Only reads that
- generate no alignments to the human genome
- at least one alignment to the SARS-CoV-2 genome (but see below); for paired-end reads, one member read aligning to the SARS-CoV-2 genome is sufficient to satisfy this criterion

are written to the output FASTQ files (`--outputFASTQ` for Nanopore data, `--outputFASTQ_R1` and `--outputFASTQ_R2` for paired-end Illumina data)

### Outputting reads that do not generate alignments to the SARS-CoV-2 genome
If you want to output reads that do not align to the SARS-CoV-2 genome (i.e. if you only want to filter out reads aligning to the human genome), you can specify the parameter `--requireAlignment2SARSCov2 0`.

### Nanopore data
`perl filterHumanGenome.pl --inputFASTQ /path/to/input.fastq --outputFASTQ /path/to/output.fastq`

### Illumina data
`perl filterHumanGenome.pl --inputFASTQ_R1 /path/to/input_R1.fastq --inputFASTQ_R2 /path/to/input_R2.fastq --outputFASTQ_R1 /path/to/output_R1.fastq --outputFASTQ_R2 /path/to/output_R2.fastq`

## Installation
1. Clone this repository to a directory of your choice, e.g. by
```
git clone https://github.com/NUM-GenSurv/HumanReadRemover.git
```
2. Download the combined human + SARS-CoV-2 reference genome from OSF, extract it into the directory just created by the cloning process, unzip the downloaded file, and check that the md5sum is equal to `774032a66a7ac8db71f89b9f0c426b7f`:

```
cd HumanReadRemover
wget -O GRCh38_full_analysis_set_plus_decoy_hla_and_SARS-CoV-2.fa.gz https://osf.io/gtbre/download 
gunzip GRCh38_full_analysis_set_plus_decoy_hla_and_SARS-CoV-2.fa.gz
md5sum GRCh38_full_analysis_set_plus_decoy_hla_and_SARS-CoV-2.fa # check that output is equal to 774032a66a7ac8db71f89b9f0c426b7f
```

3. Make sure that `filterHumanGenome.pl` can find the `minimap2` and `samtools` binaries (see next step).

### Path setup
`filterHumanGenome.pl` makes use of [minimap2](https://github.com/lh3/minimap2) and [samtools](http://www.htslib.org/). In order to call the corresponding executables, `filterHumanGenome.pl` requires the paths to the corresponding executables. There are three ways to specify these:
- via the PATH: `filterHumanGenome.pl` will examine your PATH for the presence of these executables
- via the file `paths.ini`: you can specify the paths to these executables in the file `paths.ini`, located in the same directory as `filterHumanGenome.pl`
- via the command line: when calling, you can specify the parameters `--samtools_bin` and `--minimap2_bin`.

Command-line specified paths take precedence over paths specified in the file `paths.ini`, and paths specified in `paths.ini` take precendence over paths found in the PATH.




