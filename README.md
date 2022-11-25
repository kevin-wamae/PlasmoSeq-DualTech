
# A snakemake pipeline for variant calling from _P. falciparum_ short amplicon reads
### _The pipeline is still at it's infancy stage_



We sequenced an Illumina sequencing library on the Oxford Nanopore MinION (ONT) to evaluate the cost of this approach.
- PCR amplicons from _Plasmodium falciparum_ drug resistance markers (_ama1, k13, dhps, dhfr and mdr1_) were generated in duplicate.
- Illumina sequencing libraries were generated using KAPA reagents and KAPA indexes.
- Finally, ONT sequence libraries were generated using just one set of ONT adapters and sequenced on the ONT using the Flow Cell R9.4.1. So cannot demultiplex the sequences into individual samples, hence, further analyses were done at the population level.


### Below are the project dependencies:

#### &nbsp;&nbsp;&nbsp;&nbsp; Package management
- [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) - an open-sourcepackage management system and environment management system that runs on various platforms, including Windows, MacOS, Linux.

#### &nbsp;&nbsp;&nbsp;&nbsp; Workflow management
- [snakemake](https://anaconda.org/bioconda/snakemake) - a workflow management system that aims to reduce the complexity of creating workflows by providing a fast and comfortable execution environment, together with a clean and modern specification language in python style.


#### &nbsp;&nbsp;&nbsp;&nbsp; Bioinformatics tools (packages)
- [fastqc](https://anaconda.org/bioconda/fastqc) - a tool for a quality control tool for high throughput sequence data
- [multiqc](https://anaconda.org/bioconda/multiqc) - a tool for aggregating bioinformatics analysis reports across many samples and tools
- [porechop](https://anaconda.org/bioconda/porechop) - a tool for finding and removing adapters from Oxford Nanopore reads. Adapters on the ends of reads are trimmed off, and when a read has an adapter in its middle, it is treated as chimeric and chopped into separate reads. Porechop performs thorough alignments to effectively find adapters, even at low sequence identity.
- [cutadapt](https://anaconda.org/bioconda/cutadapt) - at tool that finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.
- [bwa](https://anaconda.org/bioconda/bwa) - an aligner for short-read alignment (see [_minimap2_](https://anaconda.org/bioconda/minimap2) for long-read alignment)
- [bedtools](https://anaconda.org/bioconda/bedtools) - allows one to intersect, merge, count, complement, and shuffle genomic intervals from multiple files in widely-used genomic file formats such as BAM, BED, GFF/GTF, VCF
- [bcftools](https://anaconda.org/bioconda/bcftools) - a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF.
- [snpEff](https://anaconda.org/bioconda/snpeff) - a genetic variant annotation and effect prediction toolbox
- [SnpSift](https://anaconda.org/bioconda/snpsift) - a toolbox that allows you to filter and manipulate annotated files.

### _Where to start_
- Clone this project into your computer with Git ([_installation instructions_](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)) using the following command:
  - `git clone https://github.com/kevin-wamae/ampSeq-short-read-ONT.git`

### Directory structure
- Below is the default directory structure:
    - **config/** - contains the workflow configuration
    - **env/**   - contains the Conda environment files
    - **input/** - contains fastq, adapters and genome files
    - **output/** - output from the analysis is stored here
    - **workflow/** - contains the Snakemake script (snakefile) and additonal scripts
```
.
├── config
│   └── config.yaml
├── env
│   └── environment.yml
├── input
│   ├── 01_fastq
│   ├── 02_adapters
│   │   ├── illumina-TruSeq-adapters.fasta
│   │   └── illumina-indexes.txt
│   └── 03_genome
│       ├── falciparum_3D7_v51_annotations.gff
│       └── falciparum_3D7_v51_genome.fasta
├── output
└── workflow
    ├── scripts
    │   └── create_snpeff_db.sh
    └── snakefile
```

#### Running the analysis
After installing [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html), executing the following commands:
- Create the conda analysis environment and install the dependencies from the **env/environment.yml** by running the following command in your terminal:
  - `conda env create env/environment.yml`
- Activate the conda environment. This needs to be done every time you want to execute this pipeline:
  - `conda activate ampseq-analysis`
- Create the `snpEff` database by executing the following **bash** script below. This will donwload P. falciparum genome files from PlasmoDB and create and a **snpEff** database:
  - `bash workflow/scripts/create_snpeff_db.sh`
- Finally, execute the `Snakemake` pipeline by running the following command in your terminal. Replace the digit **4** in the command with the number of CPUs you wish to use.
  - `snakemake -c4`
