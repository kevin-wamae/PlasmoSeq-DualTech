# A snakemake pipeline for variant calling from _P. falciparum_ short amplicon reads
### Motivation
_PS: The pipeline is still at it's infancy stage_



We sequenced an Illumina sequencing library on the Oxford Nanopore MinION (ONT) to evaluate the cost of this approach.
- PCR amplicons from _Plasmodium falciparum_ drug resistance markers (_ama1, k13, dhps, dhfr and mdr1_) were generated in duplicate.
- Illumina sequencing libraries were generated using KAPA reagents and KAPA indexes.
- Finally, ONT sequence libraries were generated using just one set of ONT adapters and sequenced on the ONT using the Flow Cell R9.4.1.
- Hence, we cannot demultiplex the sequences into individual samples and further analyses were done at the population level.

---

### Below are the project dependencies:

#### &nbsp;&nbsp;&nbsp;&nbsp; Package management
- [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) - an open-source package management system and environment management system that runs on various platforms, including Windows, MacOS, Linux.

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

---

### Where to start
- Clone this project into your computer using Git ([_installation instructions_](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)) with the following command:
  - `git clone https://github.com/kevin-wamae/ampSeq-short-read-ONT.git`
- Navigate into the cloned directory using the following command:
  - `cd ampSeq-short-read-ONT`
 
 ---

### Directory structure
- Below is the default directory structure:
    - **config/** - contains the workflow configuration files
    - **env/**   - contains the Conda environment files
    - **input/** - contains fastq, adaptors and genome files
    - **output/** - contains the output from the analysis
    - **workflow/** - contains the Snakemake script (snakefile) and additonal scripts
```
.
├── LICENSE
├── README.md
├── config
│ └── config.yaml
├── env
│ └── environment.yml
├── input
│ ├── 01_fastq
│ │ ├── file-1_0.fastq.gz
│ │ └── file-2_1.fastq.gz
│ ├── 02_adapters
│ │ ├── illumina-TruSeq-adapters.fasta
│ │ └── illumina-indexes.txt
│ └── 03_genome
│     ├── genome.fasta
│     └── genome_annotations.gff
├── output
└── workflow
    ├── scripts
    │ └── create_snpeff_db.sh
    └── snakefile
```

---

### Running the analysis
Install [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and execute the following commands:

1 - Create the conda analysis environment and install the dependencies from the ***env/environment.yml*** by running the following command in your terminal:
  - `conda env create --file env/environment.yml`
  
2 - Activate the conda environment:
  - _**PS** - This needs to be done every time you want to execute this pipeline_:
  - `conda activate ampseq-analysis`
  
3 - Create the `snpEff` database by executing the bash script below. This script will download *P. falciparum* genome files from PlasmoDB and create and a **snpEff** database:
  - _**PS** - for this analysis, we will use genome data **release-51** from [PlasmoDB](https://plasmodb.org/plasmo/app), and we only need to run it once_:
  - `bash workflow/scripts/create_snpeff_db.sh`

4 - Finally, execute the whole `Snakemake` pipeline by running the following command in your terminal:
  - _**PS** - Replace **4** in the command with the number of CPUs you wish to use_
  - `snakemake -c4`
  
5 - Alternatively, you can execute a specific rule by running the following command in your terminal:
  - _**PS** - Replace **rule** in the command with respective rule-name from the `workflow/Snakefile`_
  - `snakemake -c4 rule` (_for example_ `snakemake -c4 qc_raw_files`)
  
  ---
  
  ### Expected output
  Below is the expected directory structure of the **output/** directory:
  - **01_snpeff_database/** - contains the snpEff database for variant calling
  - **02_qc_raw/** - contains the fastqc QC reports from the raw fastq files
  - **03_multiqc_raw/** - contains the aggregated fastqc QC reports
  - **04_trim_fastq_ont/** - contains fastq files after trimming ONT adaptors
  - **05_trim_fastq_illumina/** - contains fastq files after trimming Illumina adaptors
  - **06_qc_trimmed_files/** - contains the fastqc QC reports from the fastq files after quality trimming
  - **07_read_mapping/** - contains genome mapping files (index, bam and bed)
  - **08_variant_calling/** - contains variant calling files
```
output/
├── 01_snpeff_database
│   ├── P.falciparum
│   └── genomes
├── 02_qc_raw
├── 03_multiqc_raw
│   └── multiqc_data
├── 04_trim_fastq_ont
├── 05_trim_fastq_illumina
├── 06_qc_trimmed_filesmed
├── 07_read_mapping
│   └── genomeIndex
└── 08_variant_calling
```
