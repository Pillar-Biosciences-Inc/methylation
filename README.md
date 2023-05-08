# Pillar's methylation analysis script
An open-source, Bismark-based pipeline for processing methylation sequencing data. The pipeline is a wrapper around several open-source software and libraries, and is designed to work on Illumina generated FASTQ files. The output of the script is an excel file that summarizes several important output metrics from the pipeline.

## Installation
We recommend a conda-based installation as it is faster to setup.  
1. Download and install conda in your local environment  
2. Optional - install mamba (more stable with environment resolutions)  
`conda install -c conda-forge mamba`
3. Create a new environment and install the dependencies  
```
# Create the environment
conda create -n methylation --yes

# Activate the environment
conda activate methylation

# Install the required dependencies
mamba install -c conda-forge -c bioconda bedtools bismark parallel pandas samtools trim-galore xlsxwriter
```  

4. Confirm all the programs are discoverable (this step doesn't gurantee that they work, only that they are installed)
```
which bedtools bismark bismark_methylation_extractor parallel samtools trim_galore
```

## Prepare the genomes
Since our panel is designed to have amplicons with forward and reverse strand amplification (irrespective of the direction of the gene on the genome), we need to ensure that both forward (C->T) and reverse (G->A) strand bisulfite converted genomes references are available. The easiest way to create these references for the pipeline is to use bismark's genome preparation utility (`bismark_genome_preparation`). Feel free to supply N/2 threads for N available threads on your hardware using the `--parallel` option. The `bismark_genome_preparation` launches 2 instances of indexing and passes the input threads as is to both instances.
```
# Download the genomes
(seq 1 22; echo -e "X\nY\nM") | xargs -I CHROM -P 4 wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chrCHROM.fa.gz

# Concatenate the files and optionally remove the intermediate files
(seq 1 22; echo -e "X\nY\nM") | xargs -I CHROM sh -c "zcat chrCHROM.fa.gz > hg19.fa; rm chrCHROM.fa.gz"

# With the methylation environment activate, execute the following
bismark_genome_preparation --verbose ./
```

This will create the following directory structure:
```
 $ tree .
.
├── Bisulfite_Genome
│   ├── CT_conversion
│   │   ├── BS_CT.1.bt2
│   │   ├── BS_CT.2.bt2
│   │   ├── BS_CT.3.bt2
│   │   ├── BS_CT.4.bt2
│   │   ├── BS_CT.rev.1.bt2
│   │   ├── BS_CT.rev.2.bt2
│   │   └── genome_mfa.CT_conversion.fa
│   └── GA_conversion
│       ├── BS_GA.1.bt2
│       ├── BS_GA.2.bt2
│       ├── BS_GA.3.bt2
│       ├── BS_GA.4.bt2
│       ├── BS_GA.rev.1.bt2
│       ├── BS_GA.rev.2.bt2
│       └── genome_mfa.GA_conversion.fa
└── hg19.fa

3 directories, 15 files
```

**The directory where `Bisulfite_Genome` directory resides, is our *reference genome* directory**


**Please note that the directory where `Bisulfite_Genome` resides, is considered the *reference genome* directory by this script.**
E.g., if the full path to `Bisulfite_Genome` is `/home/user/methylation/reference/Bisulfite_Genome`, then use `/home/user/methylation/reference/` as the reference genome directory path as an input to the pipeline.

## Running the pipeline
Assuming that the pipeline is located in discoverable/$PATH directories and has the execute permission turned on, the following command syntax would launch a run:  
  
```
pillar_methylation.py -i <input_FASTQ_directories> -o <output_directory_name> -t <num_of_threads> -g <reference_genome_directory_path> -n <run_name>
```
  
Please note that the pipeline is designed to identify and pair Illumina named read files. If your file names have been altered and are no longer following Illumina's naming convention, the pipeline will fail to pair the samples correctly.  
  
Additionally, using the -b option, users can also supply a BED file containing regions of interest, for which they'd like to generate summary information via the pipeline:
```
pillar_methylation.py -i <input_FASTQ_directories> -o <output_directory_name> -t <num_of_threads> -g <reference_genome_directory_path> -n <run_name> -b <BED_file_path>
```
