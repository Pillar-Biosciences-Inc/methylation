# Pillar's methylation analysis script
An open-source, Bismark based pipeline for processing methylation sequencing data.

## Installation
We recommend a conda-based installation as it is faster to setup.  
1. Download and install conda in your local environment  
2. Optional - install mamba (more stable with environment resolutions)  
`conda install -c conda-forge mamba`
3. Create a new environment and install the dependencies  
```
conda create -n methylation  
conda activate methylation  
mamba install -c conda-forge -c bioconda bedtools bismark parallel pandas samtools trim-galore xlsxwriter
```  

4. Confirm all the programs are discoverable (this step doesn't gurantee that they work, only that they are installed)
```
which bedtools bismark bismark_methylation_extractor parallel samtools trim_galore
```

## Prepare the genomes
The pipeline requires both C->T (forward strand amplicons) and G->A (reverse strand amplicons) converted genome references available.  The following set of lines shows how one can download and create these references for this pipeline:
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

## Running the pipeline
Assuming the pipeline is discoverable/in a $PATH directory, the following is a typical pipeline execution:
```
pillar_methylation.py -i <input_directory> -o <output_directory> -t <parallel_threads> -g <reference_directory_location> -n <run_name>
```

Optionally, the user can also supply a BED file of input regions that they are interested in quantifying methylation levels for:
```
pillar_methylation.py -i <input_directory> -o <output_directory> -t <parallel_threads> -g <reference_directory_location> -n <run_name> -b <input_bed_file>
```
