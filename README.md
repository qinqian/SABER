# SABER
Selection of Amplicon Barcodes from Experimental Replicates or SABER is a NGS analysis pipeline based on and supplemental to the transgenic zebrafish system Genome Editing of Synthetic Target Arrays for Lineage Tracing (GESTALT) to study cellular phylogeny.

## Quick Start (TL;DR)
### Requirements
To run SABER you will need [snakemake](https://snakemake.readthedocs.io/en/stable/) and to be running macOS or linux. We recommend you install this via [conda](https://docs.conda.io/en/latest/miniconda.html). This has been tested on Ubuntu 18 and Ubuntu 16.
```
source ~/miniconda3/bin/activate
conda install -c bioconda -c conda-forge snakemake
```
Installing miniconda and then snakemake should take about five minutes. All other dependencies will be installed by snakemake in another conda environment.

### Running SABER
Clone or download a release. This should take 30 seconds to a minute.
```
git clone https://github.com/blachlylab/SABER.git
cd SABER
```
Run the demo.
```
# make sure you're in the SABER dir and can use snakemake
# download demo data and extract
wget https://github.com/blachlylab/SABER/releases/download/v1.0.0/demo_data.tar
mkdir demo
cd demo
tar -xf ../demo_data.tar
cd ..
# run pipeline
snakemake -j 8 --use-conda demo
```
First snakemake will spend 5-10 minutes downloading the conda environment required to run SABER (this is only done once). It will take another 5-15 minutes to run the demo set through the SABER pipeline and analysis. The SABER analysis scripts will also install some other R packages into the conda environment (this takes about 5 minutes and only needs to be done once).The conda environments are installed under .snakemake/conda/ and you could activate then by using:
```
conda activate .snakemake/conda/env-hash-here
```

On an Ubuntu 18 laptop with an 8-threaded intel cpu and 16GB of ram, installation of the conda environment took 7 minutes. Running the SABER pipeline and analysis (with installation of the R packages) took 12 minutes.
There will be a demo_output folder which contains all of the intermediate data used to perform the SABER analysis. The SABER folder under this contains the plots which should match the plots found in the folder demo/output in the current repository.

To setup your own data. Copy your fastqs into the fastq directory.
```
mkdir fastq
cp /path/to/your/fastqs* fastq
```
Run the pipeline on your data.
```
snakemake -j 8 --use-conda
```

### Configurable settings
SABER's YAML config file has some configuarable options.
```
output_prefix: SABER
use_umi: false
#variant_read_cutoff: 20000
#
# Uncomment below to set tvaf
# tvaf: 0.01
```
Output Prefix: prefixes the output directory

Use UMI: set if you used a UMI equipped PCR library

Variant Read Cutoff: When plotting the large heatmap of all samples all variants with more than 5000 supporting reads will be shown.

tVaf or Theta Vaf: variants accounting for more than 0.3% of reads in more than one sample are identified as common variants and excluded. If not provided, it will be calculated.


## Introduction
The transgenic zebrafish system, Genome Editing of Synthetic Target Arrays for Lineage Tracing (GESTALT), has emerged as a powerful tool for studying cellular phylogeny.<sup>1,2</sup>  GESTALT zebrafish carry a single germline copy of a synthetic array consisting of 10 tandem CRISPR/Cas9 target sites.  By microinjecting guide RNAs (gRNAs) targeting this array with either Cas9 mRNA or recombinant Cas9 protein into the single-cell zebrafish embryo, double strand breaks are induced within the array and then repaired by non-homologous end joining (NHEJ).  The combinatorial effect of editing the 10 sites of the array induces thousands of unique genetic barcodes.<sup>1</sup> The genetic barcoding can be used to trace cell phylogeny and in theory could be combined with conditional transgenics, mutants, or other genetic or chemical modifications to understand how these experimental conditions affect clonal diversity of the blood system. 

There are a number of limitations to the GESTALT method.  
*  The fraction of barcoded cells in the GESTALT zebrafish depends on the integrity and quantity of the reagents used and the efficiency with which the injection solution was delivered to the embryo. 
*  The NHEJ repair mechanism can produce stereotypical repair patterns, reducing the actual diversity of barcodes observed.  
*  The published bioinformatic analyses do not provide a means to identify these uninformative variants or a systematic approach to exclude samples with a low fraction of informative barcodes.  
*  Unique molecular identifiers (UMIs) have been used for sequencing error, PCR error, and PCR bias correction<sup>3</sup>, however the UMI PCR protocol more than doubles the sample preparation time and reagent cost, and the extent to which UMIs improve accuracy over standard PCR has not been demonstrated in this setting.  

By analyzing a large number of zebrafish blood samples, we were able to define thresholds for discriminating informative GESTALT barcodes from uninformative variants and then optimize these thresholds through modeling.  Using bootstrap analysis, we provide a method to rationally exclude samples with a low fraction of informative barcodes.  We also show that the results of a standard PCR protocol are nearly indistinguishable from a UMI-based PCR protocol.  We believe that this experimental method, conceptual framework and reference dataset will be useful to laboratories studying HSC clonal diversity and clonal evolution in the zebrafish system.     

## Pipeline
SABER is divided into three components:  1) core functions for processing sequence data, aligning to reference sequence and calling variants, 2) an optional module to handle UMI-based PCR amplicons and 3) functions to identify informative barcodes, select samples and perform statistical analysis (Fig. 1a). 

```
        SABER Pipeline

              fastqs
                  |
                  V
          Merge Read Pairs: pear
                  |
                  V
          Trim Reads: trimmomatic
                  |
                  V
          Filter reads: cutadapt
         /        |             \  If use_umi is true
        V         |              V
    fastqc        |            Extract UMI Sequences: umi_tools extract
       |          |                   |
       V          V                   V
    multiqc      Align to SABER/GESTALT Reference: needleall
                  |                      |  
                  |                      | If use_umi is true
                  |                      V
                  |        UMI Deduplication: umi_tools dedup
                  |                   |
                  V                   V
          SABER statistical analysis: analysis.R
```
#### Tools used
* [pear](https://cme.h-its.org/exelixis/web/software/pear/)
* [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
* [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [multiqc](https://multiqc.info/)
* [umi_tools](https://github.com/CGATOxford/UMI-tools)
* [needleall](http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needleall.html)
* [samtools](https://github.com/samtools/samtools)


## Statistical Analysis
In cellular phylogenetic terms, a GESTALT sequence variant can be considered an informative genetic barcode if it is unique to the clade of cells descended from the ancestral cell in which the barcoding was performed. We developed a method to search for variants shared between samples within the training set and classify them as informative and uninformative GESTALT sequences. After the variants are determined via CrispRVariants, we quantify the degree to which any two samples in the training set share GESTALT alleles. This Sharing Factor is calculated for each pair of samples in the training set. This value is optimized by filtering out uninformative variants and by optimizing the variant allele fraction threshold (tVaf or Theta Vaf). If a variant allele is detected above the threshold fraction in more than one sample it is considered a common variant and removed from further analysis.

SABER generates a model for minimizing inter-sample variant sharing and maximizing the fraction of informative reads (number of reads supporting informative barcodes / the total number of reads) remaining in each sample, the analysis was repeated with with tVaf = 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3 or 1.0.  Mean Sharing Factor and mean informative read fraction values are plotted for the 28 training samples analyzed under each condition (file name output).  The optimal tVaf value based on this permutation analysis in our data was tVaf =0.003 (mean Sharing Factor = 0.0065, mean informative read fraction = 0.71). At these threshold settings, variants accounting for more than 0.3% of reads in more than one sample are identified as common variants and excluded. Your optimal tVaf and informative read fraction may be different and these thresholds are configurable.

## Citations
1. McKenna, A. et al. Whole-organism lineage tracing by combinatorial and cumulative genome editing. Science 353, aaf7907, doi:papers3://publication/doi/10.1126/science.aaf7907 (2016).
2. Raj, B. et al. Simultaneous single-cell profiling of lineages and cell types in the vertebrate brain. Nat Biotechnol, doi:10.1038/nbt.4103 (2018).
3. Clement, K., Farouni, R., Bauer, D. E. & Pinello, L. AmpUMI: design and analysis of unique molecular identifiers for deep amplicon sequencing. Bioinformatics 34, i202-i210, doi:10.1093/bioinformatics/bty264 (2018).
