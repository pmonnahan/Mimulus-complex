# Mimulus complex phylogeography

## Data

Data for this project came from Josh Puzey, Jenn Coughlan, and Lila Fishman.  JC provided 11 M. decorus samples, labelled JC-S followed by 1 - 11.  LF provided one sample of a micranthus-looking putatively guttatus individual, labelled YLRM1-3.1.  All other samples were provided by JP.

## Fastq Pre-processing

Samples with multiple fastqs were concatenated such that each sample ultimately had just two fastq files: one forward and one reverse reads file.  Used TrimGalore (wrapper for CutAdapt that automatically detects adapter type) for adapter trimming and quality trimming of reads.  

A few samples (DPP2, MCN, MED84, and SCH) had multiple sets of fastq's corresponding to whether the files used sanger quality encoding or not.  Used the samples with sanger encodings.

MED84 did not have equal reads across two paired-end files so TrimGalore failed.
Synced reads across the two files using perl script found online at https://github.com/sestaton/Pairfq:
```
curl -sL git.io/pairfq_lite > pairfq_lite

chmod +x pairfq_lite

./pairfq_lite makepairs \
-f MED84.R1.sanger.fastq \
-r MED84.R2.sanger.fastq \
-fp MED84.R1.sanger.sync.fastq \
-rp MED84.R2.sanger.sync.fastq \
-fs MED84.R1.sanger.sing.fastq \
-rs MED84.R2.sanger.sing.fastq
```

Moved all pre-processed fastq files to Tier 2 storage at Mimulus/fastq/CutAdapt/

## Mapping

Aligned all samples to Mguttatus_256_v2.0.hardmasked.fa.  The script MapReads2.py was used to submit mapping jobs to MSI cluster.  This script creates a shell script which uses bwa to map, samtools to index and calculate stats, and samblaster to remove duplicate reads.

All BAMs are stored in Tier2 Mimulus/BAMs/BAMs5-18

## Calculating stats

Calculated depth genome-wide as well as restricted to gene space using 'samtools depth' as implemented in Samtools_depth.sh.  As opposed to 'samtools stats', 'samtools depth' includes sites with 0 depth.  This will likely include hardmasked sites as well as other unalignable sites, so coverage values are likely an underestimate.  

Percent duplicates were obtained via samblaster, and percent of reads mapping in proper pairs was obtained via samtools flagstats.

This info is reported in sample_stats_2.csv

## Variant Calling

Used GATK's HaplotypeCaller followed by GenotypeGVCFs for variant calling, while retaining non-variant sites.  For several samples, HaplotypeCaller was throwing the error:
```
##### ERROR MESSAGE: Somehow the requested coordinate is not covered by the read
```

This error message would point to reads that were heavily hard-clipped.  Tried several read filter (-rf) options and ultimately found the following solution:
```
-rf OverclippedRead --filter_is_too_short_value 50
```

Used the scripts GATK_HaplotypeCaller.sh and GATK_GenotypeGVCFs.sh to generate vcfs and GATK_SelectVariants.sh to filter sites.  As a first pass, the filters used were:
```
	-select "QD > 2.0" \
	-select "MQ > 40.0" \
	-select "FS < 60.0" \
	-select "SOR < 3.0" \
	-select "DP < 1500" \
	-select "DP > 50" \
	-selectType SNP -selectType MNP \
    -restrictAllelesTo BIALLELIC \
```

The depth limits were based on a simple visualization script using a Jupyter notebook (Visualize_VCF_Cutoffs.ipynb). Chose values that produced a nice normal-ish distribution of read depth.

VCFs and gVCFs are all stored on Tier 2 at Mimulus/gVCFs and Mimulus/vcfs.

## Basic Visualizations

Wrote an RMarkdown script (VCF_BasicPlots.Rmd) that will produce basic visualizations such as PCA's, heatmaps of pairwise genetic distance, and a neighbor joining tree.  

