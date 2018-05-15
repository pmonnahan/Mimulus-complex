#!/bin/bash
#PBS -l mem=96gb,nodes=1:ppn=24,walltime=06:00:00
#PBS -A brandvai
#PBS -m abe
#PBS -M pmonnaha@umn.edu
#PBS -o /home/brandvai/pmonnaha/OandE/GATK_HaplotypeCaller.o
#PBS -e /home/brandvai/pmonnaha/OandE/GATK_HaplotypeCaller.e
#PBS -q mesabi

module load java
module load gatk/3.7.0

newgrp brandvai

# Set the reference
REF="/home/brandvai/pmonnaha/mimulus_data/Mguttatus_256_v2.0.hardmasked.fa"
# Directory to BAM files
ALN_DIR="/panfs/roc/scratch/pmonnaha/Mimulus/BAMs/"
# Build the sample list from the files in this directory
SAMPLE_LIST=($(find ${ALN_DIR} -name '*.bam' | sort -V))
# Get the current sample, with the PBS_ARRAYID
CURRENT_SAMPLE=${SAMPLE_LIST[${PBS_ARRAYID}]}
# Get the sample name
SAMPLENAME=$(basename ${CURRENT_SAMPLE} | cut -f 1 -d '.')
# Set the output directory
OUTDIR="/panfs/roc/scratch/pmonnaha/Mimulus/gVCFs"

# What to do for heterozygosity?? J Kelly found synonymous diversity to be 0.03 in guttatus, but these are much more diverse samples.
mkdir -p ${OUTDIR}
export _JAVA_OPTIONS="-Xmx63g"
java -jar /panfs/roc/msisoft/gatk/3.7.0/GenomeAnalysisTK.jar\
    -T HaplotypeCaller\
    -R ${REF}\
    -I ${CURRENT_SAMPLE}\
    -o ${OUTDIR}/${SAMPLENAME}.g.vcf\
    -nct 24\
    --genotyping_mode DISCOVERY\
    --heterozygosity 0.05\
    --emitRefConfidence GVCF\
    -variant_index_type LINEAR\
    -variant_index_parameter 128000
    --activityProfileOut ${OUTDIR}/${SAMPLENAME}.activityProfile
