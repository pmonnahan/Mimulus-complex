#!/bin/bash
#PBS -l mem=22gb,nodes=1:ppn=8,walltime=12:00:00
#PBS -A brandvai
#PBS -o /home/brandvai/pmonnaha/OandE/GATK_SelectVariants.o
#PBS -e /home/brandvai/pmonnaha/OandE/GATK_SelectVariants.e
#PBS -q mesabi

module load java
module load gatk/3.7.0

#	Path to the reference
REF="/home/brandvai/pmonnaha/mimulus_data/Mguttatus_256_v2.0.hardmasked.fa"
REF_DICT="/home/brandvai/pmonnaha/mimulus_data/Mguttatus_256_v2.0.hardmasked.dict"
#	Path to directory containing VCFs
VCF_DIR="/panfs/roc/scratch/pmonnaha/Mimulus/vcfs/"
#	 Build the vcf list
VCF_LIST=($(find ${VCF_DIR} -maxdepth 1 -name '*.vcf'))

#   Get which sequence we are analyzing with the task ID
CURRENT=${VCF_LIST[${PBS_ARRAYID}]}
#	Make sure the output directory exists
OUTPUT_DIR="/panfs/roc/scratch/pmonnaha/Mimulus/vcfs/variant/"
mkdir -p ${OUTPUT_DIR}

SCAFF=$(basename ${CURRENT} | cut -d "_" -f 2,3)

export _JAVA_OPTIONS="-Xmx20g"
java -jar /panfs/roc/msisoft/gatk/3.7.0/GenomeAnalysisTK.jar\
	-T SelectVariants \
	-R ${REF} \
	-V ${CURRENT} \
	-o ${OUTPUT_DIR}/Mim_${SCAFF}_FULL_Variant_GATK.vcf \
	--excludeNonVariants