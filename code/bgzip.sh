#!/bin/bash
#PBS -l pmem=2500mb,nodes=1:ppn=8,walltime=24:00:00
#PBS -A hirschc1
#PBS -m abe
#PBS -q mesabi
#PBS -o /home/brandvai/pmonnaha/OandE/bgzip.o
#PBS -e /home/brandvai/pmonnaha/OandE/bgzip.e

module load htslib

# MODIFY THESE PATHS BELOW AS WELL AS -o AND -e paths above
CMDS=($(find /panfs/roc/scratch/pmonnaha/Mimulus/vcfs -name "*.vcf" -print)) 

for c in "${CMDS[@]}"; do
	bgzip -@ 8 $c;
	tabix -p vcf "$c.gz"
done

