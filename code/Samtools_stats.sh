#!/bin/bash
#PBS -l pmem=2700mb,nodes=1:ppn=3,walltime=3:00:00
#PBS -A brandvai
#PBS -m abe
#PBS -o /home/brandvai/pmonnaha/OandE/Samtools_stats.o
#PBS -e /home/brandvai/pmonnaha/OandE/Samtools_stats.e
#PBS -M pmonnaha@umn.edu
#PBS -q mesabi

module load parallel
module load liblzma


cd /panfs/roc/scratch/pmonnaha/Mimulus/BAMs/tmp

for f in *.bam; do echo "/panfs/roc/msisoft/samtools/1.3_gcc-4.9.2_haswell/bin/samtools stats -c 0,1000,1 "$f" > "${f}.stats"";done | parallel --jobs 3 --sshloginfile ${PBS_NODEFILE} --workdir ${PWD}
