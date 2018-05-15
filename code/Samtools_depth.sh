#!/bin/bash
#PBS -l pmem=2500mb,nodes=1:ppn=1,walltime=01:00:00
#PBS -A brandvai
#PBS -m abe
#PBS -o /home/brandvai/pmonnaha/OandE/Samtools_depth.o
#PBS -e /home/brandvai/pmonnaha/OandE/Samtools_depth.e
#PBS -M pmonnaha@umn.edu
#PBS -q mesabi

module load parallel
module load liblzma

cd /panfs/roc/scratch/pmonnaha/Mimulus/BAMs

export PARALLEL="--workdir . --env PATH --env LD_LIBRARY_PATH --env LOADEDMODULES --env _LMFILES_ --env MODULE_VERSION --env MODULEPATH --env MODULEVERSION_STACK --env MODULESHOME --env OMP_DYNAMICS --env OMP_MAX_ACTIVE_LEVELS --env OMP_NESTED --env OMP_NUM_THREADS --env OMP_SCHEDULE --env OMP_STACKSIZE --env OMP_THREAD_LIMIT --env OMP_WAIT_POLICY"

for f in *.bam; do echo "/panfs/roc/msisoft/samtools/1.7_gcc-7.2.0_haswell/bin/samtools depth -a "$f" > "${f}.all.depth" | awk '{sum+=\$3; sumsq+=\$3*\$3} END { print FILENAME,\"Average = \",sum/NR; print \"Stdev = \",sqrt(sumsq/NR - (sum/NR)**2)}'; /panfs/roc/msisoft/samtools/1.7_gcc-7.2.0_haswell/bin/samtools depth -b /home/brandvai/pmonnaha/mimulus_data/Mguttatus_256_v2.0.gene_only.bed "$f" > "${f}.genes.depth" | awk '{sum+=\$3; sumsq+=\$3*\$3} END { print FILENAME,\"Average = \",sum/NR; print \"Stdev = \",sqrt(sumsq/NR - (sum/NR)**2)}'";done | parallel --jobs 24 --sshloginfile ${PBS_NODEFILE} --workdir ${PWD}