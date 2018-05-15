#!/usr/bin/env python3
# written in perl by Levi Yant, streamlined by Jeff DaCosta, translated to python3 and adjusted to lates GATK bests practises
# by Christian Sailer, 20 May 2016

import os, sys, argparse, subprocess

print()

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description='This script uses the fastq files that have the adapters trimmed off (output of step P02). '+
                                 'The script creates one SLURM shell scripts to map the short reads to a reference using bwa mem per input'+
                                 ' file and sends it to the NBI SLUM cluster. If the output directory does not exist, it is created '+
                                 'automatically.')

parser.add_argument('-fastqdir', type=str, metavar='fastq_ca', default='fastq_ca/', help='REQUIRED: Full path to directory with input fastq files')
parser.add_argument('-R', type=str, metavar='reference_path', default='/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta', help='REQUIRED: Full path to reference genome')
parser.add_argument('-o', type=str, metavar='aligned_dir', default='aligned/', help='Realtive path to output directory [aligned/]')
parser.add_argument('-oe', type=str, metavar='OandEdirectory', default='aligned/', help='Realtive path to output/error directory [aligned/]')
parser.add_argument('-c', type=int, metavar='num_cores', default='12', help='Number of requested cores for job [5]')
parser.add_argument('-mem', type=str, metavar='memory', default='16000', help='Total memory for each job (Mb) [16000]')
parser.add_argument('-time', type=str, metavar='time', default='2-00:00', help='Time for job [2-00:00]')
parser.add_argument('-print', type=str, metavar='print', default='false', help='If changed to true then shell files are printed to screen and not launched [false]')
args = parser.parse_args()


# gather list of input fastq files
in_fastq_list1 = []
in_fastq_list2 = []
for file in os.listdir(args.fastqdir):
    if file.endswith('.fq.gz'):
        if file.strip(".fq.gz")[-1] == "1":
            in_fastq_list1.append(file)
        elif file.strip(".fq.gz")[-1] == "2":
            in_fastq_list2.append(file)
in_fastq_list1.sort()
in_fastq_list2.sort()

print('\nFound '+str(len(in_fastq_list1))+' input fastq files')
for in_fastq in in_fastq_list1:
    print('\t'+in_fastq)
print('\n\nCreating shell files and sending jobs to NBI SLURM cluster\n\n')

# check if output directory exists, create it if necessary
if os.path.exists(args.o) is False:
    os.mkdir(args.o)


#loop through input fastq files
count = 0
for i, fq1 in enumerate(in_fastq_list1):
    fastq_basename = fq1.replace('.fq.gz', '')
    name = fastq_basename.split("_")[0]
    #write slurm shell file
    sh_file = open(args.o + fastq_basename + '.sh', 'w')
    sh_file.write('#!/bin/bash -e\n'+
                  '#PBS -l mem=12gb,nodes=1:ppn=4,walltime=24:00:00\n'+
                  '#PBS -A brandvai\n'+
                  '#PBS -q mesabi\n'+
                  '#PBS -o ' + args.oe + fastq_basename + '.out\n'+ 
                  '#PBS -e ' + args.oe + fastq_basename + '.err\n'+
                  'module load bwa\n'+
                  'module load samtools\n'+
                  r'bwa mem -t 4 -R "@RG\\tID:' + name + r'\\tSM:' + name + r'\\tLB:' + name + '" ' + args.R + ' ' + args.fastqdir + fq1 + ' ' + args.fastqdir + in_fastq_list2[i] + ' | samtools view -Sbh -@ 4 - | samtools sort -n -@ 4 -T ' + args.o + name + '.tmp > ' + args.o + name + '.TMP.bam\n' +
                  'samtools index ' + args.o + name + '.TMP.bam\n' +
                  'samtools flagstat ' + args.o + name + '.TMP.bam > ' + args.o + name + '.flagstats\n'
                  'samtools view -h ' + args.o + name + '.TMP.bam | /home/hirschc1/pmonnaha/software/speedseq/bin/samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 | samtools view -bS -f 0x2 -@ 4 - | samtools sort -@ 4 -T ' + args.o + name + ' > '+
                  args.o + name + '.sort.bam\n'+
                  'samtools index ' + args.o + name + '.sort.bam\n'+
                  'samtools idxstats ' + args.o + name + '.sort.bam > ' + args.o + fastq_basename + '.idxstats\n'+
                  'samtools flagstat ' + args.o + name + '.sort.bam > ' + args.o + name + '.flagstat\n'+
                  'printf "\\nFinished\\n\\n"\n')

    #check if slurm shell file should be printed or sent to NBI SLURM 
    if args.print == 'false':
        #send slurm job to NBI SLURM cluster
        cmd = ('qsub ' + args.o + fastq_basename + '.sh')
        p = subprocess.Popen(cmd, shell=True)
        sts = os.waitpid(p.pid, 0)[1]
    else:
        file = open(args.o+fastq_basename+'.sh','r')
        data = file.read()
        print(data)
    count +=1
    
#if appropriate, report how many slurm shell files were sent to NBI SLURM 
if args.print == 'false':
    print('\nSent '+str(count)+' jobs to the NBI SLURM  cluster\n\nFinished!!\n\n')
