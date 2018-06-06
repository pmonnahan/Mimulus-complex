import os, sys, argparse, subprocess

# create variables that can be entered in the command line
parser = argparse.ArgumentParser()
parser.add_argument('-dir', type=str, metavar='', required=True, help='parent directory containing fastq files')
parser.add_argument('-o', type=str, metavar='', required=True, help='text nameing the outfolder of concatenated files')
args = parser.parse_args()

# create objects to work with
name = []
folders = []
inlist = []
outlist = []
outfolder = [args.o]
baseoutname = []

file_dict={}
# populate name list
for dir, subdirs, files in os.walk(args.dir):
    for file in files:
        if file.endswith(".fastq.gz") and len(file.split("_")) > 1:  #
            if (file.split("_")[-1].split(".")[0] == "R1" or file.split("_")[-1].split(".")[0] == "R2"):
                name = file.split("_")[0]
                Read = file.split("_")[-1].split(".")[0]
                file_dict.setdefault(name, [[], []])
                file_dict[name][int(Read[1]) - 1].append(os.path.join(dir, file))



sh = open(args.o + 'concat_fastqs.txt', 'w') #file containing zcat commands to be run with gnu parallel
sh1 = open(args.o + 'fastq_list.txt', 'w') #file containing concatenated fastqs that can be used for TrimGalore 
for name, fastqs in file_dict.items():
    print(name)
    for i, read in enumerate(fastqs):
        sh.write('zcat ')
        for fastq in read:
            sh.write(fastq + " ")
        sh.write("> /panfs/roc/scratch/pmonnaha/Mimulus/fastq/" + name + "_concat_R" + str(i + 1) + ".fastq\n")
    sh1.write("/panfs/roc/scratch/pmonnaha/Mimulus/fastq/" + name + "_concat_R1.fastq /panfs/roc/scratch/pmonnaha/Mimulus/fastq/" + name + "_concat_R2.fastq /panfs/roc/scratch/pmonnaha/Mimulus/fastq/CutAdapt/" + name + "_concat_R1.CA.fastq /panfs/roc/scratch/pmonnaha/Mimulus/fastq/CutAdapt/" + name + "_concat_R2.CA.fastq\n")   
sh.close()
sh1.close()
