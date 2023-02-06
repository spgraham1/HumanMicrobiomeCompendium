import os
import sys

#the distnanceToNextSample calculation is time and computationally intenstive,
#even while running in parallel. To expedite further, the function takes
#only a single sample size to analyze at a time. I wrote this python script 
#to automatically generate slurm scripts to run the R function on all regions 
# and all sample sizes

if (len(sys.argv) < 5):
        print("Rerun script with arguments sampleSizes.txt, repetitions, regionList, and prefix (directory to write to)")
        quit()
filename = sys.argv[1]
reps = sys.argv[2]
#hours = sys.argv[3]
regionsFile = sys.argv[3]
dir = sys.argv[4]
f = open(filename,"r")
regions = open(regionsFile,"r")

mem = "64g"
threads = 32

if not os.path.exists(dir):
        os.mkdir(dir)


if not os.path.exists("results"):
        os.mkdir("results")

lines = f.readlines()
regionList = regions.readlines()

baseDir = "/home/blekhman/graha880/HumanMicrobiomeCompendium" + dir
#each line of this file contains a sample size that we want to query
for j in regionList:
    mem = "32g"
    threads = 32

    os.chdir(baseDir)
    j = j.rstrip()
    regionAbbr = list(filter(lambda c: c.isupper(), j))
    regionAbbr = ''.join(regionAbbr)

    regionDir = regionAbbr + "_d2s"
    if not os.path.exists(regionDir):
            os.mkdir(regionDir)

    os.chdir(regionDir)
    if not os.path.exists("results"):
            os.mkdir("results")

    for i in lines :
            #strip string of leading and trailing whitespace
            i = i.strip()
            #concatenate index with whatever we want the output file to be called

            output = f"d2s_{i}_{regionAbbr}.sh"
            if (int(i) * int(reps)) < 10000:
                    hours = 4
            elif (int(i) * int(reps)) < 1000000:
                    mem="48g"
                    hours = 8
            elif (int(i) * int(reps)) < 10000000:
                    hours = 14
                    mem = "64g"
            else:
                    hours = 20
                    mem = "72g"

            data = ""
            data += "#!/bin/bash -l"+f"\n#SBATCH --time={hours}:00:00"+f"\n#SBATCH --ntasks={threads}"+f"\n#SBATCH --mem={mem}"+"\n#SBATCH --mail-type=ALL"+"\n#SBATCH --mail-user=graha880@umn.edu"
            data += f"\n\n" + f"\ncd /home/blekhman/graha880/d2s/{regionDir}/results"
            data += "\nmodule load R"
            if threads == 64:
                    call = f"Rscript /home/blekhman/graha880/d2s/distanceToNextSample.r /home/blekhman/graha880/d2s/transformed_final_compendium.tsv  {i} {reps} {j} /home/blekhman/graha880/d2s/sample_metadata.tsv {threads}"
            else:
                    call = f"Rscript /home/blekhman/graha880/d2s/distanceToNextSample.r /home/blekhman/graha880/d2s/transformed_final_compendium.tsv  {i} {reps} {j} /home/blekhman/graha880/d2s/sample_metadata.tsv"
            data+=f"\necho \"{call}\""
            data+=f"\n{call}"
            out = open(output, 'w')
            out.write(data)


