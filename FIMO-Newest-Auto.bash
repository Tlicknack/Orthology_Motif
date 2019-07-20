#!/bin/bash
#PBS -k o
#PBS -M tlicknac@asu.edu
#PBS -m abe
#PBS -l nodes=2:ppn=4,walltime=08:00:00,vmem=10gb

cd /N/dc2/scratch/tlicknac/Paramecium_-200to0/Newest_All-Aurelias_FIMO_Results/Ptet/

input_dir="/N/dc2/scratch/tlicknac/Paramecium_-200to0/Newest_All-Aurelias_MEME_Results-Markov"

for filename in $input_dir/*; do 
        if [[ $filename == *[0-9][0-9][3-9][0-9][0-9]* ]]; then
                file=${filename}/meme.txt
                out_dir=$(basename "$filename")
                mkdir $out_dir
                cd $out_dir
                /N/u/tlicknac/Carbonate/meme/bin/fimo $file ~/Paramecium_FASTA/ptetraurelia_mac_51.fa
                cd ..
        fi
done
