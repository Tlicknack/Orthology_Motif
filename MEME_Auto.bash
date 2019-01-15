#!/bin/bash
#PBS -k o
#PBS -M tlicknac@asu.edu
#PBS -m abe
#PBS -l nodes=1:ppn=4,walltime=0:10:00,vmem=5gb

input_dir="/N/dc2/scratch/tlicknac/test_Upstream"	
	
for filename in $input_dir/*; do 
	output_file=$(basename "$filename"); 
	output_dir="/N/dc2/scratch/tlicknac/test_MEME_Results/"; 
	output="$output_dir""$output_file";  
	/N/u/tlicknac/Carbonate/software/meme/install/bin/meme $filename -dna -oc $output -nostatus -time 18000 -maxsize 60000  -nmotifs 5 -minw 6 -maxw 15
done
