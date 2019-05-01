#!/bin/bash
#PBS -k o
#PBS -M tlicknac@asu.edu
#PBS -m abe
#PBS -l nodes=1:ppn=4,walltime=0:10:00,vmem=5gb

#Plan: Create background file for every fasta file of upstream regions

input_dir="/N/dc2/scratch/tlicknac/Newest_All-Aurelias_Upstream"        

for filename in $input_dir/*.fasta; do
	mv $filename ../
	var=`basename $filename | sed -r 's/_promoter.fasta//'`
	var2="${var}_background.fasta"
	var3=`basename $var2`
	cat *.fasta > "$var3"
	mv "../$(basename ${filename})" ./
	mv $var3 ../test_Background/
done
