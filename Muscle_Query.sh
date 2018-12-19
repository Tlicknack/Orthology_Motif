#!/bin/bash

seq=/N/u/tlicknac/Carbonate/Paramecium_Upstream_Sequences/*

cd /N/u/tlicknac/Carbonate/Paramecium_Upstream_Alligned/

for fasta in $seq
do
        echo $fasta
        outfile="$(basename $fasta)"
        muscle -in $fasta -out $outfile
done
