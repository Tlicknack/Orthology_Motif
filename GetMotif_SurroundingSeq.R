#DONE
#This script will: 
  #look for DNA motifs within the intergenic region of a genome of interest 
  #extract all upstream nt's to the start codon and downstream until the whole sequence == 150nts
  #create a fasta file of upstream-motif sequences for all instances of that motif
#INPUT: 
  #Fasta file of motifs with some ID
  #Fasta file of intergenic sequence from some organism of interest .... see Make_Intergenic_FASTA.R
  #Gff-tab file of intergenic DNA for species X .... See R script to produce this file
    #Layout- scafold_number   intergenic_start    intergenic_end    geneID_5+     geneID_3+     geneID_5-     geneID_3-
#OUTPUT:
  #Fasta file for each motif of all upstream sequences

#Need to do:
  #Get distance to start codon
    #Do this by returning value of end position - length of entire intergenic
  #Fix vpos problem
  #Get orientation ... low priority

#source("http://bioconductor.org/biocLite.R")
#biocLite("seqinr")
library("seqinr")

getIntergenic = function(fasta, vpos, size){
  vsequences = c()                                                        #initiate vector of sequences
  for(l in 1:length(vpos)){                                               #iterate through vector of positions of motif
    sequence = substr(fasta[1], (vpos[l]-(size/2)), (vpos[l]+(size/2)))      #extract surrounding region of motif
    sequence_annot = paste("[5']", sequence, "[3']", sep="")
    vsequences = append(vsequences, sequence_annot)                             #add it to vector of sequences
  }
  return(vsequences)
}


strandOrient = function(motif, position, seq){
  v_orient = c()
  reverse_motif = c2s(rev(comp(s2c(motif))))      #get reverse complementary sequence... convert string to list of char, then do comp and reverse, then convert back to string
  for(x in 1:length(position)){
    len_motif = (nchar(motif)-1)
    forward = substr(seq, position[x], (position[x]+len_motif))

    if(forward == motif){
      v_orient = append(v_orient, "+")
    }else {
      v_orient = append(v_orient, "-")
    }
  }
  return(v_orient)
}


findMotifs <- function(motif, seq){                                                                                         #Grep to find motif. Switch if interest in PWM or Dinucleotide Matrix
  vpos = c()                    #initiate vector of position
  res = gregexpr(motif, seq)    #get coordinates of regular expression match ... will match multiple times per sequence
  v = unlist(res)
  if(v[1] > (-1)){          #gregexpr returns -1 if it doesnt find motif in seq
    vpos = v                        #vector of positions
  }
  reverse_motif = c2s(rev(comp(s2c(motif))))      #get reverse complementary sequence... convert string to list of char, then do comp and reverse, then convert back to string
  rev_res = gregexpr(reverse_motif, seq)
  rev_v = unlist(rev_res)
  if(rev_v[1] > (-1)){
    vpos = append(vpos, rev_v)
  }
  return(vpos)
}


findIntergenic = function(intergenic_fasta, lmotifs, intergenic_size){ #take as an input the intergenic fasta from makeIntergenicFasta, and a list of motifs uploaded with seqinr's read.alignment
  total_counter = 0
  x = c("Motif", "Surrounding_Sequence(>100bp)", "Surrounding_Genes", "Position(Intergenic_FASTA)", "Distance(to ATG)", "Strand")              #Column names of final data frame
  final_seq_file = matrix(ncol=6)
  colnames(final_seq_file) = x
  for(i in 1:lmotifs[[1]]){                                             #iterate through fasta of motifs [[1]] which is number of sequences
    seq_counter = 0                                             #init a counter of sequences for each motif
    motif = as.character(lmotifs[[3]][i])                                             #sequence of i'th motif
    seq_count = length(intergenic_fasta)                                #get length of intergenic fasta file
    for(sequ in 1:seq_count){                                           #interate through length of sequences in intergenic fasta
      inter_seq = as.character(intergenic_fasta[[sequ]][1])                 #save the sequence of the current intergenic region
      vPosOfMotifs = findMotifs(motif, inter_seq)                         #RETURNS VECTOR OF POSITION OF MOTIFS

      if(length(vPosOfMotifs) > 0){
        v_extracted_seq = getIntergenic(intergenic_fasta[[sequ]], vPosOfMotifs, intergenic_size)          #RETURNS VECTOR OF SEQUENCES AROUND MOTIF
        strand = strandOrient(motif, vPosOfMotifs, inter_seq)                                             #RETURNS + or - FOR STRAND ORIENTATION
      } else{
        v_extracted_seq = NULL
      }
      if(length(v_extracted_seq > 1)){                      #Only add to final table if extraction worked. without this, you get error
          for(ex_seq in 1:length(v_extracted_seq)){         #iterate through number of extracted sequences (should be equal to vpos)
            seq_counter = seq_counter + 1                   #add 1 to the sequence counter so each new sequence is added to the next row in the temp table
            
            inter_length = nchar(inter_seq)
            len_to_motif = (vPosOfMotifs[ex_seq]+nchar(motif))
            seq_to_start = substr(inter_seq, len_to_motif, inter_length)
            dist_to_start = nchar(seq_to_start)
            
            #calculate distance to start codon within this loop to avoid making another function and returning vector 
            final_seq_file[(seq_counter+total_counter),1] = motif                                       #add values to the final output table 
            final_seq_file[(seq_counter+total_counter),2] = v_extracted_seq[ex_seq]
            final_seq_file[(seq_counter+total_counter),3] = getAnnot(intergenic_fasta[[sequ]])        #Modify this to simplify. 5+, 5-, 3+, 3- in separate columns?
            final_seq_file[(seq_counter+total_counter),4] = vPosOfMotifs[ex_seq]
            final_seq_file[(seq_counter+total_counter),5] = dist_to_start
            final_seq_file[(seq_counter+total_counter),6] = strand[ex_seq]
            tmp_row = c(1,2,3,4, 5, 6)                                                                    #create a temp row 
            final_seq_file = rbind(final_seq_file, tmp_row)                                         #bind the temp row to the existing table ... without this, we get out of bounds script
        }
      }         
    }
    total_counter = total_counter + seq_counter
  }
  #exportSeq(final_seq_file)
  return(final_seq_file)
}


#MAIN
lmotifs = read.alignment("/home/tlicknac/Desktop/Paramecium_Genome_Data/Paramecium_FASTA/Weibo_Motifs.fa", format="fasta", forceToLower = TRUE)
#lmotifs = read.alignment("/home/tlicknac/Desktop/Paramecium_Genome_Data/Jeff_motifs.fa", format="fasta", forceToLower = TRUE)               #fasta of motifs                        #CHANGE                                                                                           #CHANGE
#lmotifs = read.alignment("/home/tlicknac/Desktop/Paramecium_Genome_Data/tatabox.fa", format="fasta", forceToLower = TRUE) 
intergenic_size = 200
intergenic_fasta = read.fasta("/home/tlicknac/Desktop/Paramecium_Genome_Data/Paramecium_FASTA/ptet-intergenic.fasta", as.string=TRUE)
intergenic_seq = findIntergenic(intergenic_fasta, lmotifs, intergenic_size)                                                                   #save table as variable 
#write.table(intergenic_seq, file="Tata_Taattaa_surrounding_seqs.csv", sep=",")                                                                #make a file, HARDCODE
#write.table(intergenic_seq, file="Jeff_Motifs_surrounding_seqs.csv", sep=",")                                                                
write.table(intergenic_seq, file="Weibo_Motifs_surrounding_seqs.csv", sep=",")                                                                

