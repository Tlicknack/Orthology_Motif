#DONE

#This script will use aligned upstream FASTA's and count mismatches at each site, averaging them to the total number of sites:
  #i.e. mismatches per nt in orthologous sequences
  
library("seqinr")

#MAIN
input_list = list.files("/N/dc2/scratch/tlicknac/All_Aurelias_Upstream_Aligned/", "*.fasta", all.files=FALSE, recursive=FALSE)  #get all fasta files in a directory   
lmismatches = list()

for(i in 1:length(input_list)) {                                                      #iterate through fasta file names
  cat("Onto file number: ", i, "\n")
  seq_dist = 0
  path = "/N/dc2/scratch/tlicknac/All_Aurelias_Upstream_Aligned/"                      #save path to files as a string       
  path_i = paste(path, input_list[[i]], sep="")                                       #concatenate path and i'th filename
  input_fasta = read.alignment(path_i, format="fasta", forceToLower = TRUE)           #input the fasta files using alignment reader
  #[[1]] = number of seqs in the current fasta file,    [[2]] = list of geneIDs,                          [[3]] = list of sequence strings
  #[[1]] has 1 dimension,                               [[2]][i] allows you to go acress i'th geneID,     [[3]][i] allows you to access the i'th seq
  current_file = input_list[[i]]
  current_seqs = input_fasta[[3]]
  number_of_seqs = length(current_seqs)
  total_mismatches = 0
  len_alignment = nchar(current_seqs[[1]])
  vmismatches = c()
  
  for(pos in 1:len_alignment){
    vpos = c()                                                      #initiate vector of positions
    
    for(seq in 1:number_of_seqs){
      vpos = append(vpos, substr(current_seqs[seq], pos, pos))      #append vector of positions with current element
    }
    elements = sort(table(vpos), decreasing=T)                      #create a count for each element in vector... "-" "a" "t" "c" "g"
    most_common = names(elements[1])                                #return character that is most abundant ... CAN USE TO BUILD CONSENSUS SEQUENCE 
    number_of_mismatches = length(which(vpos != most_common))       #returns integer of positions that differ from consensus
    vmismatches = append(vmismatches, number_of_mismatches)         #append each integer to vector
    mismatch_per_site = sum(vmismatches) / len_alignment            #calculate mismatches per site
    
  }
  lmismatches[[i]] = c(input_list[[i]], mismatch_per_site)
}
final_df = t(data.frame(lmismatches))
colnames(final_df) = c("File", "Mismatch/Site")
write.csv(final_df, file="All_Aurelia_12ntMotif_MismatchPerSite.csv", sep="\t", row.names=F)
