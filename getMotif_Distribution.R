#This script will take a regular expression and scan across the 15 genomes, returning hits of where it was found



getPositionInfo = function(regex_positions, scaf, fasta){  
  lfeatures = list()
  
  if(fasta == "biaurelia_V1-4_assembly_v1.fa"){                                                #loading INTERGENIC GFF's for simpler searching
    gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/pbi-intergenic.tab", sep="\t", as.is=T, header=T)
  } 
  if(fasta == "ppentaurelia_mac_87_v1.0.fa"){
    gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/ppent-intergenic.tab", sep="\t", as.is=T, header=T)
  }
  if(fasta == "primaurelia_Ir4-2_assembly_v1.fasta"){
    gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/pprim-intergenic.tab", sep="\t", as.is=T, header=T)
  }
  if(fasta == "pquadecaurelia_mac_NiA_v1.0.fa"){
    gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/pquad-intergenic.tab", sep="\t", as.is=T, header=T)
  }
  if(fasta == "ptredecaurelia_209_AP38_filtered.fa"){
    gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/ptre-intergenic.tab", sep="\t", as.is=T, header=T)
  }
  if(fasta == "pnovaurelia_mac_TE_v1.0.fa"){
    gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/pnov-intergenic.tab", sep="\t", as.is=T, header=T)
  }
  if(fasta == "sexaurelia_AZ8-4_assembly_v1.fasta"){
    gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/psex-intergenic.tab", sep="\t", as.is=T, header=T)
  }
  if(fasta == "pjenningsi_mac_M_v1.0.fa"){
    gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/pjen-intergenic.tab", sep="\t", as.is=T, header=T)
  }
  if(fasta == "psonneborni_mac_ATCC_30995_v1.0.fa"){
    gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/pson-intergenic.tab", sep="\t", as.is=T, header=T)
  }
  if(fasta == "pdecaurelia_mac_223_v1.0.fa"){
    gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/pdec-intergenic.tab", sep="\t", as.is=T, header=T)
  }
  if(fasta == "pdodecaurelia_mac_274_v1.0.fa"){
    gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/pdodec-intergenic.tab", sep="\t", as.is=T, header=T)
  }
  if(fasta == "poctaurelia.fasta"){
    gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/poct-intergenic.tab", sep="\t", as.is=T, header=T)
  }
  if(fasta == "pseptaurelia_mac_38_v1.0.fa"){
    gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/psept-intergenic.tab", sep="\t", as.is=T, header=T)
  }
  if(fasta == "ptetraurelia_mac_51.fa"){
    gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/ptet-intergenic.tab", sep="\t", as.is=T, header=T)
  }
  if(fasta == "caudatum_43c3d_assembly_v1.fa"){
    gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/pcaud-intergenic.tab", sep="\t", as.is=T, header=T)
  }
  
  for(pos in regex_positions){
    row_hit = which(gff$end_position > pos)[1]                                            #which row in intergenic gff has end_pos larger than position of regex_pos? Keep only first
    
    if(pos > gff$start_position[row_hit]){                                                                        #check the start_pos to make sure the regex_pos is larger (i.e. it is intergenic)
      unique_id = paste(pos, scaf, assembly, sep="_")
      lfeatures[[unique_id]] = c(gff$X5..strand[pos], gff$X5..strand.1[pos])                                      #put in here genes in 5+ and 5- as well as distance to both
      names(lfeatures[[unique_id]]) = c("Positive_Strand_5prime", "Negative_Strand_5prime'")
    }
  }
}


#MAIN
fasta_directory = "/N/u/tlicknac/Carbonate/Paramecium_FASTA/"                                                         #location of FASTAs
FASTAs = list.files(fasta_directory, recursive=F)                                                                     #List all files in directory
FASTAs = FASTAs[grep("*.fa", FASTAs)]                                                                                 #Rm all files that arent fasta files *this will capture .fasta and .fa

regex_file = read.table("/N/dc2/scratch/tlicknac/12nt_All_Regexs_Evalues.csv", header=T, sep=",",as.is=T)             #sort and run a handful if wanted

firstloop=TRUE

for(assembly in FASTAs){                                                                                              #Iterate through fasta file names
  cat("Starting with species: ", assembly, "\n")
  a = paste(fasta_directory, assembly, sep="")                                 
  fasta = read.fasta(a, as.string=T)                                                                                  #read in fasta 1 at a time... must preserve memory
  
  for(scaf in fasta){                                                                                                 #Iterate through each scaffold
    cat("Onto scaffold: ", getName(scaf), " of species: ", assembly, "\n")
    scaf_seq = as.character(scaf[1])                                                                                  #get sequence
    
    for(regex_row in 1:nrow(regex_file)){
      regex_seq = regex_file$Regex[regex_row]
      regex_evalue = regex_file$E_Value[regex_row]
      regex_positions = unlist(gregexpr(pattern=regex_seq, text=scaf_seq, ignore.case=T))                                     #gregexpr returns list of positions with lengths
      position_info = getPositionInfo(regex_positions, scaf, assembly)                                                   #passing function a vector of positions (integers) the scaf ID and fasta name

    }
  }
}
