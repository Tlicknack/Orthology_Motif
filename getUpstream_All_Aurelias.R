#THis program will get 200bp upstream of paraorthologous genes in all 14 paramecium spp
#INPUT:
  #fasta file for assemblies (all 14 spp)
  #modified gff that contains only genic regions
  #poff table of paraorthologs
#OUTPUT:
  #fasta files of upstream sequences

#source("http://bioconductor.org/biocLite.R")
#biocLite("seqinr")
#install.packages("seqinr")
library("seqinr")

extractUpstream = function(sp, geneID, lgff, lfasta){            #Pass species name, geneID, list of (modified) GFFs and list of FASTAs for that family
  grep_geneID = which(lgff[[sp]]$V5 == geneID, arr.ind=TRUE)      #go into the 5th column of the GFF for the current species, return the row number which matches the current geneID
  
  if(length(grep_geneID) > 0){
    row_geneID = lgff[[sp]][grep_geneID, ]                          #return the contents of the entire row where the above match was found
    scaf = as.character(row_geneID$V1)                              #return the name of the scaffold holding the gene 
    scafSeq = as.character(lfasta[[sp]][scaf])                      #return the entire sequence of the scaffold which holds the current gene
  
    if(row_geneID$V4 == "+"){                                       #get Gene's starting position, depending upon its strand 
      gene_start = row_geneID$V2  
      upstream = substr(scafSeq, gene_start-200, gene_start)  
    } 
    if(row_geneID$V4 == "-"){
      gene_start = row_geneID$V3
      upstream = substr(scafSeq, gene_start, gene_start+200)
      upstream = c2s(rev(comp(s2c(upstream))))                      #must reverse complement the gene on the minus strand. Seqinr has a strange way of doing this
    }
    return(upstream)
  }
}

#MAIN                                                                                 #Run in home directory, output into scratch
#prepare poff, fastas, gffs
poff = read.table("/N/u/tlicknac/Carbonate/Paramecium_POFF/all_aurelias-cut-mod.poff", header=T)  #use poff table with modified pcaud and pmultimic names
#poff = read.table("Paramecium_POFF/all_aurelias.flipped", header=T)
poff = poff[,-2]                     #remove number of genes

vsp = colnames(poff)
vsp = vsp[-1]             #remove WGD from names

lfasta = list()
lfasta[["pbi"]] = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/biaurelia_V1-4_assembly_v1.fa", as.string=TRUE)          
lfasta[["ppent"]] = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/ppentaurelia_mac_87_v1.0.fa", as.string=TRUE)
lfasta[["pprim"]] = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/primaurelia_Ir4-2_assembly_v1.fasta", as.string=TRUE)
lfasta[["pquad"]] = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/pquadecaurelia_mac_NiA_v1.0.fa", as.string=TRUE)
lfasta[["ptre"]] = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/ptredecaurelia_209_AP38_filtered.fa", as.string=TRUE)
lfasta[["pnov"]] = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/pnovaurelia_mac_TE_v1.0.fa", as.string=TRUE)
lfasta[["psex"]] = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/sexaurelia_AZ8-4_assembly_v1.fasta", as.string=TRUE)
lfasta[["pjen"]] = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/pjenningsi_mac_M_v1.0.fa", as.string=TRUE)
lfasta[["pson"]] = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/psonneborni_mac_ATCC_30995_v1.0.fa", as.string=TRUE)
lfasta[["pdec"]] = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/pdecaurelia_mac_223_v1.0.fa", as.string=TRUE)
lfasta[["pdodec"]] = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/pdodecaurelia_mac_274_v1.0.fa", as.string=TRUE)
lfasta[["poct"]] = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/poctaurelia.fasta", as.string=TRUE)
lfasta[["psept"]] =read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/pseptaurelia_mac_38_v1.0.fa", as.string=TRUE)
lfasta[["ptet"]] = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/ptetraurelia_mac_51.fa", as.string=TRUE)
lfasta[["pmultimic"]] = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/multimicronucleatum_MO3c4_assembly_v1.fa", as.string=TRUE)
lfasta[["pcaud"]] = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/caudatum_43c3d_assembly_v1.fa", as.string=TRUE)

lgff = list()
for(spp in vsp){
  gff_name = paste("/N/u/tlicknac/Carbonate/Paramecium_GFF/", spp, "-gene.tab", sep="")
  gff_read = read.table(gff_name, sep="\t")
  lgff[[spp]] = gff_read
}
lgff[["pcaud"]] = caud_gfff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/pcaud-oldgene.tab", sep="\t")

#loop through poff
for(i in 1:nrow(poff)) {                                              
  fname = paste(poff$WGD[i], "_promoter", ".fasta", sep="")                     #unique name for each WGD family
  cat("Parsing line ", i, "\n")                         #fun printed message
  vupstream = c()                                                       #create empty vector to store all upstream sequences for genes in this row
#loop through species names  
  for(sp in vsp) {                                          
    sp_family = as.character(poff[i,sp])                                #get proteinIDs from each row (based on iterator) and column (based on species name which is the header)
    vsplit = strsplit(sp_family, ",")[[1]]                              #some columns have two geneIDs, must split them on the comma that separates them 
#loop through genes in pair of genes    
    for(geneID in vsplit){                                              #iterate through 1 or 2 geneIDs (if there's only 1, then it's fine too)
      
      if( (is.na(geneID) == F)&(geneID != ".")&(geneID != "*") ){                                #only proceed if the geneID != '*' ... *'s are used when no ortho-paralog were found. It returns no upstream sequence
        cat("Now extracting sequence upstream of ", geneID, " (which is in", sp, ")\n")     #might as well remove*'s before running proceeding function                
        upstreamSeq = extractUpstream(sp, geneID, lgff, lfasta)
        cat("Sequence extracted: '", upstreamSeq, "'\n")
        
        if(is.null(upstreamSeq) == F){
          vupstream[geneID] = upstreamSeq                                 #add the upstream sequence into a vector with the position given by the geneID. can easily call that position w/ geneID
        }
      }
    }
  }
  setwd("/N/dc2/scratch/tlicknac/Newest_All-Aurelias_Upstream")
  
  if(length(names(vupstream)) >4){
    write.fasta(as.list(vupstream), names(vupstream), file.out = fname)
    cat("line", i, "done.\n-----------------------\n")                  #print that row i is complete 
  } else{
    cat("The row has too few genes!\n---------------------\n")
  }

}
