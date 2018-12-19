#DONE

#This script will get 200bp upstream of paraorthologous genes in all 14 paramecium spp

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

extractUpstream = function(sp, geneID, lGFF, lfasta, fname){      #Pass species name, geneID, list of (modified) GFFs and list of FASTAs for that family
  
  grep_geneID = which(lGFF[[sp]]$V5 == geneID, arr.ind=TRUE)      #go into the 5th column of the GFF for the current species, return the row number which matches the current geneID
  #cat("Grep= ", grep_geneID, "\n")
  row_geneID = lGFF[[sp]][grep_geneID, ]                          #return the contents of the entire row where the above match was found
  if(row_geneID[,4] == "+"){                                      #get Gene's starting position, depending upon its strand 
    gene_start = row_geneID[,2]   
  } else{
    gene_start = row_geneID[,3]
  }
  scaf = as.character(row_geneID[,1])                             #return the name of the scaffold holding the gene 
  #cat("scaf= ", scaf, ".....", "gene start= ", gene_start, ".....", "strand= ", row_geneID[,4], "\n")
  scafSeq = as.character(lfasta[[sp]][scaf])                      #return the entire sequence of the scaffold which holds the current gene
  if(row_geneID[,4] == "+"){                                      #extract the correct upstream sequence. if +, then its just 200bp from start of gene, if -, then its 200 from end of gene
    upstream = substr(scafSeq, gene_start-200, gene_start)
  } else {
    upstream = substr(scafSeq, gene_start, gene_start+200)
    upstream = c2s(rev(comp(s2c(upstream))))                      #must reverse complement the gene on the minus strand. Seqinr has a strange way of doing this
  }
  return(upstream)
  #cat("Upstream sequence for ", geneID, "has been extracted successfully. \n")
}


getPOFFid = function(vsp, POFF, lGFF, lfasta) {                           #Pass vector of species names, POFF of family of interest, list of GFFs and FASTAs for that family
  
  for(i in 1:nrow(POFF)) {                                                #use i to iterate through each row number
    fname = paste(names(POFF[4]), "_", i, ".fasta", sep="")               #return file name based on the species name in column 4 
    cat("Parsing line ", i, "\n")           #fun printed message
    vupstream = c()                                                       #create empty vector to store all upstream sequences for genes in this row
    
    for(sp in vsp) {                                                      #iterate through species names
      sp_family = as.character(POFF[i,sp])                                #get proteinIDs from each row (based on iterator) and column (based on species name which is the header)
      vsplit = strsplit(sp_family, ",")[[1]]                              #some columns have two geneIDs, must split them on the comma that separates them 
      
      for(geneID in vsplit){                                              #iterate through 1 or 2 geneIDs (if there's only 1, then it's fine too)

        
        if( (geneID != "NA")&(geneID != ".")&(geneID != "*")){                                                #only proceed if the geneID != '*' ... *'s are used when no ortho-paralog were found. It returns no upstream sequence
          cat("Now extracting sequence upstream of ", geneID, " (which is in", sp, ")\n")     #might as well remove*'s before running proceeding function                
          upstreamSeq = extractUpstream(sp, geneID, lGFF, lfasta, fname)
          cat("Sequence extracted: '", upstreamSeq, "'\n")
          vupstream[geneID] = upstreamSeq                                 #add the upstream sequence into a vector with the position given by the geneID. can easily call that position w/ geneID
        }
      }
    }
    if(length(vupstream) >= 4) {                                          #only create a fasta of these sequences if there are at least 4. MEME will not be effective with 3 or less
      writeVectorAsFasta(vupstream, fname)
      cat("line", i, "done.\n-----------------------\n")                  #print that line row is complete ... I'll keep this print-out
    }
  }
}


writeVectorAsFasta <- function(vseq, fname){                          #Pass vector of sequences from extractUpstream and name of file. This will write to fasta file as output
  #cat("FNAME: ", fname, "\n")
  cat("", file=fname)                                                 #reset the file to empty if it already exists
  #print(vseq)
  for(id in names(vseq)){                                           #iterate through geneIDs which are names for each corresponding sequence
    cat("Writing ", id ," to file '", fname, "'")
    seq_len = nchar(vseq[[id]])                                     #return length of sequence for each geneID
    
    if(seq_len > 10){                                               #only write to fasta if the length is greater than 10 (just in case something slipped by)
      cat(">", id, "\n", vseq[[id]], "\n", sep="", file=fname, append=T)       
      #write to fasta file. >geneID \n "sequence of that geneID". name the file after the species family and row number, append until this row is finished (fname will change next row)
    }
  }
}

#MAIN                                                                                 #Run in home directory, output into scratch
poff = read.table("Paramecium_POFF/all_aurelias.poff", header=T, sep="\t")
vsp = colnames(poff)[4:(ncol(poff)-1)]
lfasta = list()
lfasta[["pbi"]] = read.fasta("Paramecium_FASTA/biaurelia_V1-4_assembly_v1.fa", as.string=TRUE)          
lfasta[["ppent"]] = read.fasta("Paramecium_FASTA/ppentaurelia_mac_87_v1.0.fa", as.string=TRUE)
lfasta[["pprim"]] = read.fasta("Paramecium_FASTA/primaurelia_Ir4-2_assembly_v1.fasta", as.string=TRUE)
lfasta[["pquad"]] = read.fasta("Paramecium_FASTA/pquadecaurelia_mac_NiA_v1.0.fa", as.string=TRUE)
lfasta[["ptre"]] = read.fasta("Paramecium_FASTA/ptredecaurelia_209_AP38_filtered.fa", as.string=TRUE)
lfasta[["pnov"]] = read.fasta("Paramecium_FASTA/pnovaurelia_mac_TE_v1.0.fa", as.string=TRUE)
lfasta[["psex"]] = read.fasta("Paramecium_FASTA/sexaurelia_AZ8-4_assembly_v1.fasta", as.string=TRUE)
lfasta[["pjen"]] = read.fasta("Paramecium_FASTA/pjenningsi_mac_M_v1.0.fa", as.string=TRUE)
lfasta[["pson"]] = read.fasta("Paramecium_FASTA/psonneborni_mac_ATCC_30995_v1.0.fa", as.string=TRUE)
lfasta[["pdec"]] = read.fasta("Paramecium_FASTA/pdecaurelia_mac_223_v1.0.fa", as.string=TRUE)
lfasta[["pdodec"]] = read.fasta("Paramecium_FASTA/pdodecaurelia_mac_274_v1.0.fa", as.string=TRUE)
lfasta[["poct"]] = read.fasta("Paramecium_FASTA/poctaurelia.fasta", as.string=TRUE)
lfasta[["psept"]] =read.fasta("Paramecium_FASTA/pseptaurelia_mac_38_v1.0.fa", as.string=TRUE)
lfasta[["ptet"]] = read.fasta("Paramecium_FASTA/ptetraurelia_mac_51.fa", as.string=TRUE)

lgff = list()
for(spp in vsp){
  gff_name = paste("Paramecium_GFF/", spp, "-gene.tab", sep="")
  gff_read = read.table(gff_name, sep="\t", as.is=T)
  lgff[[spp]] = gff_read
}

getPOFFid(vsp, poff, lgff, lfasta)
