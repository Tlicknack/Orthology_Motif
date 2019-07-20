#This script will take FIMO's output and connect features to it.
  #E.g If motif instance is located in exon, it will be tagged as such

#lgff = list()
#lgff[["pbi"]] = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/pbi-intergenic.tab", header=T, sep="\t")
#lgff[["psex"]] = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/psex-intergenic.tab", header=T, sep="\t")
#lgff[["ptet"]] = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/ptet-intergenic.tab", header=T, sep="\t")


gffGene = read.table("~/Paramecium_GFF/psex-gene.tab", as.is=T)
fpkm = read.table("/N/u/tlicknac/Carbonate/Paramecium_FPKM/psex-fpkm_prot.tab", as.is = T, header = T)
fimo_file = "/N/dc2/scratch/tlicknac/Paramecium_-200to0/Newest_All-Aurelias_FIMO_Results/Psex/WGD00008_promoter.fasta/fimo_out/fimo.tsv"
fimo = read.table(fimo_file, header=T)
fimo$intergenic = NA
fimo$geneIDleft = NA
fimo$geneFPKMleft = NA
fimo$geneDISTANCEleft = NA
fimo$geneORIENTleft = NA
fimo$geneIDright = NA
fimo$geneFPKMright = NA
fimo$geneDISTANCEright = NA
fimo$geneORIENTright = NA
species = "psex"

setwd("/N/dc2/scratch/tlicknac/Paramecium_-200to0/Newest_All-Aurelias_Motif-Distributions/")
#Split motifs
for(motif_id in unique(fimo$motif_alt_id)){
  id_rows =  fimo[which(fimo$motif_alt_id == motif_id),]
  motif_pos = round(as.numeric(summary(fimo_row$start, fimo_row$stop)[4]), digits=0)
  
  for(i in 1:nrow(id_rows)){
    fimo_row = id_rows[i,]
    gffScaf = gffGene[which(gffGene$V1 == fimo_row$sequence_name),]
    leftGene = gffScaf[which(gffScaf$V2 < motif_pos),][length(which(gffScaf$V2 < fimo_row$start)),]
    
    # If these motifs have no gene to their left, then create a fake df with NAs
    if(nrow(leftGene) == 0){   
      leftGene = as.data.frame(matrix(rep(NA, 5), ncol=5, nrow=1))
      colnames(leftGene) = c("V1", "V2", "V3", "V4", "V5")
      intergenic = NA  #motif is on end of scaffold, and is neither intergenic nor genic
    }
    
    rightGene = gffScaf[which(gffScaf$V3 > motif_pos),][1,]
    
    # If these motifs have no gene to their right, then create a fake df with NAs
    if(nrow(rightGene) == 0){   
      rightGene = as.data.frame(matrix(rep(NA, 5), ncol=5, nrow=1))
      colnames(rightGene) = c("V1", "V2", "V3", "V4", "V5")
      intergenic = NA  #motif is on end of scaffold, and is neither intergenic nor genic
    }
    
    ##  This prevents the NA's from disrupting the range of positions in the leftGene$V3:rightGene$V2 operation
    if((is.na(leftGene$V1) == F) & (is.na(rightGene$V1) == F)){
      
      if(all(c(fimo_row$start:fimo_row$stop) %in% c(leftGene$V3:rightGene$V2))){  #Motifs are intergenic if all of their postions are in between closest genes
        intergenic = "Y"
      }
      if(identical(leftGene, rightGene)){  #identical left and right means motif is genic.. this should override the above if() which returns intergenic = "Y" when leftGene == rightGene
        intergenic = "N"
      }
    } else{  #if motif is missing either a leftGene or rightGene, then its on the end
      intergenic = "end"
    }
    ##
    id_rows$intergenic[i] = intergenic
    
    if(is.na(leftGene$V1) == F){  #only give info for leftGene if there is one
      id_rows$geneIDleft[i] = leftGene$V5
      id_rows$geneFPKMleft[i] = fpkm[which(fpkm$gene_id == leftGene$V5),"FPKM"]
      
      if(leftGene$V4 == "+"){
        id_rows$geneORIENTleft[i] = "+"
        id_rows$geneDISTANCEleft[i] = leftGene$V2 - motif_pos
      }
      if(leftGene$V4 == "-"){
        id_rows$geneORIENTleft[i] = "-"
        id_rows$geneDISTANCEleft[i] = leftGene$V3 - motif_pos
      }
    } else{
      id_rows$geneORIENTleft[i] = NA
      id_rows$geneDISTANCEleft[i] = NA
    }
    if(is.na(rightGene$V1) == F){ #same for rightGene
      id_rows$geneIDright[i] = rightGene$V5
      id_rows$geneFPKMright[i] = fpkm[which(fpkm$gene_id == rightGene$V5),"FPKM"]
      
      if(rightGene$V4 == "+"){
        id_rows$geneORIENTright[i] = "+"
        id_rows$geneDISTANCEright[i] = rightGene$V2 - motif_pos
      }
      if(rightGene$V4 == "-"){
        id_rows$geneORIENTright[i] = "-"
        id_rows$geneDISTANCEright[i] = rightGene$V3 - motif_pos
      }
    } else{
      id_rows$geneORIENTright[i] = NA
      id_rows$geneDISTANCEright[i] = NA
    }
    write.table(id_rows, file=paste( gsub(pattern=".fasta", replacement="", x=as.character(strsplit(fimo_file, "/")[[1]][9])) , "_", fimo_row$motif_id, "_", species, ".txt", sep=""))
    #write.table(id_rows, file="WGD00008_promoter_CAACTGTTTGGG_pbi.txt")
  }
}


