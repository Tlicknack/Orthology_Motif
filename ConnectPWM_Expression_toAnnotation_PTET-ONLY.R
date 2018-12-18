#DONE

#This script will combine annotation and expression values in vegetative growth for high-score motifs

##NOTE: Written only for ptet. Make


library(seqinr)

ptet_fa = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/ptetraurelia_mac_51.fa", as.string=T)
ptet_fpkm = read.table("/N/u/tlicknac/Carbonate/Paramecium_FPKM/ptet-fpkm_prot.tab", header=T, sep="\t")

annotated_motif = read.table("pwm-distribution-histone-GCCAATCACAGT_With-Annotation.csv", header =T, as.is = T, sep = ",")
intergenic_motif = annotated_motif[which(annotated_motif$Intergenic == "Y"),]

df = data.frame(matrix(ncol=2, nrow = nrow(intergenic_motif)))
colnames(df) = c("FPKM_FivePrimePlus", "FPKMFivePrimeMinus")
intergenic_motif = cbind(intergenic_motif, df)

for(i in 1:nrow(intergenic_motif)){
  row = intergenic_motif[i,]
    
  if(row$Strand == "+"){
    gene = row$FivePrimeGenePlusStrand
    fpkm_row = ptet_fpkm[which(ptet_fpkm$gene_id == gene),]
    fpkm = as.numeric(fpkm_row$FPKM)
    if(length(fpkm) > 0){
      intergenic_motif$FPKM_FivePrimePlus[i] = fpkm
    }
  }
    
  if(row$Strand == "-"){
    gene = row$FivePrimeGeneMinusStrand
    fpkm_row = ptet_fpkm[which(ptet_fpkm$gene_id == gene),]
    fpkm = as.numeric(fpkm_row$FPKM)
      
    if(length(fpkm) > 0){
      intergenic_motif$FPKMFivePrimeMinus[i] = fpkm
    }
  }
  scaffold = row$Scaffold
  scaf_seq = as.character(ptet_fa[which(getName(ptet_fa) == scaffold)])
  
  if(row$Strand == "+"){
    motif_seq = substr(scaf_seq, row$Position, nchar(row$Consensus_Motif)+row$Position)
    intergenic_motif$Consensus_Motif[i] = toupper(motif_seq)
  } 
  if(row$Strand == "-"){
    motif_seq = substr(scaf_seq, row$Position-nchar(row$Consensus_Motif), row$Position)
    intergenic_motif$Consensus_Motif[i] = toupper(c2s (rev (comp (s2c (motif_seq) ) ) ) )
  }
  cat("Writing to line: ", i, "\n", sep="")
}
