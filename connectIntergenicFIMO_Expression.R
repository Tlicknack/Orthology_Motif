# This script will do everything to determine if a motif has evidence of functionality or not
# We will work with FIMO outputs from intergenic-fasta's


library(seqinr)
fpkm_ptet = read.table("~/Paramecium_FPKM/ptet-fpkm_prot.tab", header=T)
fastaInter_ptet = read.fasta("~/Paramecium_FASTA/Intergenic_FASTA/ptet-intergenic.fasta", forceDNAtolower = T, as.string = T)
fimo_file = "/N/dc2/scratch/tlicknac/Paramecium_-200to0/Newest_All-Aurelias_FIMO_Results/Ptet-Intergenic/"

setwd("/N/dc2/scratch/tlicknac/Paramecium_-200to0/Newest_All-Aurelias_FIMO_Results/Ptet-Intergenic/")

for(wgd in list.files("./")[150:155]){
  tmp = paste(getwd(), "/", wgd, "/fimo_out/fimo.tsv", sep="")
  fimo_spp = read.table(tmp, header=T, row.names = NULL)
  
  for(motif in unique(fimo_spp$motif_alt_id)){
    motif_rows = fimo_spp[which(fimo_spp$motif_alt_id == as.character(motif)),]
    motif_rows$FivePlusFPKM = NA
    motif_rows$FivePlusDistance  = NA
    motif_rows$FiveMinusFPKM  = NA
    motif_rows$FiveMinusDistance  = NA
    
    for(i in 1:nrow(motif_rows)){
      current_row = motif_rows[i,]
      scaf_len = nchar(as.character(fastaInter_ptet[which(getName(fastaInter_ptet) == current_row$sequence_name)]))
      fivePlus = strsplit(as.character(current_row$sequence_name), "_")[[1]][3]
      fiveMinus = strsplit(as.character(current_row$sequence_name), "_")[[1]][5]
      
      if(fivePlus != "NA"){
        motif_rows$FivePlusFPKM[i] = fpkm_ptet[which(fpkm_ptet$gene_id == fivePlus),10]  
        motif_rows$FivePlusDistance[i] = scaf_len - current_row$start
      } 
      
      if(fivePlus == "NA"){
        motif_rows$FivePlusFPKM[i] = NA
        motif_rows$FivePlusDistance[i] = NA
      }
      
      if(fiveMinus != "NA"){
        motif_rows$FiveMinusFPKM[i] = fpkm_ptet[which(fpkm_ptet$gene_id == fiveMinus),10]  
        motif_rows$FiveMinusDistance[i] = current_row$stop
      } 
      
      if(fiveMinus == "NA"){
        motif_rows$FiveMinusFPKM[i] = NA
        motif_rows$FiveMinusDistance[i] = NA
      }
    }
    setwd("/N/dc2/scratch/tlicknac/Paramecium_-200to0/Newest_All-Aurelias_Motif-Distributions/Ptet-Intergenic/")
    write.table(motif_rows, col.names = T, file = paste(gsub(".fasta","", wgd), "-", current_row$motif_id, ".txt" , sep=""))
    setwd("/N/dc2/scratch/tlicknac/Paramecium_-200to0/Newest_All-Aurelias_FIMO_Results/Ptet-Intergenic/")
  }
}
#This will return a directory of files that look like this:
# "motif_id" "motif_alt_id" "sequence_name" "start" "stop" "strand" "score" "p.value" "q.value" "matched_sequence" "FivePlusFPKM" "FivePlusDistance" "FiveMinusFPKM" "FiveMinusDistance"
# "4524" "TGTTWTAA" "MEME-2" "scaffold51_18|5'+_PTET.51.1.P0180253_5'-_PTET.51.1.P0180252_3'+_NA_3'-_NA" 1 8 "+" 14.1429 2.65e-05 0.107 "tgttttaa" 839.881 620 345.462 8

#Statistical analysis of motifs
ptet_fpkm_summary = summary(fpkm_ptet$FPKM)

for(motif_file in list.files("./", strsplit(strsplit(fimo_file, "/")[[1]][6], ".fasta")[[1]][1])){
  fimo_expr = read.table(paste(working_dir, motif_file, sep=""))
  median_expr = summary(c(fimo_expr$FivePlusFPKM, fimo_expr$FiveMinusFPKM))[3]
  mean_expr_sig = as.numeric(t.test(c(fimo_expr$FivePlusFPKM, fimo_expr$FiveMinusFPKM), fpkm_ptet$FPKM)[3])
  dist_expr_cor = as.numeric(cor.test(c(fimo_expr$FivePlusFPKM, fimo_expr$FiveMinusFPKM), c(fimo_expr$FivePlusDistance, fimo_expr$FiveMinusDistance))[3])
  
  out_tab = data.frame(matrix(c(motif_file, "ptet", median_expr, mean_expr_sig, dist_expr_cor),nrow = 1, ncol=5))
  colnames(out_tab) = c("MotifFile", "Species", "MedianExpr", "PvalueMeanExpr", "PvalueCorExprDist")
  write.table(out_tab, file="Motif-Expression-Significance.txt")
}











