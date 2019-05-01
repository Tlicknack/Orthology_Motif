#This script will create a fasta file of Paramecium Upstream regions with a Sharp Promoter
  #Additionally, a background file with all other upstream regions from genes with broader promoters
#NOTE: This uses only a single TSR table as an input. I will later modify this to include all TSR data.
  #For now, if a single gene from Pdec has a sharp promoter, I'm assuming it's paralogs and orthologs do as well

#INPUT:
  #POFF table for orthology relationships
  #TSR table with promoter shapes
  #FASTA file of genome assemblies
  #GFF file of genome annotations
#OUTPUT:
  #Two FASTA files for upstream sequences from genes with broad and sharp promoters

library("seqinr")

getUpstream = function(poff_genes, lgff, lfasta){
  vupstreams = c()
  vPoffGenes = c()
  
  if(length(poff_genes) > 0){#####   make sure we matched a POFF row to a TSR gene feature
    poff_genes = unlist(strsplit(poff_genes, ","))                                              #split the POFF row on a comma
    wgd = poff_genes[1]                                                                         #get wgd code 
    poff_genes = poff_genes[-1]                                                                 #remove wgd code
    
    if(is.na(any(poff_genes == "NA")) == F){
      poff_genes = poff_genes[-which(poff_genes == "NA")]  #cleanup  
    }
    if(is.na(any(poff_genes == ".")) == F){
      poff_genes = poff_genes[-which(poff_genes == ".")]   #cleanup 
    }
    if(any(is.na(poff_genes))){
      poff_genes = poff_genes[-which(is.na(poff_genes))]   #cleanup 
    }

      for(poff_gene in poff_genes){ #
      actual_gff_row = 0  #default. If none of the POFF row's genes were found, this is unchanged
      actual_species = 0  #default. We need to save species name to get correct FASTA
      ##
      for(species in names(lgff)){                              ## Find the actual GFF file for that species, and return the row that matches the poff gene
        gff = as.data.frame(lgff[species])                                  #get the annotation for the current species
        gff_row = gff[which(regexpr(poff_gene, gff[1:nrow(gff),5]) > 0),]   #get the row in the GFF that matches the current POFF gene
        
        if(is.na(gff_row[1,1]) == F){ 
          actual_gff_row = gff_row  #Matching gff row becomes actual
          actual_species = species  #matching species becomes actual
          cat(poff_gene, " found for species ", species, " \n", sep="")
        }
        cat("Species: ", species, " was checked. \n")
      }
      ##
      
      ###
      if(actual_species != 0){
        sp_fasta = lfasta[[actual_species]]
        scaf_seq = as.character(sp_fasta[as.character(actual_gff_row[1,1])])
        
        if(as.character(actual_gff_row[1,4]) == "+"){
          upstream = substr(scaf_seq,  as.numeric(actual_gff_row[1,2]) - 200,  as.numeric(actual_gff_row[1,2]))
        }
        if(as.character(actual_gff_row[1,4]) == "-"){
          upstream = c2s(rev(comp(s2c(substr(scaf_seq,  as.numeric(actual_gff_row[1,3]),  as.numeric(actual_gff_row[1,3]) + 200))))) 
        }
      }
      vPoffGenes = append(vPoffGenes, poff_gene)
      vupstreams = append(vupstreams, upstream)
    } #
  }#####
  return(list(vupstreams, vPoffGenes))
}

poff = read.table("/N/u/tlicknac/Carbonate/Paramecium_POFF/all_aurelias-cut-mod.poff", header=T)  #use poff table with modified pcaud and pmultimic names
poff = poff[,-2]     

vsp = colnames(poff)
vsp = vsp[-1]

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

tsr = read.table("/N/dc2/scratch/tlicknac/Data-TSR/Pdec_TSR_Merged.txt", header=T, as.is=T)
genic_tsr = tsr[which(is.na(tsr$featureID)==F),]
genic_tsr = genic_tsr[-which(genic_tsr$tsrWdth > 100)]  #reads got clustered together, this isnt a real signal from a single gene
rm(tsr)

#genic_tsr = genic_tsr[which(genic_tsr$tsrMSI != 1),]  #if R is pooping out on cluster
#genic_tsr = genic_tsr[which(genic_tsr$tsrMSI == 1),]

#MAIN
lsharp_seqs = list()
lsharp_annots = list()
lbroad_seqs = list()
lbroad_annots = list()

#Loop through TSR data
for(tsr_row in 1:nrow(genic_tsr)){
  tsr_gene = gsub("223.1.G", "223.1.P", genic_tsr$featureID[tsr_row])                         #change gene to proteinID
  tsr_shape = as.numeric(genic_tsr$tsrMSI[tsr_row])                                           #get the shape of that TSR
  poff_genes = as.character(unlist(unname(poff[which(regexpr(tsr_gene, poff$pdec) > 0),])))   #find row in POFF table that matches TSR gene feature
  
  if(tsr_shape == 1){
    lsharp = getUpstream(poff_genes, lgff, lfasta)
    lsharp_seqs = append(lsharp_seqs, lsharp[[1]])
    lsharp_annots = append(lsharp_annots, lsharp[[2]])
  } 
  if(tsr_shape < 1){
    lbroad = getUpstream(poff_genes, lgff, lfasta)
    lbroad_seqs = append(lbroad_seqs, lbroad[[1]])
    lbroad_annots = append(lbroad_annots, lbroad[[2]])
  }
  cat("Row ", tsr_row, " done. \n")
}
write.fasta(lsharp_seqs, unlist(lsharp_annots), file.out="Sharp-TSR_Upstream.fasta")
write.fasta(lbroad_seqs, unlist(lbroad_annots), file.out="Broad-TSR_Upstream.fasta")
