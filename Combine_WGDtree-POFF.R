#Combine POFF table of Paralogs and orthologs (limited to 2 genes per species) with WGD.tree files containing up to 8 paralogs per species from WGD3

poff_table = read.table("/N/u/tlicknac/Carbonate/Paramecium_POFF/all_aurelias-cut-mod.poff", header=T, row.names = NULL, as.is = T)

tree_dir = "/N/u/tlicknac/Carbonate/Paramecium_WGD3_Trees/"
tree_files = list.files(tree_dir, pattern = "*.")

#Create list of data frames containing each paralog tree
trees = list()

for(species in tree_files){
  input = paste(tree_dir, species, sep="")
  tree_name = substr(species, start = 0, stop = 4)
  trees[[tree_name]] = read.table(input, header=T, as.is = T)
}

newdf = data.frame()

#Combine paralog trees based on WGD1 POFF table
for(i in 1:nrow(poff_table)){
  poff_row = poff_table[i,]
  
  pbi = as.character(poff_row$pbi)
  bi1 = strsplit(pbi, ",")[[1]][1]
  bi2= strsplit(pbi, ",")[[1]][2]
  
  if(bi1 != "." & is.na(bi1) == F){
    
  }
  if(bi2 != "." & is.na(bi2) == F){
    
  }
    
  pdec = as.character(poff_row$pdec)
  dec1 = strsplit(pdec, ",")[[1]][1]
  dec2 = strsplit(pdec, ",")[[1]][2]
  
  if(dec1 != "." & is.na(dec1) == F){
    dec1row = trees[["pdec"]][with( trees[["pdec"]] , trees[["pdec"]][,2] == dec1 | trees[["pdec"]][,3] == dec1 | trees[["pdec"]][,4] == dec1 | trees[["pdec"]][,5] == dec1 | trees[["pdec"]][,6] == dec1 | trees[["pdec"]][,6] == dec1 | trees[["pdec"]][,7] == dec1 | trees[["pdec"]][,8] == dec1 | trees[["pdec"]][,9] == dec1  ) ,]
    dec1hits = as.character(dec1row[grep("PDEC.223*", dec1row)])
    }
  if(dec2 != "." & is.na(dec2) == F){
    dec2row = trees[["pdec"]][with( trees[["pdec"]] , trees[["pdec"]][,2] == dec2 | trees[["pdec"]][,3] == dec2 | trees[["pdec"]][,4] == dec2 | trees[["pdec"]][,5] == dec2 | trees[["pdec"]][,6] == dec2 | trees[["pdec"]][,6] == dec2 | trees[["pdec"]][,7] == dec2 | trees[["pdec"]][,8] == dec2 | trees[["pdec"]][,9] == dec2  ) ,]
    dec2hits = as.character(dec2row[grep("PDEC.223*", dec2row)])
  }
  
  pdodec = as.character(poff_row$pdodec)
  dodec1 = strsplit(pdodec, ",")[[1]][1]
  dodec2 = strsplit(pdodec, ",")[[1]][2]
  
  pjen = as.character(poff_row$pjen)
  jen1 = strsplit(pjen, ",")[[1]][1]
  jen2 = strsplit(pjen, ",")[[1]][2]
  
  pnov = as.character(poff_row$pnov)
  nov1 = strsplit(pnov, ",")[[1]][1]
  nov2 = strsplit(pnov, ",")[[1]][2]
  
  if(nov1 != "." & is.na(nov1) == F){
    nov1row = trees[["pnov"]][with( trees[["pnov"]] , trees[["pnov"]][,2] == nov1 | trees[["pnov"]][,3] == nov1 | trees[["pnov"]][,4] == nov1 | trees[["pnov"]][,5] == nov1 | trees[["pnov"]][,6] == nov1 | trees[["pnov"]][,6] == nov1 | trees[["pnov"]][,7] == nov1 | trees[["pnov"]][,8] == nov1 | trees[["pnov"]][,9] == nov1  ) ,]
    nov1hits = as.character(nov1row[grep("PNOV.TE*", nov1row)])
  }
  if(dec2 != "." & is.na(dec2) == F){
    nov2row = trees[["pnov"]][with( trees[["pnov"]] , trees[["pnov"]][,2] == nov2 | trees[["pnov"]][,3] == nov2 | trees[["pnov"]][,4] == nov2 | trees[["pnov"]][,5] == nov2 | trees[["pnov"]][,6] == nov2 | trees[["pnov"]][,6] == nov2 | trees[["pnov"]][,7] == nov2 | trees[["pnov"]][,8] == nov2 | trees[["pnov"]][,9] == nov2  ) ,]
    nov2hits = as.character(nov2row[grep("PNOV.TE*", nov2row)])
  }
  novHits = unique(c(nov1hits, nov2hits))
  
  poct = as.character(poff_row$poct)
  oct1 = strsplit(poct, ",")[[1]][1]
  oct2 = strsplit(poct, ",")[[1]][2]
  
  ppent = as.character(poff_row$ppent)
  pent1 = strsplit(ppent, ",")[[1]][1]
  pent2 = strsplit(ppent, ",")[[1]][2]
  
  pquad = as.character(poff_row$pquad)
  quad1 = strsplit(pquad, ",")[[1]][1]
  quad2 = strsplit(pquad, ",")[[1]][2]
  
  pprim = as.character(poff_row$pprim)
  prim1 = strsplit(pprim, ",")[[1]][1]
  prim2 = strsplit(pprim, ",")[[1]][2]
  
  psep = as.character(poff_row$psept)
  sep1 = strsplit(psep, ",")[[1]][1]
  sep2 = strsplit(psep, ",")[[1]][2]
  
  psex = as.character(poff_row$psex)
  sex1 = strsplit(psex, ",")[[1]][1]
  sex2 = strsplit(psex, ",")[[1]][2]
  
  pson = as.character(poff_row$pson)
  son1 = strsplit(pson, ",")[[1]][1]
  son2 = strsplit(pson, ",")[[1]][2]
  
  ptet = as.character(poff_row$ptet)
  tet1 = strsplit(ptet, ",")[[1]][1]
  tet2 = strsplit(ptet, ",")[[1]][2]
  
  ptre = as.character(poff_row$ptre)
  tre1 = strsplit(ptre, ",")[[1]][1]
  tre2 = strsplit(ptre, ",")[[1]][2]
}



