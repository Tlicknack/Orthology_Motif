install.packages("XML")
library("XML")

getPWM = function(xml_data, i, df_pwm){
  for(q in 1:length(xml_data$motifs[i]$motif$scores[[1]])){
    vscores = c()
    
    for(qq in 1:4){
      value = as.character(xml_data$motifs[i]$motif$scores[[1]][q]$alphabet_array[[qq]][1])  #[[]].... A, C, G, T
      vscores = append(vscores, value)
    }
    df_pwm[q,] = vscores
  }
  return(df_pwm)
}

getMismatch = function(xml_data, i){
  number_of_seqs = length(xml_data$motifs[i]$motif$contributing_sites)
  vleft = c()
  vright = c()
  vmotifs = c()
  
  for(m in 1:number_of_seqs){
    left_flank = xml_data$motifs[i]$motif$contributing_sites[m]$contributing_site$left_flank                                    #left flank
    right_flank = xml_data$motifs[i]$motif$contributing_sites[m]$contributing_site$right_flank                                  #right flank
    motif = paste(as.character(xml_data$motifs[i]$motif$contributing_sites[m]$contributing_site$site), collapse="")             #motif site
    
    vleft = append(vleft, left_flank)  #append to each vector
    vright = append(vright, right_flank)
    vmotifs = append(vmotifs, motif)  
  }
  n_right = nchar(vright[which.max(nchar(vright))])
  n_left = nchar(vright[which.max(nchar(vleft))]) 
  n_motif = nchar(vmotifs[which.max(nchar(vmotifs))]) 
  
  vleft_mismatches = c()
  
  if(length(n_left) > 0){  #some motifs lack left flank
    for(l in 1:n_left){  #iterate through each element in string.... left seqs
      vpos_n = c()
    
      for(j in 1:number_of_seqs){ #iterate through each sequence
        vpos_n = append(vpos_n, substr(vleft[j], l, l))
      }
      elements = sort(table(vpos_n), decreasing=T)               #create a count for each element in vector... a" "t" "c" "g"
      most_common = names(elements[1])                                #return character that is most abundant 
      vleft_mismatches = append(vleft_mismatches, length(which(vpos_n != most_common)))
    }
  } else{
    vleft_mismatches = NA
  }
  
  vmotifs_mismatches = c()
  
  for(m in 1:n_motif){  
    vpos_n = c()
    
    for(j in 1:number_of_seqs){
      vpos_n = append(vpos_n, substr(vmotifs[j], m, m))
    }
      elements = sort(table(vpos_n), decreasing=T)               
      most_common = names(elements[1])                                
      vmotifs_mismatches = append(vmotifs_mismatches, length(which(vpos_n != most_common)))
  }
  
  vright_mismatches = c()
  
  if(length(n_right) > 0){ #some motifs lack right flank
    for(r in 1:n_right){
      vpos_n = c()
    
      for(j in 1:number_of_seqs){
        vpos_n = append(vpos_n, substr(vright[j], r, r))
      }
      elements = sort(table(vpos_n), decreasing=T)             
      most_common = names(elements[1])                               
      vright_mismatches = append(vright_mismatches, length(which(vpos_n != most_common)))
    }
  } else{
    vright_mismatches = NA
  }
  
  lmismatches = list(vleft_mismatches, vmotifs_mismatches, vright_mismatches)
  names(lmismatches) = c("left", "motif", "right")
  return(lmismatches)
}
#--
getDistancetoStart = function(xml_data, i){
  i = as.numeric(i)
  nseqs = length(xml_data$motifs[i]$motif$contributing_sites)
  vpos = c()
  
  for(j in 1:nseqs){
    position = xml_data$motifs[i]$motif$contributing_sites[j]$contributing_site$.attrs[2]
    vpos = append(vpos, position)
  }
  mean_pos = mean(as.numeric(vpos))
  dist_to_start = round(200-mean_pos, 1)
  return(dist_to_start)
}
#--
getCommonVariant = function(xml_data, i){
  i = as.numeric(i)
  nseqs = length(xml_data$motifs[i]$motif$contributing_sites) 
  vvariants = c()
  
  for(k in 1:nseqs){
    variant = paste(as.character(xml_data$motifs[i]$motif$contributing_sites[k]$contributing_site$site), collapse="")
    vvariants = append(vvariants, variant)
  }
  most_common_variant = names(sort(table(vvariants), decreasing=T)[1])
  return(most_common_variant)
}
#--

#MAIN
directory = "/N/dc2/scratch/tlicknac/12nt_MEME_Results_TSS/"                  #location of meme outputs
meme_files = list.files(directory, recursive=F)                           #vector of file names
firstOne = T                                                              #used to create first line of data frame  

wantPWM = T #########################################if TRUE, will only get PWM

for(subdir in meme_files){
  cat("Working on file: ", subdir, "\n")
  meme = paste(directory, subdir, "/meme.xml", sep="")                  #Concatenate path to file name
  data = xmlParse(meme)                                                 #Get xml file
  xml_data = xmlToList(data)                                            #Convert xml format to an R list
  
  len = length(xml_data$motifs)
  
  for(i in 1:len){  #iterate through 3 to 5 motifs 
    
    if(wantPWM == T){                              ####PWM ONLY
      df_pwm = data.frame(matrix(ncol=4))
      colnames(df_pwm) = c("A", "C", "G", "T")
      PWM = getPWM(xml_data, i, df_pwm)
      pwm_name = paste("/N/dc2/scratch/tlicknac/12nt_PWM_All/", subdir, "|", i, ".pwm", sep="")
      write.table(PWM, file=pwm_name, sep="\t", row.names = F)
      
    } else{
    most_common_variant = getCommonVariant(xml_data, i)                 #get most common variant
    
    consensus_seq = as.character(xml_data$motifs[[i]]$.attrs[2])        #get consensus sequence for each motif
    
    regularex = xml_data$motifs[i]$motif$regular_expression             #get regex for each motif
    regularex = gsub("[\n]", "", regularex)                             #trim \n's attached for some reason
    
    nseq = as.numeric(xml_data$motifs[i]$motif$.attrs[5])               #get number of sites for each motif
    
    e_value = as.numeric(xml_data$motifs[[i]]$.attrs[9])                #get e_value; located in motifs[1:len] then the 9th attribute
    
    dist_to_start = getDistancetoStart(xml_data, i)                     #FUNCTION, get mean distance to start codon
    
    backgroundAT = round(sum(as.numeric(xml_data$model$background_frequencies$alphabet_array[4]$value$text), as.numeric(xml_data$model$background_frequencies$alphabet_array[1]$value$text)), 3)
    
    lmismatches = getMismatch(xml_data, i)
    left_mean_mismatch = mean(lmismatches$left)                        #get means for mismatches
    right_mean_mismatch = mean(lmismatches$right)
    motif_mean_mismatch = mean(lmismatches$motif)
    
    if(firstOne == T){
      final_dataframe = data.frame(matrix(data=c(most_common_variant, consensus_seq, regularex, dist_to_start, backgroundAT, left_mean_mismatch, motif_mean_mismatch, right_mean_mismatch ,nseq, e_value, subdir), ncol=11))
      names(final_dataframe) = c("CommonVariant", "Consensus", "Regex", "DistanceToStart", "BackgroundAT", "LeftMismatch", "MotifMismatch", "RightMismatch", "Nseqs", "EValue", "File")
      firstOne=F
    } else{
      new_dataframe = data.frame(matrix(data=c(most_common_variant, consensus_seq, regularex, dist_to_start, backgroundAT, left_mean_mismatch, motif_mean_mismatch, right_mean_mismatch, nseq, e_value, subdir), ncol=11))
      names(new_dataframe) = c("CommonVariant", "Consensus", "Regex", "DistanceToStart", "BackgroundAT", "LeftMismatch", "MotifMismatch", "RightMismatch", "Nseqs", "EValue", "File")
      final_dataframe = rbind(final_dataframe, new_dataframe)
    }
  write.csv(final_dataframe, file="12nt_MEME_Final_File_TSS.csv",append=T, row.names=F)
  }
  }
}
