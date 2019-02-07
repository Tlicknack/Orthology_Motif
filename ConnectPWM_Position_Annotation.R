#DONE
#This will take in a table containing positions with high scores from PWM across genomes of interest
#INPUT:
  #DF with "Consensus_Motif", "POFF_Row", "Species", "Scaffold", "Position", "PWM_Score"
  #Intergenic GFFs with: "scaffold"      "start_position"        "end_position"  "5'+strand"     "5'-strand"     "3'+strand"     "3'-strand"

#OUTPUT
  #DF with the above + genic features, position to start

#MAIN

#Load intergenic gffs
lgff = list()
lgff[["pbi"]] = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/pbi-intergenic.tab", header=T, sep="\t")
lgff[["psex"]] = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/psex-intergenic.tab", header=T, sep="\t")
lgff[["ptet"]] = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/ptet-intergenic.tab", header=T, sep="\t")


#Load pwm scoring data frame
input_file = "pwm-distribution_WGD13807|1.csv"                                  ##CHange this for each new PWM file
file_directory_path = "/N/dc2/scratch/tlicknac/Data-PWM_Distribution/"
input = paste(file_directory_path, input_file, sep="")
pwm_table = read.csv(input, header=T, as.is=T)
new_column = data.frame(matrix(ncol=6, nrow=nrow(pwm_table)))
colnames(new_column) = c("Intergenic", "FivePrimeGenePlusStrand", "FivePrimeGeneMinusStrand", "DistancetoFivePlus", "DistanceToFiveMinus", "Intragenic")
pwm_table = cbind(pwm_table, new_column)              


for(i in 1:nrow(pwm_table)){ 
  pwm_row = pwm_table[i,]
  
  if(pwm_row$Species == "ptetraurelia_mac_51.fa"){
    gff_tab = lgff[["ptet"]]
  } 
  if(pwm_row$Species == "biaurelia_V1-4_assembly_v1.fa"){
    gff_tab = lgff[["pbi"]]
  } 
  if(pwm_row$Species == "sexaurelia_AZ8-4_assembly_v1.fasta"){
    gff_tab = lgff[["psex"]]
  }
  
  if(pwm_row$Scaffold %in% gff_tab$scaffold){                                   #only do attempt this if the scaffold from the PWM table has annotation in the GFF
    gff_current_scaf = gff_tab[which(gff_tab$scaffold ==  pwm_row$Scaffold),]   #all intergenic regions in GFF
    hits = which(pwm_row$Position > gff_current_scaf$start_position)            
    firstrow = gff_current_scaf[length(hits),][1,]                              #take the last row where the position is larger than the start position...
    
    #Some hits won't be within intergenic gff coordinates... theyre on either ends of scaffold
    if(is.na(firstrow$scaffold)){ 
       
      if(pwm_row$Position < gff_current_scaf$start_position[1]){
        intergenic = "Y"
        five_plus = as.character(gff_current_scaf$X5..strand[1])
        five_minus = as.character(gff_current_scaf$X5..strand.1[1])
        cat("Hit at the beginning of scaffold: ", gff_current_scaf$scaffold[1], "\n")
        scaf_begin = T
        scaf_end = F
        intragenic = NA
      }
      if(pwm_row$Position > gff_current_scaf$end_position[nrow(gff_current_scaf)]){
        intergenic = "Y"
        five_plus = as.character(gff_current_scaf$X5..strand[nrow(gff_current_scaf)])
        five_minus = as.character(gff_current_scaf$X5..strand.1[nrow(gff_current_scaf)])
        cat("Hit at the end of scaffold: ", gff_current_scaf$scaffold, "\n")
        scaf_begin = F
        scaf_end = T
        intragenic = NA
      }
    } else {
      
      #For those not on scaffold ends, firstrow should exist... 
        if(pwm_row$Position < firstrow$end_position){
          intergenic = "Y"
          five_plus =  as.character(firstrow[4]$X5..strand)
          five_minus = as.character(firstrow[5]$X5..strand)
          cat("HIT AT pwm_row POSITION ", pwm_row$Position, " BETWEEN START POSITION ", firstrow$start_position, " AND END POSITION ", firstrow$end_position, "\n")
          scaf_begin = F
          scaf_end = F
          intragenic = NA
        }
        if(pwm_row$Position > firstrow$end_position){
          intergenic = "N"
          five_plus = NA
          five_minus = NA
          cat("No hit at pwm_row position ", pwm_row$Position, " between start_position ", firstrow$start_position, " and end_position ", firstrow$end_position, "\n")
          scaf_begin = F
          scaf_end = F
          
          if(pwm_row$Position > firstrow$end_position){
            intragenic = as.character(firstrow$X5..strand)
            if(is.na(intragenic)){
              intragenic = as.character(firstrow$X3..strand.1)
            }
          }
          if(pwm_row$Position < firstrow$start_position){
            intragenic = as.character(firstrow$X5..strand.1)
            if(is.na(intragenic)){
              intragenic = as.character(firstrow$X3..strand)
            }
          }
        }
      }
    #get distance to start codon based on orientation
    if(intergenic == "Y"){
      dist_to_plus = abs( as.numeric(firstrow$end_position) - as.numeric(pwm_row$Position) ) 
      dist_to_minus = abs ( as.numeric(firstrow$start_position) - as.numeric(pwm_row$Position) )
      
      if(scaf_begin == T){
        dist_to_plus = abs( as.numeric(gff_current_scaf$end_position[1]) - as.numeric(pwm_row$Position[1]) ) 
        dist_to_minus = abs ( as.numeric(gff_current_scaf$start_position[1]) - as.numeric(pwm_row$Position[1]) )
      }
      
      if(scaf_end == T){
        dist_to_plus = abs( as.numeric(gff_current_scaf$end_position[nrow(gff_current_scaf)]) - as.numeric(pwm_row$Position) ) 
        dist_to_minus = abs ( as.numeric(gff_current_scaf$start_position[nrow(gff_current_scaf)]) - as.numeric(pwm_row$Position) )
      }
    }
    
    if(intergenic == "N"){
      dist_to_plus = NA
      dist_to_minus = NA
    }
    
    cat("Writing to line: ", i, "\n")
    pwm_table$Intergenic[i] = intergenic
    pwm_table$FivePrimeGenePlusStrand[i] = five_plus
    pwm_table$FivePrimeGeneMinusStrand[i] = five_minus
    pwm_table$DistancetoFivePlus[i] = dist_to_plus
    pwm_table$DistanceToFiveMinus[i] = dist_to_minus
    pwm_table$Intragenic[i] = intragenic
  } else{ #when pwm Position hits on a scaffold without annotation
    intergenic = NA
    five_plus = NA
    five_minus = NA
    dist_to_plus = NA
    dist_to_minus = NA
    intragenic = NA
    pwm_table$Intergenic[i] = intergenic
    pwm_table$FivePrimeGenePlusStrand[i] = five_plus
    pwm_table$FivePrimeGeneMinusStrand[i] = five_minus
    pwm_table$DistancetoFivePlus[i] = dist_to_plus
    pwm_table$DistanceToFiveMinus[i] = dist_to_minus
    pwm_table$Intragenic[i] = intragenic
  }
  #rm(five_plus, five_minus, dist_to_plus, dist_to_minus, intergenic, pwm_row, gff_current_scaf, hits, firstrow, scaf_begin, scaf_end)
}

pwm_table$X = NULL                                                  #Might need this if old row names are added as a column
length(which(is.na(pwm_table$Intergenic)))                          #gives the number of PWM hits nowhere near any annotation
#pwm_table = pwm_table[-which(is.na(pwm_table$Intergenic)),]         #this will remove all rows that didn't receive any annotation
output_file_name = gsub(".csv", "_With-Annotation.csv", input_file) #modify file name
write.csv(pwm_table, file=output_file_name, row.names=FALSE)  
