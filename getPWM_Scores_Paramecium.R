#DONE
#This script will take PSSM (written as PWMs in the script) for potential TFBSs discovered through sequence conservation
#This will apply that PSSM to each position in the genome, one scaffold at a time (- strand too) and return a data.frame/.csv with those coordinates


library(Rcpp)
sourceCpp("/N/u/tlicknac/Carbonate/getPWM_Scores.cpp")
library(seqinr)


consensusFromPWM = function(pwm_tab){
  motif = c()
  i=1
  max_cols = max.col(pwm_tab, ties.method="first")
  
  for(biggest in max_cols){
    i=i+1
    motif[i] = names(pwm_tab[1,biggest]) 
  }
  motif = paste(motif, collapse="")
  motif = substr(motif, 3, nchar(motif))
  return(motif)
}

#MAIN
fasta_directory = "/N/u/tlicknac/Carbonate/Paramecium_FASTA/"                                                         #location of FASTAs                                                                             #Rm all files that arent fasta files *this will capture .fasta and .fa
FASTAs = c("ptetraurelia_mac_51.fa", "biaurelia_V1-4_assembly_v1.fa", "sexaurelia_AZ8-4_assembly_v1.fasta")

pwm_directory = "/N/dc2/scratch/tlicknac/Newest_All-Aurelias_PWM/"

PWMs = c("WGD13807_promoter.fasta|1.pwm")        

firstloop=TRUE

for(assembly in FASTAs){                                                                                              #Iterate through fasta file names
  cat("Starting with species: ", assembly, "\n")
  aaaa = paste(fasta_directory, assembly, sep="")                                 
  fasta = read.fasta(aaaa, as.string=T, forceDNAtolower = T)                                                          #read in fasta 1 at a time... must preserve memory
  
  for(scaf in fasta){  #scaf = fasta["scaffold51_471"]                                                                #Iterate through each scaffold
    cat("Onto scaffold: ", getName(scaf), " of species: ", assembly, "\n")
    scaf_seq = as.character(scaf[1])                                                                                  #get sequence
    reverse_seq = c2s( rev( comp( s2c(scaf_seq) ) ) )
    
    for(pwm in PWMs){
      pwm_file_path = paste(pwm_directory, pwm, sep="") 
      pwm_tab = as.matrix(read.table(pwm_file_path, header=T, as.is=T, sep="\t"))
      
      motif_length = nchar(as.character(consensusFromPWM(pwm_tab)))
      vscores = getScores(pwm_tab, scaf_seq)                                          #C++ code
      vscores = head(vscores, -(motif_length))                                        #remove the last few scores on the scaffold, as they will be abornally high due to algorithm
      vscores_rev = getScores(pwm_tab, reverse_seq)
      vscores_rev = head(vscores_rev, -(motif_length))  #??will this throw off my reverse positions by motif_length??
      
      vgood = which(vscores>0)                      
      vgood_rev = which(vscores_rev>0)                                
      
      if(length(vgood)>0){
        
        for(position in vgood){
          
          if(firstloop==TRUE){
            finaldf = data.frame( t(matrix( c( substr(scaf_seq, position, position+12), pwm, assembly, getName(scaf), position, vscores[position],"+" ), byrow=T )) )
            colnames(finaldf) = c("Motif", "WGD", "Species", "Scaffold", "Position", "PWM_Score", "Strand")
            firstloop=FALSE
            
          } else{
            tmpdf = t(data.frame(c( substr(scaf_seq, position, position+12), pwm, assembly, getName(scaf), position, vscores[position], "+" )))
            colnames(tmpdf) = c("Motif", "WGD", "Species", "Scaffold", "Position", "PWM_Score", "Strand")
            finaldf = rbind(finaldf, tmpdf)
            cat("Wrote to finaldf", "\n")
          }
        }
      }
      
      if(length(vgood_rev)>0){
        
        for(rposition in vgood_rev){
          
          if(firstloop == T){
            real_pos = (nchar(reverse_seq) - rposition) + 1             #this is required to return the actual coordinates for hits on the + strand... rposition is with respect to - strand
            finaldf = data.frame( t(matrix( c( substr(reverse_seq, rposition, rposition+12), pwm, assembly, getName(scaf), real_pos, vscores_rev[rposition], "-" ), byrow=T )) )
            colnames(finaldf) = c("Motif", "WGD", "Species", "Scaffold", "Position", "PWM_Score", "Strand")
            firstloop=FALSE
            
          } else {
            real_pos = (nchar(reverse_seq) - rposition) + 1
            tmpdf = t(data.frame(c( substr(reverse_seq, rposition, rposition+12), pwm, assembly, getName(scaf), real_pos, vscores_rev[rposition], "-" )))  #all the info we need
            colnames(tmpdf) = c("Motif", "WGD", "Species", "Scaffold", "Position", "PWM_Score", "Strand")
            finaldf = rbind(finaldf, tmpdf)
          }
        }
      }
    } 
  }
}
rownames(finaldf) = NULL

wgd_name = strsplit(PWMs[1], "_")[[1]][1]
file_number = strsplit(strsplit(PWMs[1], "_")[[1]][2], "|")[[1]][16]
filename =  paste("pwm-distribution_", wgd_name, "|", file_number, ".csv", sep="")
write.csv(finaldf, file=filename)
