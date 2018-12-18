#DONE

#This script will use a PSSM to return scores for each position relating how well the PSSM matches with the proceeding sequence
  #Iterate through the + and - strand of each scaffold for each assembly
  #Returns vector of scores for each scaffold, of which only those above some threshold will be kept
  #Generates a data frame with information for each high-score hit, exports as a .csv

#INPUT:
  #C++ code- getPWM_Scores.cpp
  #PSSMs
  #FASTA files for analysis
#OUTPUT:
#Dataframe with consensus motif, species, scaffold, position, score ONLY for positions where score>0


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
fasta_directory = "/N/u/tlicknac/Carbonate/Paramecium_FASTA/"                                  #location of FASTAs                                                                             #Rm all files that arent fasta files *this will capture .fasta and .fa
FASTAs = c("ptetraurelia_mac_51.fa", "biaurelia_V1-4_assembly_v1.fa", "sexaurelia_AZ8-4_assembly_v1.fasta")

pwm_directory = "/N/dc2/scratch/tlicknac/12nt_PWM_All/"                                        #location of PWMs

PWMs = c("pbi_1849.fasta|1.pwm")  #CHANGE
#PWMs = c("high_expression_1.pwm")
motif_type = "tss"                #CHANGE
#motif_type = "histone"
motif_number = "1"                #CHANGE

firstloop=TRUE

for(assembly in FASTAs){                                                                                              #Iterate through fasta file names
  cat("Starting with species: ", assembly, "\n")
  a = paste(fasta_directory, assembly, sep="")                                 
  fasta = read.fasta(a, as.string=T)                                                                                  #read in fasta 1 at a time... must preserve memory
  
  for(scaf in fasta){                                                                                                 #Iterate through each scaffold
    cat("Onto scaffold: ", getName(scaf), " of species: ", assembly, "\n")
    scaf_seq = as.character(scaf[1])                                                                                  #get sequence
    reverse_seq = c2s(rev((comp(s2c(scaf_seq)))))
    
    for(pwm in PWMs){
      pwm_file_path = paste(pwm_directory, pwm, sep="") 
      pwm_tab = as.matrix(read.table(pwm_file_path, header=T, as.is=T, sep="\t"))
      vscores = getScores(pwm_tab, scaf_seq) #C++ code
      vscores_rev = getScores(pwm_tab, reverse_seq)
      vgood = which(vscores>0)                                #POSITION ####CHANGE CUTOFF IF NEEDED. 0 CHOSEN AS DEFAULT... only a handful of positions will meet this criterion####
      vgood_rev = which(vscores_rev>0)                                
      
      if(length(vgood)>0){
        
        for(position in vgood){
          
          if(firstloop==TRUE){
            finaldf = data.frame( t(matrix( c( consensusFromPWM(pwm_tab), pwm, assembly, getName(scaf), position, vscores[position],"+" ), byrow=T )) )
            colnames(finaldf) = c("Consensus_Motif", "POFF_Row", "Species", "Scaffold", "Position", "PWM_Score", "Strand")
            firstloop=FALSE
            
          } else{
            tmpdf = t(data.frame(c( consensusFromPWM(pwm_tab), pwm, assembly, getName(scaf), position, vscores[position], "+" )))  #all the info we need
            colnames(tmpdf) = c("Consensus_Motif", "POFF_Row", "Species", "Scaffold", "Position", "PWM_Score", "Strand")
            finaldf = rbind(finaldf, tmpdf)
            cat("Wrote to finaldf", "\n")
          }
        }
      }
      
      if(length(vgood_rev>0)){
        
        for(rposition in vgood_rev){
          
          if(firstloop == T){
            rposition = nchar(reverse_seq)-rposition+1    #this is needed to get position correct.hit on reverse comp isnt the position of actually scaf
            finaldf = data.frame( t(matrix( c( consensusFromPWM(pwm_tab), pwm, assembly, getName(scaf), rposition, vscores_rev[rposition],"+" ), byrow=T )) )
            colnames(finaldf) = c("Consensus_Motif", "POFF_Row", "Species", "Scaffold", "Position", "PWM_Score", "Strand")
            firstloop=FALSE
            
          } else {
            tmpdf = t(data.frame(c( consensusFromPWM(pwm_tab), pwm, assembly, getName(scaf), rposition, vscores_rev[rposition], "-" )))  #all the info we need
            colnames(tmpdf) = c("Consensus_Motif", "POFF_Row", "Species", "Scaffold", "Position", "PWM_Score", "Strand")
            finaldf = rbind(finaldf, tmpdf)
          }
        }
      }
    } 
  }
}
rownames(finaldf) = NULL

txt = "pwm-distribution-"
file = strsplit(strsplit(PWMs[1], "_")[[1]][2], ".fasta")[[1]][1]
filename = paste(txt, motif_type, "_", file, "-", motif_number, ".csv", sep="")   #use motif_type and motif_number 
write.csv(finaldf, file=filename)


#Create data frame with 3 columns
  #scaffold number, position, intergenic?
