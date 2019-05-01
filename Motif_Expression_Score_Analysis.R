#This script will do the bulk of the work trying to relate expression data to DNA motifs discovered through MEME

intergenic_file = "/home/tlicknac/Desktop/Paramecium_Genome_Data/Data/WGD06719_intergenic-expression_1.csv"
pwm_table1 = read.csv(intergenic_file, header=T, as.is = T)

#Look at summary stats
nrow(pwm_table1)
plus1 = pwm_table1[which(pwm_table1$Strand == "+"),][11]
minus1 = pwm_table1[which(pwm_table1$Strand == "-"),][12]
vdists1 = c(as.numeric(plus1$DistancetoFivePlus), as.numeric(minus1$DistanceToFiveMinus))
summary(vdists1)

#Look at summary stats for PWM scores
summary(pwm_table1$PWM_Score)

#Means differences between motif absence/presence
rnaseq = read.table("")

#Correlation between score and distance
df_score_dist_plus = data.frame(matrix(ncol = 2))
df_count = 1
for(i in which(pwm_table1$Strand == "+")){
  score_dist = c(pwm_table1$PWM_Score[i], pwm_table1$DistancetoFivePlus[i])
  df_score_dist_plus[df_count,] = score_dist
  df_count = df_count+1
}
df_score_dist_plus = df_score_dist_plus[-which(is.na(df_score_dist_plus$X2)),]

df_score_dist_minus = data.frame(matrix(ncol = 2))
df_count = 1
for(i in which(pwm_table1$Strand == "-")){
  score_dist = c(pwm_table1$PWM_Score[i], pwm_table1$DistanceToFiveMinus[i])
  df_score_dist_minus[df_count,] = score_dist
  df_count = df_count+1
}
df_score_dist_minus = df_score_dist_minus[-which(is.na(df_score_dist_minus$X2)),]

df_score_dist = rbind(df_score_dist_plus, df_score_dist_minus)
cor.test(x=df_score_dist$X1, y=df_score_dist$X2)

#Correlation between distance and expression




#change this when have time...  
pwm_tab_RNAseq = read.csv("/home/tlicknac/Desktop/Paramecium_Genome_Data/Data/pwm-distribution-duplicate_15708-1_With-CellCycleSignificance.csv", header=T, as.is=T)
five_plus = pwm_tab_RNAseq[which(is.na(pwm_tab_RNAseq$CellCycle_Significance_FivePrimePlus) == F),]
five_minus = pwm_tab_RNAseq[which(is.na(pwm_tab_RNAseq$CellCycle_Significance_FivePrimeMinus) == F),]

yes_plus = five_plus[which(five_plus$CellCycle_Significance_FivePrimePlus == "YES"),]        #get rows where the plus-strand gene is DE
yes_minus = five_minus[which(five_minus$CellCycle_Significance_FivePrimeMinus == "YES"),]    #same for minus-strand gene
no_plus = five_plus[which(five_plus$CellCycle_Significance_FivePrimePlus == "NO"),]
no_minus = five_minus[which(five_minus$CellCycle_Significance_FivePrimeMinus == "NO"),]

yes_combined = rbind(yes_plus, yes_minus)  #combine for purposes of counting
no_combined = rbind(no_plus, no_minus)

nrow(yes_combined)        #significantly expressed during cell cycle
nrow(no_combined)         #significantly expressed during cell cycle

pwm_scores_yes = c(yes_plus$PWM_Score[which(yes_plus$Strand == "+")], yes_minus$PWM_Score[which(yes_minus$Strand == "-")])  #vector of pwm scores for all with PWM on correct strand
pwm_scores_no = c(no_plus$PWM_Score[which(no_plus$Strand == "+")], no_minus$PWM_Score[which(no_minus$Strand == "-")])

summary(pwm_scores_yes)
summary(pwm_scores_no)
