#This script will do some preliminrary analysis of high-score hits with their genomic context

input_file = "/home/tlicknac/Desktop/Paramecium_Genome_Data/Data/Motif/pwm-distribution-evalue_7605-1_With-Expression.csv"
pwm_table1 = read.csv(input_file, header=T, as.is = T)

#Look at summary stats for distances to start codon
plus1 = pwm_table1[which(pwm_table1$Strand == "+"),][11]
minus1 = pwm_table1[which(pwm_table1$Strand == "-"),][12]
vdists1 = c(as.numeric(plus1$DistancetoFivePlus), as.numeric(minus1$DistanceToFiveMinus))
summary(vdists1)

#Look at summary stats for PWM scores
summary(pwm_table1$PWM_Score)

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

