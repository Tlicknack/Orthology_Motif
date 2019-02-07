#FOr motif File without annotation

setwd("~/Desktop/Paramecium_Genome_Data/Data/")

dfmotif = read.csv("Newest_MEME_Final_File.csv", header=T)
hist(as.numeric(format(dfmotif$EvalueLog, scientific=T)), base=10, breaks=200, main="Distribution of E-values for All Motifs", xlab="log(E-value)")
summary(nchar(as.character(dfmotif$CommonVariant)))

#Get rid of obviously bad motifs
dfmotif = dfmotif[-which(dfmotif$EvalueLog > 0.05),]

#AT content of promoters
which(dfmotif$BackgroundAT < 0.624)  #this is the lowest we can go... ~62.5% AT 

#Low E values
summary(dfmotif$EvalueLog)
upper_25th = dfmotif[which(dfmotif$EvalueLog < as.numeric(summary(dfmotif$EvalueLog)[4])),]
nrow(upper_25th)
hist(as.numeric(format(upper_25th$EvalueLog, scientific=T)), base=10, breaks=200, main="Distribution of E-values for Most Significant Motifs", xlab="log(E-value)")
summary(nchar(as.character(upper_25th$CommonVariant)))


#Single variant motifs
single_variant = dfmotif[which(dfmotif$MotifMismatch == 0),]
nrow(single_variant)
summary(single_variant$EvalueLog)
summary(nchar(as.character(single_variant$CommonVariant)))
hist(log(as.numeric(format(single_variant$EValue, scientific=T)), base=10), breaks=200, main="Distribution of E-values for Single Variant Motifs", xlab="log(E-value)")


#Single variant and high nseqs
single_variant_highseqs = single_variant[which(single_variant$Nseqs>10),]
summary(nchar(as.character(single_variant_highseqs$CommonVariant)))


#Single variant, high nseqs, single variant flank(s)
dfhighly_conserved = single_variant_highseqs[which(single_variant_highseqs$LeftMismatch == 0 | single_variant_highseqs$RightMismatch == 0),]  #126
dfultra_conserved = single_variant_highseqs[which(single_variant_highseqs$LeftMismatch == 0 & single_variant_highseqs$RightMismatch == 0),]   #33


#Motifs found multiple times
dup_motifs = dfmotif[which(duplicated(dfmotif$Consensus)),]
dup_kmers = dfmotif[which(duplicated(dfmotif$Kmer) & nchar(as.character(dfmotif$Kmer)) > 6),]
median(as.numeric(format(dup_motifs$EValue, scientific=T)))
hist(log(as.numeric(format(dup_motifs$EValue, scientific=T)), base=10), breaks=100, main="Distribution of E-values for Duplicated Motifs", xlab="log(E-value)")

duplicates = dup_kmers[which(dup_kmers$File %in% dup_motifs$File),]

#Motifs close to start codon
utr_motifs = dfmotif[which(dfmotif$DistanceToStart < 15),]
nrow(utr_motifs)
median(utr_motifs$EValue)
hist(log(as.numeric(format(utr_motifs$EValue, scientific=T)), base=10), breaks=100, main="Distribution of E-values for UTR Motifs", xlab="log(E-value)")
dfmotif[which(dfmotif$DistanceToStart == 12),]
twleve = (dfmotif[which(dfmotif$DistanceToStart == 12),])


