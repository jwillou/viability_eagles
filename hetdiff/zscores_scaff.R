setwd("PoolSeq raptors/files to try with R scale/July_2018/")

#read in data and set some column names
data = read.table("pf_2alleles_diff_het_freq_min2_worc_splitalleles4.txt", sep="\t", header=FALSE)
colnames(data) = c("scaffold", "location", "nucleotide", "numnucl", "alleles", paste(cbind(as.data.frame(rep("noidea", 16)), as.data.frame(seq(1,16,1)))[,1], cbind(as.data.frame(rep("noidea", 16)), as.data.frame(seq(1,16,1)))[,2], sep=""), "hetdiff")

#create new column for z scores
data$zhetdiff = rep(NA, nrow(data))

#create df for storing means/sd for QC
sumdata = data.frame(scaffold    = rep(NA, length(unique(data$scaffold))), 
                     mean        = rep(NA, length(unique(data$scaffold))), 
                     sd          = rep(NA, length(unique(data$scaffold))),
                     nsnps       = rep(NA, length(unique(data$scaffold))),
                     meanhetdiff = rep(NA, length(unique(data$scaffold))),
                     length      = rep(NA, length(unique(data$scaffold))))

#find all scaffold names
scafs = unique(data$scaffold)

#set up OUTPUT df
OUTPUT = NULL

#iterate over each scaffold, generate z scores and record means/sd
for(s in 1:length(unique(data$scaffold))){
  temp = data[data$scaffold==as.character(scafs[s]),,drop=FALSE]
  temp$zhetdiff = scale(temp$hetdiff, center = TRUE, scale = TRUE)
  OUTPUT = rbind(OUTPUT, temp)
  sumdata$scaffold[s]    = as.character(scafs[s])
  sumdata$mean[s]        = mean(temp$zhetdiff, na.rm=TRUE)
  sumdata$sd[s]          = sd(temp$zhetdiff, na.rm=TRUE)
  sumdata$nsnps[s]       = nrow(temp)
  sumdata$meanhetdiff[s] = mean(temp$hetdiff)
  sumdata$length[s]      = temp$location[nrow(temp)] - temp$location[1]
}
data = OUTPUT

#find >3sd
highhetdiff = data[data$zhetdiff>=3,,drop=FALSE]
nrow(highhetdiff)
