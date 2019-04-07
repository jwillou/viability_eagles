setwd("/scratch/snyder/j/jwillou/viab_resampling/pf/")
#setwd("/Users/jannawilloughby/Dropbox/viabilityselection/Poolseq/boot_perm/")

#file name prefixes to iterate over
spp = c("ge", "ie", "pf") 

#initiate file
towrite = c(as.character(Sys.time()))
write.table(t(towrite), "results/results.txt", sep="\t", row.names=FALSE, col.names=FALSE, append=FALSE, quote=FALSE) #write time to file
d=3
#for(d in 1:length(spp)){
  ####read in data and set some column names####
  data = read.table(paste("data/", spp[d], "_2alleles_diff_het_freq_min2_worc_splitalleles4.txt", sep=""), sep="\t", header=FALSE)
  colnames(data) = c("scaffold", "location", "nucleotide", "numnucl", "alleles", paste(cbind(as.data.frame(rep("noidea", 14)), as.data.frame(seq(1,14,1)))[,1], cbind(as.data.frame(rep("noidea", 14)), as.data.frame(seq(1,14,1)))[,2], sep=""),"adultH", "juvH", "hetdiff")

  ####resampling by scaffold####
  write.table("", "results/results.txt", sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
  towrite = c(spp[d], "resampling by scaffold results")
  write.table(t(towrite), "results/results.txt", sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE) #write out file label

  #sample 80% of rows 1000 times, then estimate mean diff in H for each time - perhaps sample only one per scaffold?
  boots = NULL
  scaffs = unique(data$scaffold)
  for(r in 1:1000){
    scaff = scaffs[sample(1:length(scaffs), round(length(scaffs)*.8, 0))]
    temp  = NULL
    for(s in 1:length(scaff)){
      tt = data$hetdiff[data$scaffold==as.character(scaff[s])]
      choose = sample(1:length(tt),1)
      tt = tt[choose]
      temp = rbind(temp, tt)
    }
    boots = c(boots, mean(temp, na.rm=TRUE))
  }

  #print out histogram of results (for QC)
  pdf(paste("results/histograms/", spp[d], "_resampling_scaffold.pdf", sep=""), width=4, height=4, useDingbats=FALSE, onefile=TRUE)
  hist(boots, main=paste(spp[d], "resampling means by scaffold"), xlab="diff. in het", ylab="frequency", xlim=c(-0.03, 0.03))
  dev.off()

  #write out mean and 95% CI of diff in H
  towrite = c(paste("mean diff = ", round(mean(boots), 4), sep=""), paste("95% CI = ", round(quantile(boots, probs=c(0.025, 0.975)[1]), 4), ",", round(quantile(boots, probs=c(0.025, 0.975)[2]), 4), sep=""))
  write.table(t(towrite), "results/results.txt", sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)

#}
