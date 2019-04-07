setwd("~/Desktop/simviab")
library(RColorBrewer)
write.table(t(c("replicate", "alleleFreqA", "AAmortailityFreq", "BBmortailityFreq", "adultHet", "chickHet_premort", "chickHet_postmort", "hetDiff")), "diroutput.csv", sep=",", col.names=F, row.names=F, append=F)

####set parameters####
Nadults   = 1000                  #number of adults
Noff      = 1000                  #number of offspring produced by adults
af        = seq(0.05, 0.95, 0.01) #allele frequency of allele A
mortA     = seq(0.00, 0.95, 0.01) #mortality freq in AA individuals
mortBdiff = seq(0, 0.45, 0.15)    #reduction in BB mortality freq, relative to AA
reps      = 1                     #number of replicates

####run simulation####
for(AF in 1:length(af)){
  for(MA in 1:length(mortA)){
    for(MB in 1:length(mortBdiff)){
      if((mortA[MA]-mortBdiff[MB])<0){next} #is mortaility probability for BB is positive?
      for(r in 1:reps){
        ####assign alleles to adults in HWE####
        #generate genotypes for starting population
        adults   = matrix(nrow=Nadults, ncol=2)
        columns  = seq(1,(2),2)
        
        p = af[AF]
        #create pool of genotypes in HWE and assign to block of adults
        pool = c(rep(2, round(Nadults*p*p, 0)),                                                #homozygous A
                 rep(0, round(Nadults*(1-p)*(1-p), 0)),                                        #homozygous B
                 rep(1, Nadults-(round(Nadults*p*p, 0)+(round(Nadults*(1-p)*(1-p), 0)))))      #heterozygotes AB
        adults = cbind(adults, pool)
        
        #assign individual genotype alleles based on determined genotype
        adults[adults[,3]==2, 1:2] = 1
        adults[adults[,3]==0, 1:2] = 0
        adults[adults[,3]==1, 1] = 1
        adults[adults[,3]==1, 2] = 0
        
        #randomly pair adults (with replacement) until 1000 offspring are created
        adults = cbind(adults, seq(1, nrow(adults)))
        adults = cbind(adults, sample(c("m", "f"), nrow(adults), replace=T))
        madults = adults[adults[,5]=="m",,drop=FALSE]
        fadults = adults[adults[,5]=="f",,drop=FALSE]
        
        #create chicks with parentally inhereited genotypes
        chicks = NULL
        for(i in 1:1000){
          m = madults[sample(c(1:nrow(madults)), 1),,drop=FALSE]
          f = fadults[sample(c(1:nrow(fadults)), 1),,drop=FALSE]
          chicks = rbind(chicks, c(sample(m[,1:2], 1), sample(f[,1:2], 1)))
        }
        chicks = cbind(chicks, (as.numeric(as.character(chicks[,1])) + as.numeric(as.character(chicks[,2]))))
        
        #apply mortaility to homozygous AA individuals
        chicks = cbind(chicks, rep(1, nrow(chicks))) #add alive indicator column
        chicks[chicks[,3]==2, 4] = sample(c(0,1), nrow(chicks[chicks[,3]==2,,drop=FALSE]), prob=c(mortA[MA], (1-(mortA[MA]))), replace=T)
        
        #apply mortaility to heterozygous AB individuals
        chicks[chicks[,3]==1, 4] = sample(c(0,1), nrow(chicks[chicks[,3]==1,,drop=FALSE]), prob=c(mortA[MA]-mortBdiff[MB], (1-(mortA[MA]-mortBdiff[MB]))), replace=T)
        
        #estimate output parameters 
        #adultHet          = nrow(adults[adults[,3]==1, ,drop=FALSE])/nrow(adults)  #heterozygosity of adults
        #chickHet_premort  = nrow(chicks[chicks[,3]==1, ,drop=FALSE])/nrow(chicks)  #heterozygosity of offspring BEFORE viability selection
        #chickHet_postmort = nrow(chicks[chicks[,3]==1 & chicks[,4]==1, ,drop=FALSE])/nrow(chicks[chicks[,4]==1, ,drop=FALSE])  #heterozygosity of offspring AFTER viability selection
        #estimated het as if from pooled seq data
        adultHet          = 2*sum(as.numeric(adults[,1]), as.numeric(adults[,2]))/(nrow(adults)*2)*(1-sum(as.numeric(adults[,1]), as.numeric(adults[,2]))/(nrow(adults)*2))   #heterozygosity of adults 
        chickHet_premort  = 2*sum(as.numeric(chicks[,1]), as.numeric(chicks[,2]))/(nrow(chicks)*2)*(1-sum(as.numeric(chicks[,1]), as.numeric(chicks[,2]))/(nrow(chicks)*2))  #heterozygosity of offspring BEFORE viability selection
        chickHet_postmort = 2*sum(as.numeric(chicks[chicks[,4]==1,1]), as.numeric(chicks[chicks[,4]==1,2]))/(nrow(chicks[chicks[,4]==1,])*2)*(1-sum(as.numeric(chicks[chicks[,4]==1,1]), as.numeric(chicks[chicks[,4]==1,2]))/(nrow(chicks[chicks[,4]==1,])*2))  #heterozygosity of offspring BEFORE viability selection#heterozygosity of offspring AFTER viability selection
        
        #output: replicate, alleleFreqA, AAmortailityFreq, BBmortailityFreq, adultHet, chickHet_premort, chickHet_postmort, hetDiff
        writeout = c(r, af[AF], mortA[MA], mortA[MA]-mortBdiff[MB], adultHet, chickHet_premort, chickHet_postmort, (adultHet-chickHet_postmort))
        write.table(t(writeout), "diroutput.csv", sep=",", col.names=F, row.names=F, append=T)
        
      }
    }
  }
}

####create heatmaps####
alldata = read.table("diroutput.csv", header=T, sep=",")

minHD = -0.4
maxHD = 0.4
numb  = length(seq(minHD, maxHD, 0.1)) - 1

####BB=AA####
data = alldata[alldata$AAmortailityFreq==alldata$BBmortailityFreq,]
data = cbind(data$alleleFreqA, data$AAmortailityFreq, data$hetDiff)

#reorder for heatmap function
afs = unique(data[,1])
mor = unique(data[,2])
out = data.frame(mor = mor)
for(a in 1:length(afs)){
  as = data[data[,1]==afs[a],,drop=F]
  as = data.frame(mor = as[,2], hetDiff = as[,3])
  out = merge(x=out, y=as, by="mor")
}
rownames(out) = out$mor
out$mor = NULL
colnames(out) = paste("afs", afs, sep="_")
data = as.matrix(out)

#plot
pdf(file="dirheatmapBB=AA.pdf", width = 7, height = 7)
heatmap(data, Colv=NA, Rowv=NA, breaks=seq(minHD, maxHD, 0.1), col=brewer.pal(numb, "YlOrRd"), scale="none", bg=brewer.pal(numb, "YlOrRd"))
dev.off()

####BB=AA-0.15####
data = alldata[alldata$BBmortailityFreq==(alldata$AAmortailityFreq-0.15),]
data = cbind(data$alleleFreqA, data$AAmortailityFreq, data$hetDiff)

#reorder for heatmap function
afs = unique(data[,1])
mor = unique(data[,2])
out = data.frame(mor = mor)
for(a in 1:length(afs)){
  as = data[data[,1]==afs[a],,drop=F]
  as = data.frame(mor = as[,2], hetDiff = as[,3])
  out = merge(x=out, y=as, by="mor")
}
rownames(out) = out$mor
out$mor = NULL
colnames(out) = paste("afs", afs, sep="_")
data = as.matrix(out)

#plot
pdf(file="dirheatmapBB=AA-0.15.pdf", width = 7, height = 7)
heatmap(data, Colv=NA, Rowv=NA, breaks=seq(minHD, maxHD, 0.1), col=brewer.pal(numb, "YlOrRd"), scale="none", bg=brewer.pal(numb, "YlOrRd"))
dev.off()

####BB=AA-0.30####
data = alldata[alldata$BBmortailityFreq==(alldata$AAmortailityFreq-0.30),]
data = cbind(data$alleleFreqA, data$AAmortailityFreq, data$hetDiff)

#reorder for heatmap function
afs = unique(data[,1])
mor = unique(data[,2])
out = data.frame(mor = mor)
for(a in 1:length(afs)){
  as = data[data[,1]==afs[a],,drop=F]
  as = data.frame(mor = as[,2], hetDiff = as[,3])
  out = merge(x=out, y=as, by="mor")
}
rownames(out) = out$mor
out$mor = NULL
colnames(out) = paste("afs", afs, sep="_")
data = as.matrix(out)

#plot
pdf(file="dirheatmapBB=AA-0.30.pdf", width = 7, height = 7)
heatmap(data, Colv=NA, Rowv=NA, breaks=seq(minHD, maxHD, 0.1), col=brewer.pal(numb, "YlOrRd"), scale="none", bg=brewer.pal(numb, "YlOrRd"))
dev.off()

####BB=AA-0.45####
data = alldata[alldata$BBmortailityFreq==(alldata$AAmortailityFreq-0.45),]
data = cbind(data$alleleFreqA, data$AAmortailityFreq, data$hetDiff)

#reorder for heatmap function
afs = unique(data[,1])
mor = unique(data[,2])
out = data.frame(mor = mor)
for(a in 1:length(afs)){
  as = data[data[,1]==afs[a],,drop=F]
  as = data.frame(mor = as[,2], hetDiff = as[,3])
  out = merge(x=out, y=as, by="mor")
}
rownames(out) = out$mor
out$mor = NULL
colnames(out) = paste("afs", afs, sep="_")
data = as.matrix(out)

#plot
pdf(file="dirheatmapBB=AA-0.45.pdf", width = 7, height = 7)
heatmap(data, Colv=NA, Rowv=NA, breaks=seq(minHD, maxHD, 0.1), col=brewer.pal(numb, "YlOrRd"), scale="none", bg=brewer.pal(numb, "YlOrRd"))
dev.off()
#

