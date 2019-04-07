setwd("/scratch/snyder/j/jwillou/hetbygene")
#setwd("/Users/jannawilloughby/Dropbox/viabilityselection/Poolseq/genehetdif")

#iterate over data file for each species
species = c("ge", "ie", "pf")

for(sp in 1:length(species)){
  data = read.table(paste(species[sp], "_het_all_loci_by_gene.txt", sep=""), sep="\t", header=TRUE)
  colnames(data) = c("scaffoldID", "SNPlocation", "SNP", "hetdiff", "scaffoldID2", "geneSpos", "geneEpos", "gene")
  
  #initiate output file
  write.table(t(c("gene", "Mhetdiff", "SDhetdiff", "nsnps", "quant2.5", "quant97.5")), paste("summarized_", species[sp], "_het_all_loci_by_gene.txt", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, append=FALSE)
  
  #get list of unique genes
  genes = unique(data$gene)
  
  #iterate over each gene, record means/sd/nu/CI
  for(r in 1:length(genes)){
    #create new data frame for prot summaries
    gene = data.frame(gene = 0,          #gene name
                      Mhetdiff = 0,      #mean difference in het within each gene
                      SDhetdiff = 0,     #SD in difference in het within each gene
                      nsnps = 0,         #number of SNPs in each gene
                      quant2.5 = 0,      #lower bound of 95% CI of difference in het within each gene
                      quant97.5 = 0)     #upper bound of 95% CI of difference in het within each gene
    
    #find all snps within the gene
    temp = data[data$gene==as.character(genes[r]),,drop=FALSE]
    
    #record information
    gene$gene      = genes[r]
    gene$Mhetdiff  = mean(temp$hetdiff, na.rm=TRUE)
    gene$SDhetdiff = sd(temp$hetdiff, na.rm=TRUE)
    gene$nsnps     = nrow(temp)
    gene$quant2.5  = quantile(temp$hetdiff, probs=0.025)
    gene$quant97.5 = quantile(temp$hetdiff,  probs=0.975)
    write.table(gene, paste("summarized_", species[sp], "_het_all_loci_by_gene.txt", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE)
  }
}

