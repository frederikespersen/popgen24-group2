# Loading PLINK genotype data
library(snpMatrix)
data <- read.plink("AF.imputed.thin")
geno <- matrix(as.integer(data@.Data),nrow=nrow(data@.Data))
geno[geno==0] <- NA
geno <- geno-1
# keep only SNPs without missing data
geno <- geno[,complete.cases(t(geno))]

# Loading clustering definitions
clusters <- read.csv("clusters.csv", header=1, row.names=1)

WC84<-function(x,pop){
  # function to estimate Fst using Weir and Cockerham estimator.
  # x is NxM genotype matrix, pop is N length vector with population assignment for each sample
  # returns list with fst between population per M snps (theta) and other variables
  start.time = Sys.time()
  ###number ind in each population
  n<-table(pop)
  ###number of populations
  npop<-nrow(n)
  ###average sample size of each population
  n_avg<-mean(n)
  ###total number of samples
  N<-length(pop)
  ###frequency in samples
  p<-apply(x,2,function(x,pop){tapply(x,pop,mean)/2},pop=pop)
  ###average frequency in all samples (apply(x,2,mean)/2)
  p_avg<-as.vector(n%*%p/N )
  ###the sample variance of allele 1 over populations
  s2<-1/(npop-1)*(apply(p,1,function(x){((x-p_avg)^2)})%*%n)/n_avg
  ###average heterozygotes
  # h<-apply(x==1,2,function(x,pop)tapply(x,pop,mean),pop=pop)
  #average heterozygote frequency for allele 1
  # h_avg<-as.vector(n%*%h/N)
  ###faster version than above:
  h_avg<-apply(x==1,2,sum)/N
  ###nc (see page 1360 in wier and cockerhamm, 1984)
  n_c<-1/(npop-1)*(N-sum(n^2)/N)
  ###variance between populations
  a <-n_avg/n_c*(s2-(p_avg*(1-p_avg)-(npop-1)*s2/npop-h_avg/4)/(n_avg-1))
  ###variance between individuals within populations
  b <- n_avg/(n_avg-1)*(p_avg*(1-p_avg)-(npop-1)*s2/npop-(2*n_avg-1)*h_avg/(4*n_avg))
  ###variance within individuals
  c <- h_avg/2
  ###inbreeding (F_it)
  F <- 1-c/(a+b+c)
  ###(F_st)
  theta <- a/(a+b+c)
  ###(F_is)
  f <- 1-c(b+c)
  ###weighted average of theta
  theta_w<-sum(a)/sum(a+b+c)
  print(theta_w)
  print(unique(pop))
  return(theta_w)
}

# Initialising dataframe for storage
df <- data.frame(C=character(), pair=character(), Fst=numeric())

# Comparing clusters by Fst
for (C in names(clusters)[-1]) {
  print(C)
  C_clusters <- clusters[[C]]
  cluster_pairs <- t(combn(unique(C_clusters), 2))
  Fsts <- apply(cluster_pairs, 1, function(pair) WC84(geno[C_clusters %in% pair,], C_clusters[C_clusters %in% pair]))
  df <- rbind(df,
              data.frame(
                C=C,
                pair=apply(cluster_pairs, 1, paste, collapse=" "),
                Fst=Fsts))
  write.csv(df, "fsts.csv")              
}