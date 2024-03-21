# Population Genetics Exam 2024 - Group 2

##### Quick commands
```unix
# Navigate to group working directory
cd ~/groupdirs/SCIENCE-BIO-popgen_course-project/Group2_ArcticFox
```

---
## Data
We are provided with...
* **`introduction_AF.pdf`**: An introduction to the project and research questions
* **`AF.imputed.thin.{bed,bim,fam}`**: A PLINK dataset of 816284 biallelic SNPs from 47 arctic foxes.
* **`sample_popinfo.tsv`**: Tabular data on the region and country of origin for each sample in `AF.imputed.thin.fam`

##### Quick commands
```unix
# (Re)download project data
cd ~/groupdirs/SCIENCE-BIO-popgen_course-project/Group2_ArcticFox
cp ~/groupdirs/SCIENCE-BIO-Popgen_Course/projects/arctic_fox/* .

# View files
ls -l
```

```
# Fix chromosome identifier in .bim-file
python bim_fix.py 
```

### Characterizing the data
```R
# Loading PLINK genotype data
library(snpMatrix)
data <- read.plink("AF.imputed.thin")
geno <- matrix(as.integer(data@.Data),nrow=nrow(data@.Data))
geno[geno==0] <- NA
geno <- geno-1

# Loading sample data
geo <- read.table("sample_popinfo.tsv",header=1)

# Sorting data according to regions
region_order <- c("Zackenberg", "Scoresbysund", "Kangerlussuaq", "Qanisartuut", "Taymyr", "Bylot_island", "Karrak_lake")
geo$Region <- factor(geo$Region, levels=region_order)
geo <- geo[order(geo$Region), ]
geno <- geno[as.integer(rownames(geo)),]
```

```R
# Checking data dimensions
dim(geno)
dim(fam)
dim(geo)
```
---

## Population substructure
TODO
- Try to subdivide populations by PCAs
    - Kmeans clustering? 3,4,5 k
    - Fst matrix
- Relatedness plot?

### PCA with PLINK
In UNIX:
```unix
plink --bfile AF.imputed.thin --pca 50 --maf 0.15 --geno 0 --out pca/pca_results
```

In R:
```R
# Loading PCA results
eigenvals <- read.table("pca/pca_results.eigenval")
pcs <- read.table("pca/pca_results.eigenvec")[3:23]
names(pcs) <- paste("PC",1:20, sep="")
names(eigenvals) <- paste("PC",1:47, sep="")


# Loading sample data
geo <- read.table("sample_popinfo.tsv",header=1)

# Formatting text
geo$Region <- gsub("_", " ", geo$Region)

# Sorting data according to regions
region_order <- c("Zackenberg", "Scoresbysund", "Kangerlussuaq", "Qanisartuut", "Taymyr", "Bylot island", "Karrak lake")
geo$Region <- factor(geo$Region, levels=region_order)
geo <- geo[order(geo$Region), ]
pcs <- pcs[as.integer(rownames(geo)),]

# Extracting importance of PCs
pca_importance <- eigenvals / sum(eigenvals)
names(pca_importance) <- c("Importance")
pca_importance$PC <- 1:47
```

```R
library(ggplot2)
library(scales)

# Plot PC importance
p <- ggplot(pca_importance, aes(x = PC,y = Importance)) +
  geom_line(col="#7E0605") +
  geom_point(col="#7E0605") +
  labs(x = "Principal component", y = "Explained variance") +
  theme_classic() +
  theme(text = element_text(family = "serif")) + 
  scale_x_continuous(limits = c(0.6, 47.6), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.001, 0.1), expand = c(0, 0), labels=percent_format())

ggsave("pca_importance.png", p, width = 5, height = 3, dpi = 300)
```

```R
# Plot PCs
library(ggplot2)
library(ggpubr)

PC1_explained <- round(pca_importance[1,1]*100, 1)
PC2_explained <- round(pca_importance[2,1]*100, 1)
PC3_explained <- round(pca_importance[3,1]*100, 1)

region_palette <- c("#7E0605", "#C00E38", "#FF0000", "#FF7D7B", "#B82EAA", "#13501B", "#17A238")

p12 <- 
  ggplot(pcs, aes(x=PC1, y=PC2, color=geo$Region)) +
  geom_point() +
  labs(x = paste0('PC1 (', PC1_explained, '% of variance)'),
       y = paste0('PC2 (', PC2_explained, '% of variance)'),
       color = "Sample region") +
  theme_classic() +
  theme(text = element_text(family = "serif")) +
  scale_color_manual(values = region_palette)

p13 <- 
  ggplot(pcs, aes(x=PC1, y=PC3, color=geo$Region)) +
  geom_point() +
  labs(x = paste0('PC1 (', PC1_explained, '% of variance)'),
       y = paste0('PC3 (', PC3_explained, '% of variance)')) +
  theme_classic() +
  theme(text = element_text(family = "serif")) +
  scale_color_manual(values = region_palette)

p23 <- 
  ggplot(pcs, aes(x=PC2, y=PC3, color=geo$Region)) +
  geom_point() +
  labs(x = paste0('PC2 (', PC2_explained, '% of variance)'),
       y = paste0('PC3 (', PC3_explained, '% of variance)')) +
  theme_classic() +
  theme(text = element_text(family = "serif")) +
  scale_color_manual(values = region_palette)

combined_plot <-
  ggarrange(p12, p13, p23, ncol = 3,
            common.legend = TRUE, legend = "right")

ggsave("pca_plink.png", combined_plot, width = 10, height = 3, dpi = 300)
```

### Clustering with PCA from PLINK
```R
# Loading PLINK genotype data
library(snpMatrix)
data <- read.plink("AF.imputed.thin")
geno <- matrix(as.integer(data@.Data),nrow=nrow(data@.Data))
geno[geno==0] <- NA
geno <- geno-1

# Removing SNPs missing data
geno <- geno[,colSums(is.na(geno))==0]

# Loading PCA fit
pcs <- read.table("pca/pca_results.eigenvec")[3:23]
names(pcs) <- paste("PC",1:20, sep="")
```

```R
WC84<-function(x,pop){
     # function to estimate Fst using Weir and Cockerham estimator.
     # x is NxM genotype matrix, pop is N length vector with population assignment for each sample
     # returns list with fst between population per M snps (theta) and other variables
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
     list(F=F,theta=theta,f=f,theta_w=theta_w,a=a,b=b,c=c,total=c+b+a)
}
```

```R
# Assigning samples to clusters
k = 7
pcs_max = 5
clusters <- kmeans(pcs[1:pcs_max], k)$cluster

# Comparing clusters by Fst
cluster_pairs <- t(combn(unique(clusters), 2))
fsts <- apply(cluster_pairs, 1, function(pair) WC84(geno[clusters %in% pair,], clusters[clusters %in% pair]))
names(fsts) <- apply(cluster_pairs, 1, paste, collapse="_")
lapply(fsts, function(x) x$theta_w)
```