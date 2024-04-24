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
```

```R
# Checking data dimensions
dim(geno)
dim(fam)
dim(geo)
```
---

## Population substructure

### PCA with PLINK
In UNIX:
```unix
mkdir pca
plink --bfile AF.imputed.thin --pca 47 --geno 0 --out pca/pca_results
```

In R:
```R
# Loading PCA results
eigenvals <- read.table("pca/pca_results.eigenval")
pcs <- read.table("pca/pca_results.eigenvec")[3:49]
names(pcs) <- paste("PC",1:47, sep="")

# Extracting importance / explained variance of PCs
pca_importance <- eigenvals / sum(eigenvals)
names(pca_importance) <- c("Importance")
pca_importance$PC <- 1:47

# Loading sample data and formatting text
geo <- read.table("sample_popinfo.tsv",header=1)
geo$Region <- gsub("_", " ", geo$Region)

# Sorting data according to regions
region_order <- c("Zackenberg", "Scoresbysund", "Kangerlussuaq", "Qanisartuut", "Taymyr", "Bylot island", "Karrak lake")
geo$Region <- factor(geo$Region, levels=region_order)
geo <- geo[order(geo$Region), ]
pcs <- pcs[as.integer(rownames(geo)),]
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

ggsave("pca_plot.png", combined_plot, width = 10, height = 3, dpi = 300)
```

### Clustering with PCA from PLINK
```R
# Loading PCA fit
pcs <- read.table("pca/pca_results.eigenvec")[3:49]
names(pcs) <- paste("PC",1:47, sep="")

# Loading sample data and formatting text
geo <- read.table("sample_popinfo.tsv",header=1)
geo$Region <- gsub("_", " ", geo$Region)
names(pcs) <- geo$Sample

# Sorting data according to regions
clustering_order <- c(35, 22, 13, 39, 3, 42, 43, 41, 45, 44, 1, 5, 46, 47, 23, 14, 20, 36, 29, 17, 26, 11, 4, 40, 8, 32, 6, 33, 30, 27, 21, 37, 9, 24, 15, 18, 16, 38, 19, 34, 7, 31, 25, 28, 12, 2, 10)
```

```R
# Scanning over m and k
num_samples <- length(pcs)
clustering_frequency <- matrix(0, nrow = num_samples, ncol = num_samples)
for (replicate in 1000) {
  for (m in 1:6) {
    for (k in 2:8) {
      clustering <- kmeans(pcs[clustering_order,][1:m], k)$cluster
      clustering_frequency <- clustering_frequency + outer(clustering, clustering, "==")
    }
  }
}
clustering_frequency <- clustering_frequency / max(clustering_frequency)

# Plotting clustering frequency
region_palette <- c("#7E0605", "#B82EAA", "#17A238", "#13501B", "#C00E38", "#FF0000", "#FF7D7B")
region_colours <- region_palette[as.numeric(factor(geo$Region[clustering_order], levels=unique(geo$Region[clustering_order])))]
par(cex = 1, family = "serif")
heatmap(clustering_frequency, Rowv=NA, Colv=NA, revC=TRUE, scale="none", RowSideColors=region_colours, ColSideColors=region_colours, labRow=geo$Sample[clustering_order], labCol=geo$Sample[clustering_order])
```

```R
# Assigning clusters manually
clusters <- data.frame(
  row.names = row.names(geo[clustering_order,]),
  sample = geo$Sample[clustering_order],
  K = c(rep(1,4), rep(2,13), rep(3,9), rep(4,10), rep(5,5), rep(6,6))
)
clusters <- clusters[order(as.numeric(rownames(clusters))),]

# Merging clusters
clusters$K12 <- clusters$K
clusters$K12[clusters$K12 == 1 | clusters$K12 == 2] <- 12
clusters$K24 <- clusters$K
clusters$K24[clusters$K24 == 2 | clusters$K24 == 4] <- 24
clusters$K56 <- clusters$K
clusters$K56[clusters$K56 == 5 | clusters$K56 == 6] <- 56
write.csv(clusters, "clusters.csv")

# Clustering frequency cutoff mask
par(mfrow = c(3,3), mar=c(rep(1.4,4)))
box_dims <- tapply(1:47, clusters$K[clustering_order], function(x) c(ytop=1-(min(x)-1)/47, ybot=1-max(x)/47, xrig=max(x)/47, xlef=(min(x)-1)/47))
marg_fix <- function(x) -(1-(x/0.5))/(47*6)
for (t in seq(0.1, 0.9, 0.1)) {
  image(t<clustering_frequency[,nrow(clustering_frequency):1], axes=FALSE)
  title(paste("Clustering frequency > ",t))
  for (box in box_dims) {
    rect(xleft = box["xlef"] + marg_fix(box["xlef"]),
    xright = box["xrig"] + marg_fix(box["xrig"]),
    ybottom = box["ybot"] + marg_fix(box["ybot"]),
    ytop = box["ytop"] + marg_fix(box["ytop"]),
    col = NA, border = "black", lty=3)
  }
}
```

### Fst

```Unix
# Calculating Fst
r cmd batch fst_calc.r fst_calc.out &
```

```R
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

C_clusters <- clusters[["K124_56"]]
  cluster_pairs <- t(combn(unique(C_clusters), 2))
  Fsts <- apply(cluster_pairs, 1, function(pair) WC84(geno[C_clusters %in% pair,], C_clusters[C_clusters %in% pair]))
```

```R
library(tidyverse)

# Load the CSV file and split the "pair" column by space delimiter
fsts <- read.csv("fsts.csv", header=1, row.names=1) %>%
  separate(pair, into = c("c1", "c2"), sep = " ") %>%
  mutate(c1 = as.character(c1), c2 = as.character(c2))

# Add reverse pairs
fsts_ <- fsts
fsts_$c1 <- fsts$c2
fsts_$c2 <- fsts$c1
fsts <- rbind(fsts, fsts_) %>% arrange(C, c1, c2)
```

```R
# Shape to tables and save as csv
results <- fsts %>%
  group_by(C) %>%
  nest() %>%
  mutate(data = map(data, ~ {
    .x %>%
      pivot_wider(names_from=c2,values_from=Fst) %>%
      mutate(across(everything(), ~replace(., is.na(.), 0))) %>%
      select(-1) %>%
      select(sort(names(.)))
  })) %>%
  mutate(data = map2(data, C, ~ {
    write.csv(.x, file = paste0("fst/fsts_", .y, ".csv"))
    .x 
  }))
```

```R
# Plotting Fst table as tiles
fsts$Fst_fix <- fsts$Fst # Fix for colouring
fsts$Fst_fix[fsts$Fst_fix <= 0.05] <- 0.05
fsts$Fst_fix[fsts$Fst_fix >= 0.15] <- 0.15
fsts$Fst_fix[is.na(fsts$Fst_fix)] <- 0
fsts$c2 <- factor(fsts$c2, levels = c("6", "5", "4", "3", "2", "1"))


ggplot(data = fsts[fsts$C=="K",], aes(x = c1, y = c2, fill = Fst_fix, label = round(Fst, 3))) +
  geom_tile(color = "black") +
  geom_text(color = "black", size = 3, family = "serif") +
  scale_fill_gradient(low = "white", high = "#7E0605", breaks=c(0,0.05,0.10,0.15)) +
  theme_classic() +
  theme(text = element_text(family = "serif"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
```


## Effect of human colonization

### Inbreeding with PLINK

In UNIX:
```unix
mkdir inbreed
cd inbreed
# Remove missing genotypes
plink --bfile ../AF.imputed.thin --geno 0 --allow-extra-chr --make-bed --out ArcFox0
# Compute heterozygosity per individual
plink --bfile ArcFox0 --allow-extra-chr --het small-sample --out ArcFox0
```
In R:
```R
library(ggplot2)

# Load het table
F<-read.table("ArcFox0.het", h=T)
# Store samples information
samples<-read.table('../sample_popinfo.tsv', h=T)
# Sort samples according to the regions
F<-cbind(F, samples)
F<-F[order(F$Region),]
F<-F[order(F$Country),]

inb_p<-ggplot(F, aes(x=IID, y=F, fill=Region)) +
  geom_bar(stat="identity")+
  theme_classic()+
  theme(plot.title = element_text(size=14,hjust = 0.5), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  scale_x_discrete(limits=F$IID)+xlab('Samples')+ylab('F (Inbreeding coefficient)')+
  ggtitle('Inbreeding coefficient', )+scale_fill_brewer(palette="Dark2")
```
