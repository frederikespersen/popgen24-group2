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

- Filter some 0.15 maf
     - Does it make a difference?
- Nice plot
    - Sort legend by location
    - Color by country
- Map
    - Matching a color

- Try to subdivide populations by PCAs
    - Kmeans clustering? 3,4,5 k
    - Fst matrix

### PCA
```R
# Number of missing samples per site
nMis <- colSums(is.na(geno))

# Only keep sites with 0 missing samples.
geno <- geno[,nMis==0]

# Perform PCA
pca <- prcomp(geno, scale=T, center=T)

# Show summary
summary(pca)

# Extract importance of PCs.
pca_importance <- summary(pca)$importance

# Extract percentage of the variance that is explained by PC1-3
PC1_explained <- round(pca_importance[2,1]*100, 1)
PC2_explained <- round(pca_importance[2,2]*100, 1)
PC3_explained <- round(pca_importance[2,3]*100, 1)

# Extract the PCs
pcs <- as.data.frame(pca$x)

# Plot
palette(c("#7E0605", "#C00E38", "#FF0000", "#FF7D7B", "#B82EAA", "#13501B", "#17A238"))
par(mfrow=c(2,2))
plot(pca_importance[2,], type='b', xlab='PC', ylab='Proportion of variance', las=1,
	pch=19, col='darkred', bty='L', main='Proportion of variance explained per PC')
plot(pcs$PC1, pcs$PC2, col=geo$Region, pch=19, las=1, bty='L',
     main='PCA on 47 arctic foxes',
     xlab=paste0('PC1 (', PC1_explained, '% of variance)'),
     ylab=paste0('PC2 (', PC2_explained, '% of variance)'))
legend('topright', legend=levels(geo$Region), col=1:length(levels(geo$Region)), pch=19)
plot(pcs$PC1, pcs$PC3, col=geo$Region, pch=19, las=1, bty='L',
     main='PCA on 47 arctic foxes',
     xlab=paste0('PC1 (', PC1_explained, '% of variance)'),
     ylab=paste0('PC3 (', PC3_explained, '% of variance)'))
plot(pcs$PC2, pcs$PC3, col=geo$Region, pch=19, las=1, bty='L',
     main='PCA on 47 arctic foxes',
     xlab=paste0('PC2 (', PC2_explained, '% of variance)'),
     ylab=paste0('PC3 (', PC3_explained, '% of variance)'))
```