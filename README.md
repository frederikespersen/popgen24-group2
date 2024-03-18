# Population Genetics Exam 2024 - Group 2

##### Quick commands
```unix
# Navigate to group working directory
cd ~/groupdirs/SCIENCE-BIO-popgen_course-project/Group2
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
cd ~/groupdirs/SCIENCE-BIO-popgen_course-project/Group2
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

# Loading geographic data
geo <- read.table("sample_popinfo.tsv",header=1)
```

```R
# Quality check

```
---

## Population substructure

### PCR
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
plot(pca_importance[2,], type='b', xlab='PC', ylab='Proportion of variance', las=1,
	pch=19, col='darkred', bty='L', main='Proportion of variance explained per PC')

# Extract percentage of the variance that is explained by PC1 and PC2
PC1_explained <- round(pca_importance[2,1]*100, 1)
PC2_explained <- round(pca_importance[2,2]*100, 1)

# Extract the PCs
pcs <- as.data.frame(pca$x)

# Plot
plot(pcs$PC1, pcs$PC2, col=geo$Region, pch=19, las=1, bty='L',
     main='PCA on 29 wild-born chimpanzees',
     xlab=paste0('PC1 (', PC1_explained, '% of variance)'),
     ylab=paste0('PC2 (', PC2_explained, '% of variance)'))
legend('topright', legend=levels(geo$Region), col=1:length(levels(geo$Region)), pch=19)
```