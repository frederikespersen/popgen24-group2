We started by finding the number of clusters to test for. 
We choose minimum K of 3, and the choice of maximum K was done looking at the number of regions in `sample_popinfo.tsv` with following code: 

```bash
# Finding number of maximum cluster
tail -n +2 sample_popinfo.tsv | cut -f1 | sort | uniq | wc -l
```

After, we ran admixture 10 times for each K in 3-7.

```bash
# Assumed number of ancestral populations
for K in {3..7}
do
   mkdir admix_K${K}
   
   for i in {1..10}
   do
      # Run admixture with seed i
      admixture -s ${i} AF.imputed.thin.bed ${K} > admix_K${K}/AF.imputed.thin.K${K}_run${i}.log
   
      # Rename the output files
      cp AF.imputed.thin.${K}.Q admix_K${K}/AF.imputed.K${K}_run${i}.Q
      cp AF.imputed.thin.${K}.P admix_K${K}/AF.imputed.K${K}_run${i}.P
      echo "Run ${i} done"
   done
   echo "Runs for ${K} done"
done
```

Afterwards, we found the best run for each K with (Shown for K=3):

```bash
# Show the likelihood of all the 10 runs (in a sorted manner):
cd admix_K3
grep ^Loglikelihood: *K${K}*log | sort -k2
```

Taking the best run for each K, we did evaluation of the model with `evalAdmix`.

```bash
# Running evaluation of admixture output to find optimal K
# Assumed K number of ancestral populations

# K = 3, best run = 1
./evalAdmix -plink AF.imputed.thin -fname AF.imputed.K3_run1.P -qname AF.imputed.K3_run1.Q -o K3.output.corres.txt

# K = 4, best run = 10
./evalAdmix -plink AF.imputed.thin -fname AF.imputed.K4_run10.P -qname AF.imputed.K4_run10.Q -o K4.output.corres.txt

# K = 5, best run = 7
./evalAdmix -plink AF.imputed.thin -fname AF.imputed.K5_run7.P -qname AF.imputed.K5_run7.Q -o K5.output.corres.txt

# K = 6, best run = 8
./evalAdmix -plink AF.imputed.thin -fname AF.imputed.K6_run8.P -qname AF.imputed.K6_run8.Q -o K6.output.corres.txt

# K = 7, best run = 4
./evalAdmix -plink AF.imputed.thin -fname AF.imputed.K7_run4.P -qname AF.imputed.K7_run4.Q -o K7.output.corres.txt


# Copying visFuns.R to admixture folder
cp -r visFuns.R /science/groupdirs/jmz230/SCIENCE-BIO-popgen_course-project/Group2_ArcticFox/admixture/
```

Evaluating the fit of admixture analysis with correlation matrix (done for each `K${K}.output.corres.txt`)

```R
# Read in plotting functions
source("visFuns.R")

# Read in the population info
popinfo <- read.table("sample_popinfo.tsv", header = TRUE)
pop = as.vector(popinfo[,1])

# Read in the output 
r <- as.matrix(read.table("K3.output.corres.txt"))
plotCorRes(cor_mat = r, pop = pop, title = "Correlation of residuals (K=3)", max_z=0.15, min_z=-0.15)
```

Following code was used to make admixture plot of K equal to 3, 4, 5, 6, and 7 taking the best run out 10.

```R
# Margins and colors
par(mar=c(7,3,2,1), mgp=c(2,0.6,0))
palette(c("#E69F00", "#56B4E9", "#D55E00", "#999999", "#66CC00", "#CC0066", "#9999FF"))

# Load sample names
popinfo <- read.table("sample_popinfo.tsv", header = TRUE)
region_names <- popinfo$Region

# Read sample ancestral proportions
snp_k <- as.matrix(read.table("AF.imputed.K7_run4.Q"))

# Ensure row names of snp_k correspond exactly to region_names
rownames(snp_k) <- popinfo$Sample

# Get the order of region_names
order_indices <- order(match(region_names, unique(region_names)))

# Reorder the rows of snp_k
snp_k_sorted <- snp_k[order_indices, ]

# Barplot with sorted data using the given colors
barplot(t(snp_k_sorted), col=c(7,6,5,4,3,2,1), 
        names.arg=region_names[order_indices], cex.names=0.8,
        border=NA, main="K=7 - Run 4", las=2, ylab="Ancestry proportion")
```
