We started by finding the number of clusters to test for. 
We choose minimum K of 3, and the choice of maximum K was done looking at the number of regions in `sample_popinfo.tsv` with following code: 

```bash
# Finding number of maximum cluster
tail -n +2 sample_popinfo.tsv | cut -f1 | sort | uniq | wc -l
```
Then we wanted to evaluate which K is best for admixture. Therefore, we ran admixture 1 time for each K

```bash
# For admixture evaluation - Running admixture (K in 3-7)

# Assumming K number of ancestral populations
for K in {3..7}
do
   # Run admixture with seed 1
   admixture -s 1 AF.imputed.thin.bed ${K} > AF.imputed_K${K}_run1.log

done

# Show the likelihood of all the runs (in a sorted manner):
grep ^Loglikelihood: *K${K}*log | sort -k2
```
Afterwards, `evalAdmix` was used on each K file.

```bash
# Copying visFuns.R to admixture folder
cp -r visFuns.R /science/groupdirs/jmz230/SCIENCE-BIO-popgen_course-project/Group2_ArcticFox/admixture/

# Running evaluation of admixture output to find optimal K
# Assumed K number of ancestral populations
for K in {3..7}
do
   # Run evalAdmix on admixture output files
   ./evalAdmix -plink AF.imputed.thin -fname AF.imputed.thin.${K}.P -qname AF.imputed.thin.${K}.Q -o K${K}.output.corres.txt

done
```

