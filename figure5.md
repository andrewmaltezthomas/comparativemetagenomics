# Analysis for figure 5

This tutorial shows how we conducted differential abundance analysis (KEGG modules and pathways) on TD and BM shotgun data from the HMP. 

**Differential abundance analysis of KEGG modules and pathways**

**In the command line**. Filter and join HUMAnN analysis files for KEGG modules and pathways:
```{r eval=FALSE, engine='bash'}
# Create list (a file called dir_name.txt) with subjects identifiers -> which will be used to download subjects HUMAnN analysis results from HMP and naming of folders to store these results files
for i in "SRS013234" "SRS013502" "SRS013506" "SRS013705" "SRS013711" "SRS013825" "SRS016349" "SRS016529" "SRS016533" "SRS016600" "SRS017080" "SRS017127" "SRS017209" "SRS017215" "SRS017533" "SRS017537" "SRS017713" "SRS018145" "SRS018439" "SRS021496"; do 
echo $i >> dir_name.txt;
done

# Download HMP files (for each subject)
cat dir_name.txt | while read name; do 
wget http://downloads.hmpdacc.org/data/HMMRC/"$name"_vs_KEGG_v54.tar.bz2;
done

# Uncompress downloaded files
cat dir_name.txt | while read name; do mkdir "$name"; done
cat dir_name.txt | while read name; do
tar -xf "$name"_vs_KEGG_v54.tar.bz2 -C ./"$name";
done

# Select pathways: a) the pathway must be present in all samples, and b) the coverage must be >= 0.7 for at least one sample
cat SRS*/*_vs_KEGG_v54_04a-mpt-cop-nul-nve-nul-nve-xpe.txt | awk '$2 >= 0.7 {print $1}' | sort | uniq -c | awk -F " " '$1 == 20 {print $2}' | sed "$ d" > all_samples_greater0.7.ko.txt

# For each subject, get the abundance of the selected pathways -> tabular output follows [pathway_ko, pathway_abundance, subject]
cat all_samples_greater0.7.ko.txt | while read p; do
 cat dir_name.txt | while read s; do
  awk -v pat="$p" -v name="$s" '$0 ~ pat {print $0"\t"name}' "$s"/"$s"_vs_KEGG_v54_04b-mpt-cop-nul-nve-nul-nve.txt;
 done
done > Abundance_ko.txt

# Select modules: a) the module must be present in all samples, and b) the coverage must be >= 0.7 for at least one sample
cat SRS*/*_vs_KEGG_v54_04a-mpm-cop-nul-nve-nul-nve-xpe.txt | awk '$2 >= 0.7 {print $1}' | sort | uniq -c | awk -F " " '$1 == 20 {print $2}' | sed "$ d" > all_samples_greater0.7.mo.txt

# For each subject, get the abundance of the selected modules -> tabular output follows [module_mo, pathway_abundance, subject]
cat all_samples_greater0.7.mo.txt | while read m; do
 cat dir_name.txt | while read s; do
  awk -v pat="$m" -v name="$s" '$0 ~ pat {print $0"\t"name}' "$s"/"$s"_vs_KEGG_v54_04b-mpm-cop-nul-nve-nul-nve.txt; 
 done
done > Abundance_mo.txt
```

**In R**. Calculate differential abundance statistics:
```{r eval=FALSE}
library(curl)
library(tidyr)
library(gplots)

## Subjects
buccalmucosa <- c("SRS013506", "SRS013711", "SRS013825", "SRS016349", "SRS016533", "SRS016600", "SRS017080", "SRS017127", "SRS017215", "SRS017537")
tonguedorsum <- c("SRS013234", "SRS013502", "SRS013705", "SRS016529", "SRS017209", "SRS017533", "SRS017713", "SRS018145", "SRS018439", "SRS021496")

## Create folders to store HMP files about shotgun metadata with sequencing depth
system("mkdir -p Read_abundancy/Buccal")
system("mkdir -p Read_abundancy/Tongue")

## Download metadata for BM
for (i in buccalmucosa) { 
  curl_download(url = paste("http://downloads.hmpdacc.org/data/HMSCP/", i, "/", i, "_metric.txt.bz2", sep = ""), 
                destfile = paste("Read_abundancy/Buccal/", i, "_metric.txt.bz2", sep = ""))
}

## Unzip files for BM
system("bzip2 -d Read_abundancy/Buccal/*.bz2")

## Extract sequencing depth for BM
system("cat Read_abundancy/Buccal/*_metric.txt | grep -v 'Sample Name' | cut -f2 > Read_abundancy/buccal_mucosa_reads.txt")

## Download metadata for TD
for (i in tonguedorsum) { 
  curl_download(url = paste("http://downloads.hmpdacc.org/data/HMSCP/", i, "/", i, "_metric.txt.bz2", sep = ""), 
                destfile = paste("Read_abundancy/Tongue/", i, "_metric.txt.bz2", sep = ""))
}

## Unzip files for TD
system("bzip2 -d Read_abundancy/Tongue/*.bz2")

## Extract sequencing depth for TD
system("cat Read_abundancy/Tongue/*_metric.txt | grep -v 'Sample Name' | cut -f2 > Read_abundancy/tongue_dorsum_reads.txt")

## Calculate normalization factor
bm_reads <- read.table("Read_abundancy/buccal_mucosa_reads.txt")
td_reads <- read.table("Read_abundancy/tongue_dorsum_reads.txt")
norm_factor <- mean(bm_reads$V1)/mean(td_reads$V1)

## KEGG Pathways ##

## Load KO abundance data
ko <- read.table("Abundance_ko.txt", header = F, stringsAsFactors = F)
ko <- as.data.frame(ko)

## Transform dataset
ko <- spread(data = ko, key = V3, value = V2)
colnames(ko)[1] <- "KO"

## Reordering columns: KO -> BM Samples -> TD Samples
ko <- ko[,c("KO", buccalmucosa, tonguedorsum)]

## Extract matrices according to sample site
buccal_mucosa_ko <- ko[,which(colnames(ko) %in% buccalmucosa)]
tongue_dorsum_ko <- ko[,which(colnames(ko) %in% tonguedorsum)]

## Normalize
tongue_dorsum_ko[-1,] <- tongue_dorsum_ko[-1,] * norm_factor

## Compute p-values by wilcox test and FDR correction
pvalues <- rep(0, nrow(ko))
for (i in 1:nrow(ko)) {
  pvalues[i] <- wilcox.test(as.numeric(buccal_mucosa_ko[i,-1]), as.numeric(tongue_dorsum_ko[i,-1]))$p.value
}
pvalues <- p.adjust(pvalues, "fdr")
ko$Pvalue <- pvalues

## KEGG Modules ##

## Load MO abundance data
mo <- read.table("Abundance_mo.txt", header = F, stringsAsFactors = F)
mo <- as.data.frame(mo)

## Transform dataset
mo <- spread(data = mo, key = V3, value = V2)
colnames(mo)[1] <- "MO"

## Reordering columns: MO -> BM Samples -> TD Samples
mo <- mo[,c("MO", buccalmucosa, tonguedorsum)]

## Extract matrices according to sample site
buccal_mucosa_mo <- mo[,which(colnames(mo) %in% buccalmucosa)]
tongue_dorsum_mo <- mo[,which(colnames(mo) %in% tonguedorsum)]

## Normalize
tongue_dorsum_mo <- tongue_dorsum_mo * norm_factor

## Compute p-values by wilcox test and FDR correction
pvalues <- rep(0, nrow(mo))
for (i in 1:nrow(mo)) {
  pvalues[i] <- wilcox.test(as.numeric(buccal_mucosa_mo[i,-1]), as.numeric(tongue_dorsum_mo[i,-1]))$p.value
}
pvalues <- p.adjust(pvalues, "fdr")
mo$Pvalue <- pvalues

## Select significant modules by p-values
significant_modules <- mo[which(mo$Pvalue < 0.05),]

## Get modules descriptions by KEGG REST API
for (i in 1:nrow(significant_modules)) {
  mod_desc_url <- paste("http://rest.kegg.jp/list/md:", significant_modules$MO[i], sep = "")
  mod_desc <- read.table(mod_desc_url, sep = "\t")
  rownames(significant_modules)[i] <- as.character(unlist(mod_desc[2]))
}

## Edit modules descriptions
for (i in grep("=>", rownames(significant_modules))) {
  rownames(significant_modules)[i] <- 
    unlist(strsplit(rownames(significant_modules)[i], split = " =>"))[1]
}
for (i in which(nchar(rownames(significant_modules)) > 47)) {
  rownames(significant_modules)[i] <- gsub(", [A-Za-z0-9 ]+$", "", rownames(significant_modules)[i])
}

# Plot Heatmap
heatmap.2(t(log10(significant_modules[,c(-1,-22)])),
          RowSideColors = c(replicate(10, "#67a9cf"), replicate(10, "#fc8d59")),
          cexRow = 0.7, cexCol = 0.8,
          dendrogram = "row",
          margins = c(16, 8), lhei = c(3, 8),
          trace = "none", density.info = "none")
```
