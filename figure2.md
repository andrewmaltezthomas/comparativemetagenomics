# Analysis for figure 2

**All of the code in this page is meant to be run in R unless otherwise stated.**
These are the commands we used to calculate and generate boxplots for gene and genera richness. We would like to point out that it's not possible to reproduce the steps used to calculate the total number of genera analyzed via MG-RAST because of changes to their online platform, whereby [hierarchical .tsv abundance files were no longer being generated.](https://groups.google.com/d/topic/mg-rast/ft-jFL3MTrk/discussion)
```{r eval=FALSE}
library("curl")
library("gdata")
library("reshape2")

## Subjects
buccalmucosa <- c("SRS013506", "SRS013711", "SRS013825", "SRS016349", "SRS016533", "SRS016600", "SRS017080", "SRS017127", "SRS017215", "SRS017537")
tonguedorsum <- c("SRS013234", "SRS013502", "SRS013705", "SRS016529", "SRS017209", "SRS017533", "SRS017713", "SRS018145", "SRS018439", "SRS021496")

## Colors auxiliary variables
BuccalMucosa_Col = "#67a9cf"
TongueDorsum_Col = "#fc8d59"

###------------------------------------------ ## MetaPhlAn ### ---------------------------------------

## Get MetaPhlAn data and Prepare
metaphlan <- read.table("HMP.ab.txt", header = T, row.names = 1, stringsAsFactors = F)
metaphlan <- metaphlan[,c(buccalmucosa, tonguedorsum)]
metaphlan <- metaphlan[-c(which(rowSums(metaphlan) == 0)),]
genera.i <- grep("^(\\w+\\|){5}[A-Za-z\\_]+$", rownames(metaphlan))
metaphlan <- metaphlan[genera.i,]
write.table(metaphlan, "Metaphlan_filtered_genus.tsv", quote = F, sep = "\t", row.names = T, col.names = T)

## Load MetaPhlAn filtered (genus only) data
table <- read.table("Metaphlan_filtered_genus.tsv", header = T, stringsAsFactors = F)

## Create array to store genera counts
total_genera_meta <- table[1,]
rownames(total_genera_meta) <- "counts"

## Compute and store genera counts
samples <- colnames(table)
for (sample in samples) {
  number_genera <- length(which(table[,sample] > 0))
  total_genera_meta[1,sample] <- number_genera
}

## Asses means differences thorugh t-test
p.metaphlan <- wilcox.test(as.numeric(total_genera_meta[,buccalmucosa]), as.numeric(total_genera_meta[,tonguedorsum]))
p.metaphlan <- t.test(as.numeric(total_genera_meta[,buccalmucosa]), as.numeric(total_genera_meta[,tonguedorsum]))

## Prepare text to plot the test p-value
p.metaphlan.text <- paste("P = ", round(p.metaphlan$p.value, 3), collapse = "",  sep = "")

##----------------------------------------- ## Gene counts ### ---------------------------------------

## Download sample summary data (with gene counts)
curl_download("http://downloads.hmpdacc.org/data/HMGI/sample_summary.xls", destfile = "sample_summary.xls")

## Load sample summary data
xls <- read.xls("sample_summary.xls", sheet = 1)

## Subset gene counts from sample summary data
gene.count <- subset(xls, SAMPLE_ID %in% samples)[,1:2]

## Correct GENE_COUNT column to numeric (remove ",")
gene.count$GENE_COUNT <- as.numeric(gsub(",", "", gene.count$GENE_COUNT))

## Asses means differences thorugh t-test
rownames(gene.count) <- gene.count$SAMPLE_ID
p.functional <- wilcox.test(gene.count[buccalmucosa,]$GENE_COUNT, gene.count[tonguedorsum,]$GENE_COUNT)

## Prepare text to plot the test p-value
p.functional.text <- "P <0.001"

## Prepare Genera data to plot
groups <- c(rep("Bucccal.Mucosa", 10), rep("Tongue.Dorsum", 10))
alpha_diversity <- rbind(total_genera_meta, groups)
alpha_diversity <- as.data.frame(t(total_genera_meta))
alpha_diversity$Group <- ""
alpha_diversity[buccalmucosa,]$Group <- "Buccal.Mucosa"
alpha_diversity[tonguedorsum,]$Group <- "Tongue.Dorsum"

## Prepare Genes data to plot
functional_diversity <- gene.count
functional_diversity$Group <- ""
functional_diversity[buccalmucosa,]$Group <- "Buccal.Mucosa"
functional_diversity[tonguedorsum,]$Group <- "Tongue.Dorsum"

# Plot the genera data 
boxplot(counts ~ Group, data = alpha_diversity ,
        main = "MetaPhlAn", xlab = "Sample Type", ylab = "Number of Genera", names = c("BM", "TD"), vertical = T,
        yaxt = "n", xaxt = "n")
stripchart(counts ~ Group, data = alpha_diversity, pch = 21, bg = c("black"), col = c(BuccalMucosa_Col, TongueDorsum_Col), cex = 1.1,
           vertical = T, add = T)

## Configure axes: x -> Sample Type, y -> Number of Genera
axis(side = 1, at = c(1,2), labels = c("BM", "TD"), mgp = c(-0.2, 0.1, 0), tck = -0.02, cex.axis = 0.8)
axis(side = 2, at = seq(from = 25, to = 55, by = 5), mgp = c(-0.2, 0.1, 0), tck = -0.02, cex.axis = 0.8)

## Add p-value
text(2, min(alpha_diversity$counts), p.metaphlan.text, cex = 0.7)

## Plot the genes data
boxplot(GENE_COUNT ~ Group, data = functional_diversity, 
        main = "Functional\nRichness", xlab = "Sample Type", ylab = "Number of Genes", names = c("BM", "TD"), vertical = T,
        yaxt = "n", xaxt = "n")
stripchart(GENE_COUNT ~ Group, data = functional_diversity, pch = 21, bg = c("black"), col = c(BuccalMucosa_Col, TongueDorsum_Col), cex = 1.1,
           vertical = T, add = T)

## Configure axes: x -> Sample Type, y -> Number of Genes
axis(side = 1, at = c(1,2), labels = c("BM", "TD"), mgp=c(-0.2, 0.1, 0), tck=-0.02, cex.axis = 0.8)
axis(side = 2, at = seq(from = 50000, to = 250000, by = 50000), mgp=c(-0.2, 0.1, 0), tck=-0.02, cex.axis = 0.8) 

## Add p-value
text(2, min(functional_diversity$GENE_COUNT), p.functional.text, cex = 0.7)
```

