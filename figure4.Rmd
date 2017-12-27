# Analysis for figure 4

This tutorial shows how we conducted differential abundance analysis (taxons) on TD and BM shotgun data from the HMP. 

**Differential abundance analysis of taxa**

**In R**. Get and prepare MetaPhlAn data:
```{r eval=FALSE}
## Subjects
buccalmucosa <- c("SRS013506", "SRS013711", "SRS013825", "SRS016349", "SRS016533", "SRS016600", "SRS017080", "SRS017127", "SRS017215", "SRS017537")
tonguedorsum <- c("SRS013234", "SRS013502", "SRS013705", "SRS016529", "SRS017209", "SRS017533", "SRS017713", "SRS018145", "SRS018439", "SRS021496")

## Get MetaPhlAn data and Prepare
metaphlan <- read.table("HMP.ab.txt", header = T, row.names = 1, stringsAsFactors = F)
metaphlan <- metaphlan[,c(buccalmucosa, tonguedorsum)]
metaphlan <- metaphlan[-c(which(rowSums(metaphlan) == 0)),]
genera.i <- grep("^(\\w+\\|){5}[A-Za-z\\_]+$", rownames(metaphlan))
metaphlan <- metaphlan[genera.i,]
write.table(metaphlan, "Metaphlan_filtered_genus.tsv", quote = F, sep = "\t", row.names = T, col.names = T)
```

**In the command line**. Linear discriminant analysis (LDA) using LEfSe on MetaPhlAn taxon table:
```{r eval=FALSE, engine='bash'}
# Tabular data for class information
FirstLine="bodysite \t BuccalMucosa \t BuccalMucosa \t BuccalMucosa \t BuccalMucosa \t BuccalMucosa \t BuccalMucosa \t BuccalMucosa \t BuccalMucosa \t BuccalMucosa \t BuccalMucosa \t TongueDorsum \t TongueDorsum \t TongueDorsum \t TongueDorsum \t TongueDorsum \t TongueDorsum \t TongueDorsum \t TongueDorsum \t TongueDorsum \t TongueDorsum"

# Add classes and format subject descriptions
sed -i 's/SRS013506/'"$FirstLine"'\nid\tSRS013506/g' Metaphlan_filtered_genus.tsv
sed -i 's/SRS//g' Metaphlan_filtered_genus.tsv

# Format the table (making sure the class information is in the first row, and subject information in the second)
lefse/format_input.py Metaphlan_filtered_genus.tsv merged_abundance_table.lefse -c 1 -u 2 -o 1000000

# Run Lefse with relaxed parameters to retrieve LDAs and P values for most taxa (will be used to plot all genera)
lefse/run_lefse.py -l 0.1 -a 0.99 -w 0.99 merged_abundance_table.lefse merged_abundance_table_relaxed.lefse.out

# Run Lefse with strict parameters to retrieve LDAs and P values for only most significant taxa (will be used to plot only biomarker genera)
lefse/run_lefse.py -l 2.0 -a 0.02 -w 0.02 merged_abundance_table.lefse merged_abundance_table_significant.lefse.out

# Plot cladogram with only most significant taxa (below command is for those who desire to plot in command line, we used the galaxy version)
lefse/plot_cladogram.py --dpi 300 --format png merged_abundance_table_significant.lefse.out lefse_biomarkers_cladogram.png
lefse/plot_cladogram.py --dpi 300 --format pdf merged_abundance_table_significant.lefse.out lefse_biomarkers_cladogram.pdf
```

**In R**. Plot the figure:
```{r eval=FALSE}
## Colors auxiliary variables
BuccalMucosa_Col = "#67a9cf"
TongueDorsum_Col = "#fc8d59"

## Load LEfSe result data
lda.data <- read.delim("merged_abundance_table_relaxed.lefse.out", header = F, col.names = c("Genera", "value", "Group", "lda", "pvalue"))

## Remove non marker taxons
lda.data <- subset(lda.data, Group != "")

## Keep only genera to plot
lda.genera.i <- grep("^(\\w+\\.){5}[A-Za-z\\_]+$", lda.data$Genera)
lda.data <- lda.data[lda.genera.i,]

## Create a colors list for each point according to groups
lda.colors <- replicate(nrow(lda.data), BuccalMucosa_Col)
lda.colors[lda.data$Group != "BuccalMucosa"] <- TongueDorsum_Col

## Separate Buccal Mucosa and Tongue Dorsum points in the graph
lda.data$lda[lda.data$Group == "BuccalMucosa"] <- lda.data$lda[lda.data$Group == "BuccalMucosa"] * (-1)

## Make sure pvalue is numeric for R
lda.data$pvalue <- as.numeric(as.character(lda.data$pvalue))

## Manual observation of data to be highlighted in the graph with Genera names (significant genera names and max LDA scores)
subset(lda.data, pvalue < 0.05 & Group == "BuccalMucosa")[,c("Genera", "lda", "pvalue")]
subset(lda.data, pvalue < 0.05 & Group == "TongueDorsum")[,c("Genera", "lda", "pvalue")]
range(lda.data$lda)

## Plot the data
par(mfrow=c(1,1))
plot(lda.data$pvalue ~ lda.data$lda, main = "LEfSe Analysis", xlab = "LDA Score", ylab = "", 
     axes = F, bty = 'n',
     col = lda.colors, pch = 16, log ="y", xlim=c(-6,7), ylim = rev(range(lda.data$pvalue,na.rm = T)))

## Configure axes: x -> LDA score, y -> p-value
axis(1, -6:7, pos = 1.01, cex.axis = 0.8)

axis(2, at = c(0.0001,10^(-0),10^(-1),10^(-2),10^(-3)),
     labels = c("", "", expression(10^-1), expression(10^-2), expression(10^-3)),
     pos = 0, cex.axis = 0.8, las = 2)

## Sign the P-value axis
text(-0.7, 0.00031, "P-value", cex = 0.7)

## Mark significant P threshold in the graph
abline(h = 0.05, lty = 2)
text(-5.5, 0.04, "P = 0.05", cex = 0.7)

## Add genera names to highlight their LDA scores points
text(-5.260745, 0.0005, expression(italic("Streptococcus")), cex = 0.6)
text(-4.448734, 0.0008, expression(italic("Gemella")), cex = 0.6)
text(-2.958624, 0.01, expression(italic("Propionibacterium")), cex = 0.6)
text(5.11203, 0.006501702, expression(italic("Prevotella")), cex = 0.6)
text(4.608213, 0.0003, expression(italic("Veillonella")), cex = 0.6)

## Add arrows to sign the Buccal Mucosa and Tongue Dorsum groups
mtext("Buccal Mucosa enriched", 1, at = -4, padj = 4, cex = 0.8)
arrows(-2, 3.5, -6.5, 3.5, xpd = T, cex = 0.8, length = 0.1)

mtext("Tongue Dorsum enriched", 1, at = 5, padj = 4, cex = 0.8)
arrows(3, 3.5, 7.5, 3.5, xpd = T, cex = 0.8, length = 0.1)
```
