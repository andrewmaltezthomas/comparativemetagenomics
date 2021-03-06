# Analysis for figure 1

**All of the code in this page is meant to be run on the command line unless otherwise stated.**
These are the commands we used to generate sample/sequence rarefaction curves and calculate different alpha diversity metrics using qiime:
```{r eval=FALSE, engine='bash'}
# Convert table to biom format
biom convert -i OTU-Table-16S.tsv -o OTU-Table-16S.biom --table-type="OTU table" --to-json --process-obs-metadata taxonomy

# Download phylogenetic tree 
wget http://downloads.hmpdacc.org/data/HMQCP/rep_set_v35.tre.gz

# Unzip it
gzip -d rep_set_v35.tre.gz

# Get OTU ids present in the analyzed samples only
tail -n +2 OTU-Table-16S.tsv | cut -f1 > OTUs_to_keep.txt

# Prune the phylogenetic tree 
filter_tree.py -i rep_set_v35.tre -t OTUs_to_keep.txt -o pruned_rep_set_v35.tree

# Filter low abundant otus
filter_otus_from_otu_table.py -i OTU-Table-16S.biom -n 2 -s 3 -o OTU-Table-16S.filtered.biom

# Get the stats for the filtered table
biom summarize-table -i OTU-Table-16S.filtered.biom -o OTU-Table-16S.filtered.summary

# Rarefy OTU table
multiple_rarefactions_even_depth.py -i OTU-Table-16S.filtered.biom -d 1802 -o Rarefied_Tables/

# Calculate alpha diversity on each of the rarefied tables
alpha_diversity.py -i Rarefied_Tables/ -m observed_species,PD_whole_tree,shannon,equitability -t pruned_rep_set_v35.tree -o Alpha_diversity/

# Collate different alpha diversity tables
collate_alpha.py -i Alpha_diversity/ -o Alpha_Diversity_Collated/

# Parametric testing of alpha diversity
for i in Alpha_Diversity_Collated/*; 
do
 compare_alpha_diversity.py -i $i -m Qiime.map -c Description -o $i.statistics; 
done

for i in Alpha_Diversity_Collated/*.statistics; 
do
 measure=$(echo $i | cut -f2 -d"/" | cut -f1 -d"."); 
 echo $measure; grep -A 1 "p-value" $i/Description_stats.txt; 
done

# Make parameters file
echo "alpha_diversity:metrics observed_species" > alpha_params.txt 

# Rarefaction
alpha_rarefaction.py -i OTU-Table-16S.filtered.biom -m Qiime.map -e 1802 -p alpha_params.txt -o Alpha_Rarefaction
```

**In R**. Pre-process and plot figure 1A:
```{r eval=FALSE}
library(vegan)

## Subjects
buccalmucosa <- c("SRS013506", "SRS013711", "SRS013825", "SRS016349", "SRS016533", "SRS016600", "SRS017080", "SRS017127", "SRS017215", "SRS017537")
tonguedorsum <- c("SRS013234", "SRS013502", "SRS013705", "SRS016529", "SRS017209", "SRS017533", "SRS017713", "SRS018145", "SRS018439", "SRS021496")

## Colors auxiliary variables
BuccalMucosa_Col = "#67a9cf"
TongueDorsum_Col = "#fc8d59"

## Load OTUs abunance matrix data
OTUsAbundance <- read.delim("OTU-Table-16S.tsv", TRUE, "\t")

## Rename rows (with OTUs identifiers)
rownames(OTUsAbundance) <- OTUsAbundance$OTU_ID

## Remove first (OTUs identifiers) and last (taxonomy) columns
toRemove <- c(-1, -(ncol(OTUsAbundance)))
OTUsAbundance <- OTUsAbundance[,toRemove]

## Compute species accumulation curve (SAC) for Buccal Mucosa
SpeciesCollect1 <- specaccum(t(OTUsAbundance[,buccalmucosa]), method = "rarefaction")

## Compute species accumulation (SAC) curve for Tongue Dorsum
SpeciesCollect2 <- specaccum(t(OTUsAbundance[,tonguedorsum]), method = "rarefaction")

## Set plot area parameters
par(mar = c(2.3, 2.45, 1.9, 0.75), oma = c(0, 0, 0, 0), mgp = c(0, 0.2, 0))

## Plot SAC for Buccal Mucosa
plot(SpeciesCollect1, ci.type = "bar", col = BuccalMucosa_Col, ci.col = BuccalMucosa_Col, lwd = 2,
     main = "Sampling Effort", xlab = "", ylab = "",
     cex.main = 0.7, cex.axis = 0.6, cex.lab = 0.6, tck = -0.02)

## Add SAC for Tongue Dorsum
plot(SpeciesCollect2, ci.type = "bar", col = TongueDorsum_Col, ci.col = TongueDorsum_Col, lwd = 2,
     cex.axis = 0.6, cex.lab = 0.6, tck = -0.02,
     add = T)

## Add axes description
title(xlab = "Samples", ylab = "Number of OTUs", line = 1.3, cex.lab = 0.6)
title(main = "A", outer = T, adj = 0.01, line = -1.6, cex.main = 1.2)

## Add legends
legend("bottomright", c("Buccal Mucosa","Tongue Dorsum"),
       lty = c(1, 1), lwd = c(2.5, 2.5), col = c(BuccalMucosa_Col, TongueDorsum_Col), cex = 0.5)
```

**In R**. Pre-process and plot figure 1B:
```{r eval=FALSE}
library(reshape2)

## Subjects
buccalmucosa <- c("SRS013506", "SRS013711", "SRS013825", "SRS016349", "SRS016533", "SRS016600", "SRS017080", "SRS017127", "SRS017215", "SRS017537")
tonguedorsum <- c("SRS013234", "SRS013502", "SRS013705", "SRS016529", "SRS017209", "SRS017533", "SRS017713", "SRS018145", "SRS018439", "SRS021496")

## Colors auxiliary variables
BuccalMucosa_Col = "#67a9cf"
TongueDorsum_Col = "#fc8d59"

## Auxiliary function to compute standard error (SE)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

## Auxiliary functionto add errors bars to the plot
addErrorsBars <- function(plotObject, rowsIndexes, plotColor) {
  for (i in rowsIndexes) {
    x = plotObject$variable[i]
    y = plotObject$value[i]
    sd = plotObject$se[i]
    
    ## Adpated from http://stackoverflow.com/questions/15063287/add-error-bars-to-show-standard-deviation-on-a-plot-in-r
    segments(x, y-sd,x, y+sd, col = plotColor, lwd = 2)
    epsilon = 0.02
    segments(x-epsilon,y-sd,x+epsilon,y-sd, col = plotColor, lwd = 2)
    segments(x-epsilon,y+sd,x+epsilon,y+sd, col = plotColor, lwd = 2)
  }
}

## Load data
rarefaction <- read.delim("Alpha_Rarefaction/alpha_div_collated/observed_species.txt", header = T, row.names = 1)

## Establish sampling intervals
sequences <- seq(10, 110, 20)

## Matrix to store richness for each interval
observed_otus <- matrix(0, nrow = 6, ncol = 21)

## Rarefaction method applying
for (i in 1:length(sequences)) {
  start <- sequences[i] - 9
  observed_otus[i,1] <- rarefaction[sequences[i],1] 
  for (j in 2:21) {
    observed_otus[i,j] <- mean(rarefaction[start:sequences[i],j])
  }
}

## Transpose the data -> put the samples in the rows
observed_otus <- t(observed_otus)

## Set rows names with samples identifiers
rownames(observed_otus) <- colnames(rarefaction)[-2]

## Store iteration (first row) as columns names and remove it from data
colnames(observed_otus) <- observed_otus[1,]
observed_otus <- observed_otus[-1,]

## Convert matrix to data frame
observed_otus <- as.data.frame(observed_otus)

## Add the Group (Buccal Mucosa or Tongue Dorsum) column
observed_otus$Group <- ""
observed_otus[buccalmucosa,]$Group <- "Buccal Mucosa"
observed_otus[tonguedorsum,]$Group <- "Tongue Dorsum"

## Add samples identifiers as column
samples <- rownames(observed_otus)
observed_otus <- as.data.frame(cbind(samples, observed_otus))

## Convert data structure to a long format
observed_otus <- melt(observed_otus, id.vars = c("samples", "Group"))
# observed_otus <- observed_otus[order(observed_otus$variable),]

## Convert variable (which represents iteration) to numeric
observed_otus$variable <- as.numeric(as.character(observed_otus$variable))

## Compute standard erros and prepare data to plot
plotData <- summarySE(observed_otus, measurevar="value", groupvars=c("Group", "variable"))

## Set plot area parameters
par(mar=c(2.3,2.45,1.9,0.75), oma=c(0,0,0,0), mgp=c(0,0.2,0))

## Plot the Buccal Mucosa data
plot(subset(plotData, Group == "Buccal Mucosa")$variable,
     subset(plotData, Group == "Buccal Mucosa")$value,
     type = "l", lwd = 2, col = BuccalMucosa_Col,
     ylim = c(0, 500),
     main = "Sequencing Effort", xlab = "", ylab = "", cex.main = 0.7, cex.lab = 0.6,
     axes = F)

## Plot the Tongue Dorsum data
lines(subset(plotData, Group == "Tongue Dorsum")$variable,
      subset(plotData, Group == "Tongue Dorsum")$value,
      type = "l", lwd = 2, col = TongueDorsum_Col)

## Plot the errors bars
addErrorsBars(plotData, 1:6, BuccalMucosa_Col)
addErrorsBars(plotData, 7:12, TongueDorsum_Col)

## Add axes and the margin lines
axis(side = 1, at = plotData$variable[1:6], cex.axis = 0.6, tck = -0.02)
axis(side = 2, at = c(10, 100, 200, 300, 400, 500), cex.axis = 0.6,tck = -0.02)
box()

## Add axes labels
title(xlab = "Sequences", ylab = "Number of OTUs", line = 1.3, cex.lab = 0.6)

## Add legends
legend("bottomright", c("Buccal Mucosa","Tongue Dorsum"),
       lty = c(1, 1), lwd = c(2.5, 2.5), col = c(BuccalMucosa_Col, TongueDorsum_Col), cex = 0.5)

```

**In R**. Pre-process and plot figure 1C:
```{r eval=FALSE}
library(reshape2)

## Subjects
buccalmucosa <- c("SRS013506", "SRS013711", "SRS013825", "SRS016349", "SRS016533", "SRS016600", "SRS017080", "SRS017127", "SRS017215", "SRS017537")
tonguedorsum <- c("SRS013234", "SRS013502", "SRS013705", "SRS016529", "SRS017209", "SRS017533", "SRS017713", "SRS018145", "SRS018439", "SRS021496")

## Colors auxiliary variables
BuccalMucosa_Col = "#67a9cf"
TongueDorsum_Col = "#fc8d59"

## Load alpha diversity data
Equitability <- read.delim("Alpha_Diversity_Collated/equitability.txt")
ObservedSpecies <- read.delim("Alpha_Diversity_Collated/observed_species.txt")
PDWholeTree <- read.delim("Alpha_Diversity_Collated/PD_whole_tree.txt")
Shannon <- read.delim("Alpha_Diversity_Collated/shannon.txt")

## Remove metadata
Equitability[,1:3] <- NULL
ObservedSpecies[,1:3] <- NULL
PDWholeTree[,1:3] <- NULL
Shannon[,1:3] <- NULL

## Compute means -> Convert to long format -> Add sample group information (Buccal Mucosa or Tongue Dorsum)
## For all diversity measures

EquitabilityMeans <- sapply(Equitability, mean)
EquitabilityMeans_Melted <- melt(EquitabilityMeans)
EquitabilityMeans_Melted$group <- ""
EquitabilityMeans_Melted[buccalmucosa,]$group <- "Buccal Mucosa"
EquitabilityMeans_Melted[tonguedorsum,]$group <- "Tongue Dorsum"

SpeciesRichness <- sapply(ObservedSpecies, mean)
SpeciesRichness_Melted <- melt(SpeciesRichness)
SpeciesRichness_Melted$group <- ""
SpeciesRichness_Melted[buccalmucosa,]$group <- "Buccal Mucosa"
SpeciesRichness_Melted[tonguedorsum,]$group <- "Tongue Dorsum"

PhylogeneticDiversity <- sapply(PDWholeTree, mean)
PhylogeneticDiversity_Melted <- melt(PhylogeneticDiversity)
PhylogeneticDiversity_Melted$group <- ""
PhylogeneticDiversity_Melted[buccalmucosa,]$group <- "Buccal Mucosa"
PhylogeneticDiversity_Melted[tonguedorsum,]$group <- "Tongue Dorsum"

ShannonDiversity <- sapply(Shannon, mean)
ShannonDiversity_Melted <- melt(ShannonDiversity)
ShannonDiversity_Melted$group <- ""
ShannonDiversity_Melted[buccalmucosa,]$group <- "Buccal Mucosa"
ShannonDiversity_Melted[tonguedorsum,]$group <- "Tongue Dorsum"

## Set plot area parameters
par(mfrow = c(1, 4), mar = c(2.3, 2.45, 1.9, 0.75), oma = c(0, 0, 0, 0), mgp = c(0, 0.2, 0))

## Plot equitability data
boxplot(value ~ group, data = EquitabilityMeans_Melted,
        main = "Equitability", xlab = "", ylab = "", names = c("BM", "TD"),
        cex.main = 0.6, cex.axis = 0.6, cex.lab = 0.6,
        yaxt = "n", vertical = T)

## Add points
stripchart(value ~ group, data = EquitabilityMeans_Melted,
           pch = 21, bg = c("black"), col = c(BuccalMucosa_Col, TongueDorsum_Col), cex = 1,
           vertical = T, add = T)

## Add y axis values
axis(side = 2, at = seq(from = 0.70, to = 0.85, by = .05), tck = -0.02, cex.axis = 0.6)

## Add p-value
text(2.15, 0.68, "P = 0.05", cex = 0.55)

## Add y axis title
title(ylab = "Alpha Diversity Measure", line = 1.3, cex.lab = 0.6)

## Plot richness data
boxplot(value ~ group, data = SpeciesRichness_Melted,
        main = "Richness", xlab = "", ylab = "", names = c("BM", "TD"),
        cex.main = 0.6, cex.axis = 0.6,
        yaxt = "n", vertical = T)

## Add points
stripchart(value ~ group, data = SpeciesRichness_Melted,
           pch = 21, bg = c("black"), col = c(BuccalMucosa_Col, TongueDorsum_Col), cex = 1,
           vertical = T, add = T)

## Add y axis values
axis(side = 2, at = seq(from = 250, to = 500, by = 50), tck = -0.02, cex.axis = 0.6)

## Add p-value
text(2.15, 235, "P = 0.09", cex = 0.55)

## Plot phylogenetic diversity data
boxplot(value ~ group, data = PhylogeneticDiversity_Melted,
        main = "Phylogenetic Diversity", xlab = "", ylab = "", names = c("BM", "TD"),
        cex.main = 0.6, cex.axis = 0.6,
        yaxt = "n", vertical = T)

## Add points
stripchart(value ~ group, data = PhylogeneticDiversity_Melted,
           pch = 21, bg = c("black"), col = c(BuccalMucosa_Col, TongueDorsum_Col), cex = 1,
           vertical = T, add = T)

## Add y axis values
axis(side = 2, at = seq(from = 8, to = 18, by = 2), tck = -0.02, cex.axis = 0.6)

## Add p-value
text(2.15, 8.3, "P = 0.8", cex = 0.55)

## Plot shannon diversity data
boxplot(value ~ group, data = ShannonDiversity_Melted,
        main = "Shannon Diversity", xlab = "", ylab = "", names = c("BM", "TD"),
        cex.main = 0.6, cex.axis = 0.6,
        yaxt = "n", vertical = T)

## Add points
stripchart(value ~ group, data = ShannonDiversity_Melted,
           pch = 21, bg = c("black"), col = c(BuccalMucosa_Col, TongueDorsum_Col), cex = 1,
           vertical = T, add = T)

## Add y axis values
axis(side = 2, at = seq(from = 5.5, to = 8, by = .5), tck = -0.02, cex.axis = 0.6)

## Add p-value
text(2.15, 5.6, "P = 0.05", cex = 0.55)
```
