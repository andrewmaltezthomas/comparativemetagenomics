# Analysis for figure 3

**All of the code in this page is meant to be run in R unless otherwise stated.**

These are the commands we used to generate beta diversity dissimilarity boxplots and PCoA plots.

**Commmand line**. Calculate beta diversity and inter/intra group distances using qiime :
```{r eval=FALSE, engine='bash'}
# Calculate beta diversity using QIIME
beta_diversity.py -i Rarefied_Tables/rarefaction_1802_0.biom -m unweighted_unifrac,weighted_normalized_unifrac,bray_curtis -t pruned_rep_set_v35.tree -o Beta_Diversity/

# Calculate intra/inter group distances using QIIME
for i in Beta_Diversity/*.txt; 
do
 distance=$(echo $i | cut -f2 -d"/" | cut -f1,2 -d"_"); 
 make_distance_boxplots.py -d $i -m Qiime.map -f Description --save_raw_data -o Distances/$distance; 
done

# Convert rarefied table to JSON biom format (to be used in R)
biom convert -i Rarefied_Tables/rarefaction_1802_0.biom --to-json -o Rarefied_Tables/rarefaction_1802_0.json.biom
```

Plot the results in R:
```{r eval=FALSE}
library("reshape2")
library("phyloseq")
library("ggplot2")
library("plyr")

## Colors auxiliary variables
BuccalMucosa_Col = "#67a9cf"
TongueDorsum_Col = "#fc8d59"

##====================================================================================================

## Load Bray-Curtis Distances data
bc <- read.table("Distances/bray_curtis/Description_Distances.txt", fill = T)

## Filter BM vs. BM and TD vs. TD
bc <- cbind(as.numeric(bc[which(bc$V1 == "Buccal.Mucosa_vs._Buccal.Mucosa"),-1]),
            as.numeric(bc[which(bc$V1 == "Tongue.Dorsum_vs._Tongue.Dorsum"),-1]))
colnames(bc) <- c("Buccal_Mucosa", "Tongue_Dorsum")
bc <- bc[complete.cases(bc),]

## Convert to long format
bc <- melt(bc)
colnames(bc) <- c("X1", "X2", "value")

## Compute p-value and prepare a string to add the p-value to the plot
p.bc <- wilcox.test(value ~ X2, data = bc)
p.text.bc <- paste("P = ", round(p.bc$p.value, 2), collapse = "", sep = "")

## Load Weighted Unifrac data
wunifrac <-read.table("Distances/weighted_normalized//Description_Distances.txt", fill = T)

## Filter BM vs. BM and TD vs. TD
wunifrac <- cbind(as.numeric(wunifrac[which(wunifrac$V1 == "Buccal.Mucosa_vs._Buccal.Mucosa"),-1]), 
                  as.numeric(wunifrac[which(wunifrac$V1 == "Tongue.Dorsum_vs._Tongue.Dorsum"),-1]))
colnames(wunifrac) <- c("Buccal_Mucosa", "Tongue_Dorsum")
wunifrac <- wunifrac[complete.cases(wunifrac),]

## Convert to long format
wunifrac <- melt(wunifrac)
colnames(wunifrac) <- c("X1", "X2", "value")

## Compute p-value and prepare a string to add the p-value to the plot
p.wunifrac <- wilcox.test(value ~ X2, data = wunifrac)
p.text.wunifrac <- paste("P = ", round(p.wunifrac$p.value, 2), collapse = "", sep = "")
p.text.wunifrac <- "P <0.001"

## Load Unweighted Unifrac data
unwunifrac <-read.table("Distances/unweighted_unifrac/Description_Distances.txt", fill = T)

## Filter BM vs. BM and TD vs. TD
unwunifrac <- cbind(as.numeric(unwunifrac[which(unwunifrac$V1 == "Buccal.Mucosa_vs._Buccal.Mucosa"),-1]), 
                    as.numeric(unwunifrac[which(unwunifrac$V1 == "Tongue.Dorsum_vs._Tongue.Dorsum"),-1]))
colnames(unwunifrac) <- c("Buccal_Mucosa", "Tongue_Dorsum")
unwunifrac <- unwunifrac[complete.cases(unwunifrac),]

## Convert to long format
unwunifrac <- melt(unwunifrac)
colnames(unwunifrac) <- c("X1", "X2", "value")

## Compute p-value and prepare a string to add the p-value to the plot
p.unwunifrac <- wilcox.test(value ~ X2, data = unwunifrac)
p.text.unwunifrac <- paste("P = ", round(p.unwunifrac$p.value, 3), collapse = "", sep = "")

# Set the plot parameters
par(mfrow = c(1, 3))

## Plot Bray-Curtis data
boxplot(value ~ X2, data = bc, 
        main = "Bray-Curtis",  xlab = "Pairwise distance comparison", ylab = "Distance", names = c("BM vs. BM", "TD vs. TD"),
        vertical = T)
stripchart(value ~ X2, data = bc, pch = 21, bg = c("black"), cex = 1.1,
           col = c(BuccalMucosa_Col, TongueDorsum_Col), vertical = T, add = T)

## Plot Weighted Unifrac data
boxplot(value ~ X2, data = wunifrac, pch = 16,
        main = "Weighted Unifrac", xlab = "Pairwise distance comparison", ylab = "Distance", names = c("BM vs. BM", "TD vs. TD"),
        vertical = T)
stripchart(value ~ X2, data = wunifrac, pch = 21, bg = c("black"), cex = 1.1,
           col = c(BuccalMucosa_Col, TongueDorsum_Col), vertical = T, add = T)

## Plot Unweighted Unifrac data
boxplot(value ~ X2, data = unwunifrac, pch = 16,
        main = "Unweighted Unifrac",  xlab = "Pairwise distance comparison", ylab = "Distance", names = c("BM vs. BM", "TD vs. TD"),
        vertical = T)
stripchart(value ~ X2, data = unwunifrac, pch = 21, bg = c("black"), cex = 1.1,
           col = c(BuccalMucosa_Col, TongueDorsum_Col), vertical = T, add = T)

##====================================================================================================
## B
##====================================================================================================

## Import the biom data
biom_table = import_biom("Rarefied_Tables/rarefaction_1802_0.json.biom", 
                         treefilename = "../Figure1/pruned_rep_set_v35.tree")

## Import the mapping file
map = import_qiime_sample_data("../Figure1/Qiime.map")

# Merge mapping file and biom table
x1 = merge_phyloseq(biom_table, map)

## List the distance methods
dist_methods <- c("bray", "wunifrac", "unifrac")

## Create a list to store ordination analysis for each method
plist <- vector("list", length(dist_methods))
names(plist) = dist_methods

## For each method, perform ordination analysis, get a ggplot object and store in the list created above
for (i in dist_methods) {
  # Compute distance matrix
  iDist <- distance(x1, method = i)
 
  # Compute ordination
  iMDS <- ordinate(x1, "PCoA", distance = iDist)
  
  ## Make plot Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  
  # Create plot, store as temp variable, p
  p <- plot_ordination(x1, iMDS, color = "Description", shape = "Description")
  
  # Add title to each plot
  p <- p + ggtitle(paste("PCoA using distance method ", i, sep = "")) 
  
  # Save the graphic to file.
  plist[[i]] = p
}

# Extract ordination data
df = ldply(plist, function(x) {
  x$data
})

## Rename columns and edit some values
names(df)[1] <- "distance"
df$distance <- as.factor(df$distance)
df$distance <- revalue(df$distance, 
                       c("bray" = "Bray-Curtis",
                         "unifrac" = "Unweighted Unifrac", 
                         "wunifrac" = "Weighted Unifrac"))
df$Description <- revalue(df$Description, 
                          c("Buccal.Mucosa" = "BuccalMucosa",
                            "Tongue.Dorsum" = "TongueDorsum"))

# Plot the figure
p <- ggplot(df, aes(Axis.1, Axis.2, color = Description, shape = Description), groups = Distance) +
  stat_ellipse(type = "t", show.legend = F) +
  scale_color_manual(name = NULL, values=c(BuccalMucosa = BuccalMucosa_Col, TongueDorsum = TongueDorsum_Col), labels = c("Buccal Mucosa", "Tongue Dorsum"))
p <- p + geom_point(size = 3, alpha = 0.9) +
  scale_shape_manual(values=c(BuccalMucosa=16, TongueDorsum=17), labels = c("Buccal Mucosa", "Tongue Dorsum"))
p <- p + facet_wrap(~distance, scales = "free")
p <- p + theme_bw() +
  theme(strip.background = element_rect(colour="black"),
        axis.title = element_text(size = rel(0.8)),
        axis.text.x = element_text(size = rel(0.7)),
        axis.text.y = element_text(size = rel(0.7)),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        legend.key = element_rect(colour = NA),
        legend.key.size = unit(0.5, "cm")) +
  ylab("PC2") + xlab("PC1") +
  guides(title = NULL, shape = guide_legend(override.aes = list(color = c(BuccalMucosa_Col, TongueDorsum_Col))), color = "none")
p
```
