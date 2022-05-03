#### install packages and load libraries ----

# packages <- c("readr", "tidyverse", "matrixStats", "ggplot2", "limma", "edgeR", "DT")
# install.packages(packages)

library(tidyverse) ## for data wraggling
library(readr) ## reads in _tsv and _csv files
library(matrixStats) ##
library(ggplot2) ## visualization
library(limma) ## for 
library(edgeR) ## for
library(DT) ## for interactive table

##### Import data (Combined logFC and master_table) -----

abundance <- read_tsv("combined_abundance.tsv.gz") ##combinelogFC data consisting of geneIDs, sample gene count and tpm
abundanace <- abundance %>% mutate(gene_id = substr(abundance$target_id, 1, 16), .before = "target_id")
colnames(abundance)[2] <- "transcript_id"
rownames(abundance) <- abundance$gene_id

master_table <- read.csv("master_table.csv")  ## read in the master table
master_table <- subset(master_table, genotype.name != "M82") #### removed M82 because it lacks respective control
colnames(master_table)[1] <- "sra_run_id" 
sample_all <- unique(master_table$sample.group) ## extracts all samples and their respective controls
control_group <- unique(master_table$respective_control)[-1] ## extracts only control
sample_group <- sample_all[!(sample_all %in% control_group)] ## extracts only samples


diffGenes.df <- NULL
TopGenes <- NULL

system.time(
  for (i in 1:length(sample_group)){
    tryCatch({
      targets <- rbind((master_table[(master_table$sample.group == sample_group[i]), ## create a table for each sample, its control, sample id
                                     c("sample.name", "sra_run_id", "sample.group")]), 
                       (master_table[(master_table$sample.group == control_group[i]), 
                                     c("sample.name", "sra_run_id", "sample.group")]))
      sampleLabels <- targets$sra_run_id ## extracts the sample ids
      samples <- targets$sample.name 
      myTPM <- as.matrix(abundance[ ,paste0(sampleLabels, "_", "tpm")]) ## extracts tpm for each target sample from abundance file
      myCount <- abundance[ ,paste0(sampleLabels, "_est_counts")] ## extracts count for each target sample
      myDGEList <- DGEList(myCount) ## converts count to a DGElist using (pkg - edgeR)
      cpm <- cpm(myDGEList)
      keepers <- rowSums(cpm > 1) >= (length(samples)/2)
      myDGEList.filtered <- myDGEList[keepers,]
      myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
      log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
      log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
      colnames(log2.cpm.filtered.norm.df) <- c("geneID", samples)
      log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                      cols = samples[1]:samples[length(samples)], # column names to be stored as a SINGLE variable
                                                      names_to = "samples", # name of that new variable (column)
                                                      values_to = "expression") # name of new variable (column) storing all the values (data)
      
      targets$sample <- gsub("-", "_", targets$sample.group)
      group <- factor(targets$sample)
      design <- model.matrix(~0 + group)
      colnames(design) <- levels(group)
      
      v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)  ## add weight and plot
      fit <- lmFit(v.DEGList.filtered.norm, design) ## fit linear model
      
      contrast.matrix <- makeContrasts(contrast = paste0(as.name(colnames(design)[1]) ," - ", as.name(colnames(design)[2])),
                                       levels=design)
      ebFit <- eBayes(contrasts.fit(fit, contrast.matrix))
      myTopHits <- topTable(ebFit, adjust ="fdr", coef=1, p.value=0.05, number=10, sort.by="logFC")
      TopGenes <- rbind(TopGenes, myTopHits)
      results <- decideTests(ebFit, method="global", adjust.method="fdr", p.value=0.05, lfc=1)
      
      diffGenes <- as_tibble(v.DEGList.filtered.norm$E[results[,1] != 0,], rownames = "geneID")
      colnames(diffGenes) <- c("geneID", sampleLabels)
      df <- data.frame(geneID = abundance$gene_id, row.names = abundance$gene_id)
      diffGenes2 <- as.data.frame(t(df %>% left_join(diffGenes, by = "geneID")))
      diffGenes.df <- dplyr::bind_rows(diffGenes.df, diffGenes2[-1,])
    }, error=function(e){})
    print(paste0(i, " is done"))
}
)

diffGenes.df[is.na(diffGenes.df)] <- 0
diffGenes.df <- as.data.frame(t(diffGenes.df))
diffGenes.df2 <- diffGenes.df %>% mutate(Sum = rowSums(as.data.frame(lapply(diffGenes.df, as.numeric))))
rownames(diffGenes.df2) <- row.names(df)
diffGenes.df2 <- subset(diffGenes.df2, Sum != 0)

write.csv(diffGenes.df2, "combined_logFC_18185.csv")

#### Genotype specific heat response -----

##### Load packages
library(RColorBrewer) #need colors to make heatmaps
library(gplots) #the heatmap2 function in this package is a primary tool for making heatmaps
library(gameofthrones) #because...why not.  Install using 'devtools::install_github("aljrico/gameofthrones")'
library(heatmaply) #for making interactive heatmaps using plotly
library(d3heatmap) #for making interactive heatmaps using D3
library(tidyverse)

#### Genotypes in pollens -----

#### import data - master_table and combine_logFC and 
master_table <- read.csv("master_table.csv")  ## read in the master table
master_table <- subset(master_table, genotype.name != "M82")
master_table <- subset(master_table, tissue == "pollen")
master_table$genotype.name <- gsub(" ", "_", master_table$genotype.name)
colnames(master_table)[1] <- "sra_run_id"  ## name sample/sra
sample_all <- unique(master_table$sample.group)
genotypes <- unique(master_table$genotype.name)
control_group <- unique(master_table$respective_control)[-1]
sample_group <- sample_all[!(sample_all %in% control_group)]

sra_sample_group <- master_table[unique(master_table$sample.group == sample_group), ]$sra_run_id

sra_control_group <- master_table[unique(master_table$sample.group == control_group),]$sra_run_id

all_sra <- c(sra_control_group, sra_sample_group)
genotypes <- master_table[unique(master_table$sra_run_id == all_sra), ]$genotype.name

combine_logFC <- read.csv("combined_logFC_18185.csv")
rownames(combine_logFC) <- combine_logFC$X

selected_logFC <- combine_logFC[, names(combine_logFC) %in% all_sra] ### logFC for all selected pollen samples 

names <- data.frame(sra_run_id = colnames(selected_logFC)) %>% left_join(master_table, by = "sra_run_id")
colnames(selected_logFC) <- ifelse(colnames(selected_logFC) %in% names$sra_run_id, 
                             names[unique(names$sra_run_id %in% colnames(selected_logFC)), ]$genotype.name,
                             names$sra_run_id)
## remove redundant genes (gene with logFC = 0 for these selected genes)

selected_logFC2 <- data.frame(lapply(selected_logFC, as.numeric))

rownames(selected_logFC2) <- rownames(selected_logFC)

selected_logFC2 <- selected_logFC2 %>% mutate(Sum = rowSums(selected_logFC2))

selected_logFC_filtered <- subset(selected_logFC2, Sum != 0)
selected_logFC_filtered <- subset(selected_logFC_filtered, select = -c(Sum))

#colnames(selected_logFC) <- colnames_selected_genotypes
t_selected_logFC_filtered <- t(selected_logFC_filtered)
clustRows <- hclust(as.dist(1-cor(t_selected_logFC_filtered, method="pearson")), method="complete") 
 
# clustColumns <- as.dist(1-cor(selected_logFC_filtered, method="spearman")) #cluster columns by spearman correlation)
library(poppr)

clustColumns <- hclust(rogers.dist(1 - as.matrix(selected_logFC_filtered)))


module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 

myheatcolors <- colorRampPalette(colors=c("blue","white","red", "green"))(100)

myheatcolors <- colorRampPalette(colors=c("blue","white","red"))(100)

heatmap.2(as.matrix(selected_logFC_filtered), 
          Rowv=as.dendrogram(clustRows), 
          #Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=rev(myheatcolors), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20)) 

#### Accross tissue ----

master_table <- read.csv("master_table.csv")  ## read in the master table
master_table <- subset(master_table, genotype.name == "Moneymaker")

tissue_table = NULL
tissue <- c("pollen", "leaf", "stem", "anther")

for (i in 1:length(tissue)) {
  tissue_tab <- master_table[master_table$tissue == tissue[i], ]
  tissue_table <- rbind(tissue_table, tissue_tab)
}

colnames(tissue_table)[1] <- "sra_run_id"  ## name sample/sra

tissue_sra <- tissue_table$sra_run_id

combine_logFC <- read.csv("combined_logFC_all.csv")
rownames(combine_logFC) <- combine_logFC$X

selected_logFC <- combine_logFC[, names(combine_logFC) %in% tissue_sra]

names <- data.frame(sra_run_id = colnames(selected_logFC)) %>% left_join(tissue_table, by = "sra_run_id")
colnames(selected_logFC)  <- ifelse(colnames(selected_logFC) %in% names$sra_run_id, 
                                    names[unique(names$sra_run_id %in% colnames(selected_logFC)), ]$tissue,
                                    names$sra_run_id)
## remove redundant genes (gene with logFC = 0 for these selected genes)

selected_logFC2 <- data.frame(lapply(selected_logFC, as.numeric))

rownames(selected_logFC2) <- rownames(selected_logFC)

selected_logFC2 <- selected_logFC2 %>% mutate(Sum = rowSums(selected_logFC2))

selected_logFC_filtered <- subset(selected_logFC2, Sum != 0)
selected_logFC_filtered <- subset(selected_logFC_filtered, select = -c(Sum))

#colnames(selected_logFC) <- colnames_selected_genotypes
t_selected_logFC_filtered <- t(selected_logFC_filtered)
clustRows <- hclust(as.dist(1-cor(t_selected_logFC_filtered, method="pearson")), method="complete") 

# clustColumns <- as.dist(1-cor(selected_logFC_filtered, method="spearman")) #cluster columns by spearman correlation)
#library(poppr)

#clustColumns <- hclust(rogers.dist(1 - as.matrix(selected_logFC_filtered)))


module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 

myheatcolors <- colorRampPalette(colors=c("blue","white","red"))(100)

heatmap.2(as.matrix(selected_logFC_filtered), 
          Rowv=as.dendrogram(clustRows), 
          #Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=rev(myheatcolors), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20)) 


##### Top DEGs annontation ----

library(gt)
library(DT)
library(plotly)
library(gprofiler2)

top_genes <- read.csv("Top10GenesAcross.csv") ### read in top_ten genes across experiment

colnames(top_genes)[1] <- "geneID"

# visualize the distribution of the top_genes
# vplot <- ggplot(top_genes) +
#   aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
#   geom_point(size=2.5) +
#   geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
#   geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
#   #geom_vline(xintercept = -10, linetype="longdash", colour="#2C467A", size=1) +
#   annotate("rect", xmin = 1, xmax = 15, ymin = -log10(0.05), ymax = 12, alpha=.2, fill="#BE684D") +
#   #annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
#   labs(title="Volcano plot (log10FC >= 1, FDR <= 0.05)",
#        #subtitle = "Cutaneous leishmaniasis",
#        #caption=paste0("produced on ", Sys.time())
#        ) +
#   theme_bw()
# vplot

top_upregulated <- subset(top_genes, logFC > 1 & adj.P.Val < 0.05) ## fi

gt(top_upregulated)  ### fetch out interactive table of top 151 genes

gost.res <- gost(top_upregulated$geneID, correction_method = "fdr", significant = TRUE, organism = 'slycopersicum')

gostplot(gost.res, interactive = T, capped = T)
mygostplot <- gostplot(gost.res, interactive = F, capped = T)

publish_gostplot(
  mygostplot, #your static gostplot from above
  highlight_terms = c("GO: 00031072", "GO:0006950", "GO:0009266", "GO:0009408", "GO:0009628", 
                      "GO:0009651", "GO:0006970", "GO:0006979", "GO:0006950", "GO:0034605"),
  filename = NULL,
  width = NA,
  height = NA)

