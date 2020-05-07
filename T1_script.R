library(pamr)
library(plyr)
library(pheatmap)
library(WGCNA)
library(biomaRt)
library(tidyr)
library(CMScaller)

CollectCollapsedData <- function(base.dir)
{
  
  file.path <- paste0(base.dir, "Alk5ca_WT_mouse_collapsed.txt")
  CollapsedData <- read.delim(file.path, stringsAsFactors=FALSE, sep='\t')
  
  return(CollapsedData)
}

data.dir <- "T1_project/Alk5_TGFb_MouseData/T1.paper_Scripts/"
df <- CollectCollapsedData(data.dir)

##############################  PAMR analysis
## Preparing data
df.2 <- df
rownames(df.2) <- c()
df.2 <- column_to_rownames(df.2, var = "selectedRowID")[,-1]

x <- as.matrix(df.2)
y <- factor(c(rep("WT",3), rep("Alk5ca", 3)))

mydata <- list(x=x, y=y, geneid=as.character(rownames(df.2)), genenames=as.character(df$group))

set.seed(367703)
#train classifier
Alk5_train <- pamr.train(mydata)

## Cross-validate the classifier
Alk5_cv <- pamr.cv(Alk5_train, mydata)
Alk5_cv
## Plot the cross-validated error curves
pamr.plotcv(Alk5_cv)

## A function to plot the genes that survive the thresholding, from the nearest shrunken centroid classifier produced by pamr.train
par("mar")
par(mar=c(1,1,1,1))
pamr.geneplot(Alk5_train, mydata, 4.98)

## Compute the confusion matrix for a particular model (threshold=4.0) 
pamr.confusion(Alk5_cv, threshold=4.98)


GeneList <- pamr.listgenes(Alk5_train, mydata, 4.98, genenames = TRUE)

GeneList <- as.data.frame(GeneList)



### Convert the GeneList to human orthologs using Biomart 
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mouse.genes = GeneList$name
human.orthologs = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mouse.genes ,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)


### Merging GeneList and human.orthologs
colnames(human.orthologs)[1] <- "name"
j <- join(GeneList, human.orthologs, by = "name", type = "inner")
colnames(j)[3] <- "Alk5ca.score"

#### Preparing NTP template

s1 <- subset(j, Alk5ca.score > 0)
dim(s1)
s2 <- subset(j, Alk5ca.score < 0)
dim(s2)

class = c(rep("Alk5ca", 16), rep("WT", 37))

probe = c(s1$HGNC.symbol, s2$HGNC.symbol)

NTP_template <- data.frame(class, probe)

##### use NTP template to classify T1 samples 

## read T1 data 
CollectCollapsed.T1.Data <- function(base.dir)
{
  
  file.path <- paste0(base.dir, "T1_collapseGenes.txt")
  Collapsed.T1.Data <- read.delim(file.path, stringsAsFactors=FALSE, sep='\t')
  
  return(Collapsed.T1.Data)
}

data.dir <- "T1_project/Microarray_Data/Feasibility cohort/New_analysis_14.04.2020/just_Maxmean/"
T1 <- CollectCollapsed.T1.Data(data.dir)
T1 <- T1[,-c(1,2)]

# adjust matrix before classification
mat_exp_adj <- ematAdjust(T1)

NTP_res <- ntp(emat=mat_exp_adj,NTP_template, seed=367707,doPlot = TRUE)

#################### Visualize the prediction result using pheatmap 

#### replace Alk5ca & WT with 2 & 1 respectively
NTP_res.copy <- NTP_res
## convert factor to character
NTP_res.copy %>% mutate_if(is.factor, as.character) -> NTP_res.copy

NTP_res.copy[NTP_res.copy$prediction == "Alk5ca",][,1] <- 2
NTP_res.copy[NTP_res.copy$prediction == "WT",][,1] <- 1

rownames(NTP_res.copy) <- c( "non-Relapse1", "Relapse1", paste0("non-Relapse", 2:6), paste0("Relapse", 2:3),
                             paste0("non-Relapse", 7:13), paste0("Relapse", 4:7), paste0("non-Relapse", 14:15), "Relapse8", "non-Relapse16",
                             paste0("Relapse", 9:10), "non-Relapse17")

NTP_res.copy[,1] <- as.numeric(NTP_res.copy[,1])

### plot
pheatmap(NTP_res.copy[,-c(3,4,5)], fontsize = 12, angle_col = 0, cluster_cols = FALSE)

