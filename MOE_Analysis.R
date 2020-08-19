#title: Filtering, Normalization and MOE data analysis
#Mohammad Hussainy 17.08.2020

rm(list=ls())
library(RColorBrewer)
library(Seurat)
library(slingshot)
library(dplyr)
library(data.table)
library(ggpubr)

# load data and corrected UMI counts
colr <- colorRampPalette(brewer.pal(12,"Paired"))
load("~/SingleCellSeqR/OlfrClusters/Seurat/_filtdata.Rda")
MOEdata <- CreateSeuratObject(counts, project = "OE", min.cells = 1,meta.data = cbind(batch = batch, expt = expt),min.features = 1000)
MOEdata <- SCTransform(MOEdata,variable.features.n = 4000)

# adding the top10 expressed genes and OR clusters (named by their Greek islands) expression matrix  

# adding top10 OR genes
data <- as.matrix(MOEdata[["SCT"]]@counts)
all.genes <- rownames(counts)
olfind <- grep(pattern = "Olfr", x = all.genes)
olfgenes <- rownames(counts)[olfind]
counts <- counts[olfgenes,colnames(data)]
counts <- counts[!rownames(counts) %in% c("Olfr856-ps1","Olfr613"),]
multigenic <- CreateSeuratObject(counts,min.cells = 1,min.features = 5)
allOlfs <- counts
sortall = apply(allOlfs,2,sort,decreasing=T)
ors <- sortall[1:10,]
rownames(ors) <- c("OR1","OR2","OR3","OR4","OR5","OR6","OR7","OR8","OR9","OR10")

data <- data[!rownames(data) %in% olfgenes,]
counts <- rbind(ors,data,allOlfs)

## Change ORs names to their corresponding enhancers (Greek islands) as they clustered in Monahana et al., 2017
myORsgenes <- read.csv(file = "~/SingleCellSeqR/bedtools/myORsgenes_Greek_Islands_Forbed2_All.xls",stringsAsFactors = F)
inter <- intersect(myORsgenes$name2,rownames(allOlfs))
mygreek <- setDT(myORsgenes,key = "name2")[J(inter)]
allgi <- allOlfs[inter,]
rownames(allgi) <- mygreek$Greek_Island
x <- vector()
y <- vector()     
for (i in 1:nrow(allgi)) {
        x[i] <- grep(substr(rownames(allgi)[i],1,4),x = rownames(allgi))
        y[i] <- substr(rownames(allgi)[i],1,4)
}
gi <- rowsum(allgi,group = x)
rownames(gi) <- unique(y)
#change the name of OR cluster that has the same name with other gene in our data
dub <- intersect(rownames(counts),rownames(gi))
rownames(gi)[rownames(gi)%in%dub] <- paste0(dub,1,"")

counts <- rbind(gi,counts)
dim(counts)

------------------------------------------
        ## Seurat Clustring and cell types identification
#Construct Seurat object and calculate the relative expression of OR in each cell
colr <- colorRampPalette(brewer.pal(12,"Paired"))
MOEdata <- CreateSeuratObject(counts = counts[rownames(counts),], project = "OE", min.cells = 1,min.features = 300)
MOEdata <- PercentageFeatureSet(MOEdata, pattern = "Olfr", col.name = "RC.Olfr")
MOEdata

## extract Olfr genes in separate seurat object
# all.genes <- rownames(olfdata)
# olfind <- grep(pattern = "Olfr", x = all.genes)
# olfgenes <- rownames(olfdata)[olfind]
# olfseurat <- olfdata[olfgenes,]

#Cell cycle genes and cc pahse of cells
findGenes.GO_Name <- function(GO.Name,your.genes){  #My function to find the cell cycle genes
        require(biomaRt)
        require(dplyr)
        ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
        
        GO.Genes <- getBM(attributes = c("external_gene_name","entrezgene_id","ensembl_gene_id", "chromosome_name","start_position","end_position","external_synonym","affy_mg_u74a"),filters = "go_parent_name",values = GO.Name,mart = ensembl)
        GO.Genes.In.YourGenes <- intersect(GO.Genes$external_gene_name,your.genes)
}
g1s <- findGenes.GO_Name(GO.Name = "G1/S transition of mitotic cell cycle",your.genes = rownames(MOEdata))
g2m <- findGenes.GO_Name(GO.Name = "G2/M transition of mitotic cell cycle",your.genes = rownames(MOEdata))

MOEdata <- CellCycleScoring(MOEdata, s.features = g1s, g2m.features = g2m)

#Normalization and scaling data, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques
MOEdata <- NormalizeData(MOEdata, normalization.method = "LogNormalize", scale.factor = 10000)
MOEdata <- FindVariableFeatures(MOEdata, selection.method = "vst", nfeatures = 4000)
MOEdata <- ScaleData(MOEdata, features = VariableFeatures(MOEdata))#, vars.to.regress = c("S.Score", "G2M.Score"))
MOEdata

# Identify the 10 most highly variable genes
top10OE <- head(VariableFeatures(MOEdata), 10)
plot1 <- VariableFeaturePlot(MOEdata)
plot2 <- LabelPoints(plot = plot1, points = top10OE, repel = TRUE)
plot2

##    Perform linear dimensional reduction
MOEdata <- RunPCA(MOEdata, features = VariableFeatures(MOEdata))
print(MOEdata[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(MOEdata, reduction = "pca")

##    Determine the ‘dimensionality’ of the dataset (selected top 15 PCs as source for heterogenity of data)
ElbowPlot(MOEdata)
#jackStraw method
#olfdata <- JackStraw(olfdata,dims = 20)
#ScoreJackStraw(olfdata,dims = 1:20,do.plot = T)

##   Cluster the cells
#KNN clustering using the top 15 PCs, and 14 nearest neighbors
MOEdata <- FindNeighbors(MOEdata, dims = 1:15,k.param = 14)
MOEdata <- FindClusters(MOEdata, resolution = 1.5)

#kmeans clustering
# clusterkmeans <- kmeans(olfdata[["pca"]]@cell.embeddings[,1:15], centers = 13,nstart = 5)
# clusterkmeans[["size"]]
# olfdata@meta.data$kmeans <- factor(clusterkmeans$cluster)

##     Run non-linear dimensional reduction (UMAP/tSNE)
MOEdata <- RunUMAP(MOEdata, dims = 1:15)
MOEdata <- RunTSNE(MOEdata, dims = 1:15)

## Visualization of cells
DimPlot(MOEdata, reduction = "umap", label = T)
DimPlot(MOEdata,reduction = "tsne",label = T)
DimPlot(MOEdata,reduction = "pca",label = T)

##  Using Known markers to identify the cell types
## HBC cells markers
VlnPlot(MOEdata,features = c("Trp63","Krt5","Krt14"))
## Sus cells markers
VlnPlot(MOEdata,features = c("Cyp2g1","Notch2","Cyp1a2","Sox2","Hey1"))
## MV cells markers
VlnPlot(MOEdata,features = c("Ascl3","Cftr","Coch"))
## qGBC cells markers
VlnPlot(MOEdata,features = c("Sox9","Cd44","Tmprss4","Kit","Hes6","Lgr5"))
## GBC cells markers
VlnPlot(MOEdata,features = c("Sox2","Ascl1","Mki67","Top2a"))
## INP.Initial cells markers 
VlnPlot(MOEdata,features = c("Neurod1","Neurog1","Top2a","Mki67"))
## INP.Early cells markers 
VlnPlot(MOEdata,features = c("Neurod1","Top2a","Lhx2"),assay = "alra") # no cell cycle genes
## INP.Mid cells markers 
VlnPlot(MOEdata,features = c("Lhx2","Ebf1","Gng8","Gap43"))
## iOSN cells markers
VlnPlot(MOEdata,features = c("Gng8","Gap43","Trib3"))
## mOSN cells markers
VlnPlot(MOEdata, features = c("Omp","Gnal","Gng13"))

DotPlot(MOEdata,features = c("Trp63","Krt5","Krt14","Cyp2g1","Notch2","Cyp1a2","Hey1","Sox2","Ascl3","Cftr","Coch","Sox9","Hes1","Hes6","Lgr5","Tmprss4","Kit",
                             "Ascl1","Mki67","Top2a","Neurod1","Neurog1","Lhx2","Gng8","Ebf1","Trib3","Gap43","Gng13","Gnal","Omp"),
        col.min = 0,cols = brewer.pal(12,"Paired")[c(1,6)])+theme(axis.text.x = element_text(angle = 45,hjust = 1))

# renamed clusters according to their cell types
new.cluster.ids <- c("HBC2","iOSN","HBC0", "iSus","qGBC","mSus","GBC","MV","INP.Mid","mOSN","HBC1","INP.Initial","INP.Early")

names(new.cluster.ids) <- levels(MOEdata)
MOEdata <- RenameIdents(MOEdata, new.cluster.ids)
MOEdata@active.ident <- factor(MOEdata@active.ident , levels = c("HBC0", "HBC1","HBC2","iSus","mSus", "MV","qGBC","GBC", "INP.Initial", "INP.Early","INP.Mid", "iOSN", "mOSN"))
MOEdata@meta.data$seurat_clusters <- MOEdata@active.ident 

DimPlot(MOEdata,reduction = "pca",pt.size = 1.5,cols = colr(20),repel = T,label = T)+
        #ggtitle("MOE Cell Populations")+
        theme(plot.title = element_text(hjust = 0.5))
# Cell types Known Markers
DotPlot(MOEdata,features = c("Trp63","Krt5","Krt14","Cyp2g1","Notch2","Cyp1a2","Hey1","Sox2","Ascl3","Cftr","Coch","Sox9","Hes1","Hes6","Lgr5","Tmprss4","Kit",
                             "Ascl1","Mki67","Top2a","Neurod1","Neurog1","Lhx2","Gng8","Ebf1","Trib3","Gap43","Gng13","Gnal","Omp"),
        col.min = 0,cols = brewer.pal(12,"Paired")[c(1,6)])+theme(axis.text.x = element_text(angle = 45,hjust = 1))
--------------------------------------------------------------------

        #Library size of corrected-UMI and before dealing with dropouts  
MOEdata[["nGenes_per_Cell"]] <- MOEdata@meta.data$nFeature_RNA
MOEdata[["nCounts_per_Cell"]] <- MOEdata@meta.data$nCount_RNA
ggarrange(VlnPlot(MOEdata,"nGenes_per_Cell",cols = colr(20))+theme_pubclean(),
          VlnPlot(MOEdata,"nCounts_per_Cell",cols = colr(20))+theme_pubclean()
          ,labels = c("a","b"),nrow = 2,common.legend = T)

#Retrieving the missing values (dropouts) in scRNA-seq data (corrected nGenes detected per cell)
MOEdata <- RunALRA(MOEdata,assay = "RNA",genes.use = rownames(MOEdata),setDefaultAssay = F)
# Distribution of genes per cell after manipulation of dropouts
library(AUCell)
normCounts <- as.matrix(MOEdata[["alra"]]@data)
cells_rankings <- AUCell_buildRankings(normCounts,plotStats = T)


#find differentially expressed genes for identified clusters using Wilcoxon Rank Sum test 
DE.genes <- FindAllMarkers(MOEdata, min.pct = 0.3, logfc.threshold = 0.5,only.pos = T,
                                  assay = "RNA")

markers <- DE.genes[which(DE.genes$p_val_adj < 0.05 & DE.genes$avg_logFC > 1),]
top10O <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#heatmap plot for top 10 DE genes in each cell type
ggarrange(DoHeatmap(MOEdata, features = top10O$gene,group.colors = colr(20)
                    ,assay = "RNA",combine = T,slot = "data",label = T,raster = F)+
                  scale_fill_gradientn(colors = c("white","blue","red"))+
                  theme(plot.title = element_text(colour = "darkblue",size = 21,hjust = 0.5,face = "bold.italic"),
                        legend.background = element_rect(fill = "grey",)),
          common.legend = T,legend = "top")+theme_pubclean()

#Markers for each cell type
HBC0 <- markers[which(markers$cluster=="HBC0"),]
HBC1 <- markers[which(markers$cluster=="HBC1"),]
HBC2 <- markers[which(markers$cluster=="HBC2"),]
iSus <- markers[which(markers$cluster=="iSus"),]
mSus <- markers[which(markers$cluster=="mSus"),]
MV <- markers[which(markers$cluster=="MV"),]
qGBC <- markers[which(markers$cluster=="qGBC"),]
GBC <- markers[which(markers$cluster=="GBC"),]
Initial <- markers[which(markers$cluster=="INP.Initial"),]
Early <- markers[which(markers$cluster=="INP.Early"),]
Mid <- markers[which(markers$cluster=="INP.Mid"),]
iOSN <- markers[which(markers$cluster=="iOSN"),]
mOSN <- markers[which(markers$cluster=="mOSN"),]

#Using AUC to test the gene signature of each cell type
cellsAUC <- AUCell_calcAUC(geneSets = list(HBC0.Signature=unique(HBC0$gene),HBC1.Signature=unique(HBC1$gene),HBC2.Signature=unique(HBC2$gene),
                                           iSus.Signature=unique(iSus$gene),mSus.Signature=unique(mSus$gene),MV.Signature=unique(MV$gene),
                                           qGBC.Signature=unique(qGBC$gene),GBC.Signature=unique(GBC$gene),Initial.Signature=unique(Initial$gene),
                                           Early.Signature=unique(Early$gene),Mid.Signature=unique(Mid$gene),Immature.Signature=unique(iOSN$gene),
                                           Mature.Signature=unique(mOSN$gene),ORs= rownames(allOlfs),OR.Clus=rownames(gi)),rankings = cells_rankings)
#plot the area under curve analysis for each gene signature
par(mfrow=c(3,5))
AUC_thres <- AUCell_exploreThresholds(cellsAUC, th = 0.05, nCores = 1,
                                      smallestPopPercent = 0.5, plotHist = TRUE, densAdjust = 2,
                                      assignCells = T, nBreaks = 100, verbose = TRUE)
selectedThresholds <- getThresholdSelected(AUC_thres)

# 0.5 threshold
selectedThresholds["Early.Signature"] <- 0.5
selectedThresholds["Mid.Signature"] <- 0.5
selectedThresholds["Immature.Signature"] <- 0.5
selectedThresholds["Mature.Signature"] <- 0.5

par(mfrow=c(4,4))
AUCell_plotTSNE(tSNE= MOEdata[["tsne"]]@cell.embeddings , exprMat=normCounts, 
                cellsAUC=cellsAUC[10,], thresholds=selectedThresholds)
plot(MOEdata[["tsne"]]@cell.embeddings,
     col = c(rep("lightblue",9),"red",rep("lightblue",3))[MOEdata@active.ident],main="Early", asp = 1, pch = 16,las=1,xlim=c(-35,40),bty="n")
AUCell_plotTSNE(tSNE= MOEdata[["tsne"]]@cell.embeddings , exprMat=normCounts, 
                cellsAUC=cellsAUC[11,], thresholds=selectedThresholds)
plot(MOEdata[["tsne"]]@cell.embeddings,
     col = c(rep("lightblue",10),"red",rep("lightblue",2))[MOEdata@active.ident],main="Mid", asp = 1, pch = 16,las=1,xlim=c(-35,40),bty="n")
AUCell_plotTSNE(tSNE= MOEdata[["tsne"]]@cell.embeddings , exprMat=normCounts, 
                cellsAUC=cellsAUC[12,], thresholds=selectedThresholds)
plot(MOEdata[["tsne"]]@cell.embeddings,
     col = c(rep("lightblue",11),"red", rep("lightblue",1))[MOEdata@active.ident],main="Immature", asp = 1, pch = 16,las=1,xlim=c(-35,40),bty="n")
AUCell_plotTSNE(tSNE= MOEdata[["tsne"]]@cell.embeddings , exprMat=normCounts, 
                cellsAUC=cellsAUC[13,], thresholds=selectedThresholds)
plot(MOEdata[["tsne"]]@cell.embeddings,
     col = c(rep("lightblue",12),"red",rep("lightblue",0))[MOEdata@active.ident],main="Mature", asp = 1, pch = 16,las=1,xlim=c(-35,40),bty="n")


saveRDS(MOEdata, file = "~/MOEdata.rds")


#Trajectory inference
slings <- slingshot(MOEdata@reductions[["pca"]]@cell.embeddings[,1:10], MOEdata@active.ident,start.clus="qGBC")
slings@lineages
lineages <- getLineages(MOEdata[["pca"]]@cell.embeddings[,1:10], clusterLabels = MOEdata@active.ident, start.clus= 'qGBC', end.clus = c('mOSN','MV','mSus',"HBC0"))
par(mfrow=c(1,1))
plot(MOEdata[["pca"]]@cell.embeddings[,1:15], col = colr(20)[MOEdata@active.ident], asp = 1, pch = 16)
lines(lineages, lwd = 3, col = 'black', show.constraints = TRUE)
legend(-28,-11, legend=levels(MOEdata@active.ident),col = colr(20),ncol = 2,text.width=2,y.intersp = 0.6,x.intersp = 0.5,
       title="Cell Types", text.font=3, pch = 19,bty = "n")#,title.adj = 0.35
legend(13,-15, legend=c("Start","End"),text.width=2,
       fill = c("black","black"),col = c("green","red"),
       title="Lineage", text.font=3, pch = c(19,19),bty = "n") #1, -13


##VolcanoPlot of possible key genes of the four detected lineages
#I modified EnhancedVolcano function obtained from https://github.com/kevinblighe/EnhancedVolcano
VolcanoPlot <- function(tf.markers, AdjustedCutoff=0.05, LabellingCutoff=0.05, FCCutoff=1, main="VolcanoPlot")
{ 
        tf.markers$Significance <- "NS"
        tf.markers$Significance[(abs(tf.markers$avg_logFC) > FCCutoff)] <- "FC"
        tf.markers$Significance[(tf.markers$p_val_adj<AdjustedCutoff)] <- "p_adjusted"
        tf.markers$Significance[(tf.markers$p_val_adj<AdjustedCutoff) & (tf.markers$avg_logFC>= FCCutoff)] <- "Upregulated"
        tf.markers$Significance[(tf.markers$p_val_adj<AdjustedCutoff) & (tf.markers$avg_logFC<= -FCCutoff)] <- "Downregulated"
        table(tf.markers$Significance)
        tf.markers$Significance <- factor(tf.markers$Significance, levels=c("NS", "FC", "p_adjusted", "Upregulated", "Downregulated"))
        
        
        plot <- ggplot(tf.markers, aes(x=avg_logFC, y=-log10(p_val_adj))) +
                
                #Add points:
                #   Colour based on factors set a few lines up
                #   'alpha' provides gradual shading of colour
                #   Set size of points
                geom_point(aes(color=factor(Significance)), alpha=1/2, size=2.5) +
                
                #Choose which colours to use; otherwise, ggplot2 choose automatically (order depends on how factors are ordered in toptable$Significance)
                scale_color_manual(values=c(NS="black", FC="grey", p_adjusted ="royalblue", Upregulated="green3", Downregulated = "red2" ), labels=c(NS="NS", FC=paste("LogFC>|", FCCutoff, "|", sep=""), p_adjusted=paste("p_adjusted<", AdjustedCutoff, sep=""), Upregulated=paste("Upregulated"),Downregulated=paste("Downregulated"))) +
                
                #Set the size of the plotting window
                theme_bw(base_size=14) +
                
                #Modify various aspects of the plot text and legend
                theme(legend.background=element_rect(),
                      plot.title=element_text(angle=0, size=16, face="bold", vjust=1),
                      
                      panel.grid.major=element_blank(), #Remove gridlines
                      panel.grid.minor=element_blank(), #Remove gridlines
                      
                      axis.text.x=element_text(angle=0, size=16, vjust=1),
                      axis.text.y=element_text(angle=0, size=16, vjust=1),
                      axis.title=element_text(size=16),
                      
                      #Legend
                      legend.position="top",            #Moves the legend to the top of the plot
                      legend.key=element_blank(),       #removes the border
                      legend.key.size=unit(0.5, "cm"),  #Sets overall area/size of the legend
                      legend.text=element_text(size=16), #Text size
                      title=element_text(size=16),       #Title text size
                      legend.title=element_blank()) +       #Remove the title
                
                #Change the size of the icons/symbols in the legend
                guides(colour = guide_legend(override.aes=list(size=5))) +
                
                #Set x- and y-axes labels
                xlab(bquote(~Log[2]~ "fold change")) +
                ylab(bquote(~-Log[10]~adjusted~italic(P))) +
                
                #Set the axis limits
                xlim(-4, 4) +
                # ylim(0, 20) +
                
                #Set title
                ggtitle(main) +
                
                #Tidy the text labels for a subset of genes
                geom_text(data=subset(tf.markers, p_val_adj<LabellingCutoff & abs(avg_logFC)>FCCutoff),
                          aes(label=rownames(subset(tf.markers, p_val_adj<LabellingCutoff & abs(avg_logFC)>FCCutoff))),
                          size=4,
                          segment.color="black", #This and the next parameter spread out the labels and join them to their points by a line
                          segment.size=0.01,
                          check_overlap=TRUE,
                          vjust=1.0) +
                
                # geom_text(data=subset(tf.markers[[1]], abs(avg_logFC)>FCCutoff),
                #           aes(label=rownames(subset(tf.markers[[1]], abs(avg_logFC)>FCCutoff))),
                #           size=2.25,
                #           #segment.color="black", #This and the next parameter spread out the labels and join them to their points by a line
                #           #segment.size=0.01,
                #           check_overlap=TRUE,
                #           vjust=1.0) +
                #Add a vertical line for fold change cut-offs
                geom_vline(xintercept=c(-FCCutoff, FCCutoff), linetype="longdash", colour="black", size=0.4) +
                
                #Add a horizontal line for P-value cut-off
                geom_hline(yintercept=-log10(AdjustedCutoff), linetype="longdash", colour="black", size=0.4)
        
        return(plot)
        
}

neuronal_lineage <- FindMarkers(MOEdata,ident.1 = "mOSN",ident.2 = "qGBC",logfc.threshold = 0.1,min.pct = 0.3)
Microvellous_lineage <- FindMarkers(MOEdata,ident.1 = "MV",ident.2 = "qGBC",logfc.threshold = 0.1,min.pct = 0.3)
Sustentacular_lineage <- FindMarkers(MOEdata,ident.1 = "iSus",ident.2 = "qGBC",logfc.threshold = 0.1,min.pct = 0.3)
HBC_lineage <- FindMarkers(MOEdata,ident.1 = "HBC0",ident.2 = "qGBC",logfc.threshold = 0.1,min.pct = 0.3)
tf_neurolineage <- list(OSN_lineage= neuronal_lineage,MV_lineage = Microvellous_lineage, Sus_lineage = Sustentacular_lineage, HBC_lineage = HBC_lineage)

plotlist <- list()
namesstep <- names(tf_neurolineage)
for(i in 1:length(tf_neurolineage)){
        plotlist[[i]] <- VolcanoPlot(tf_neurolineage[[i]], main = namesstep[i]) 
        plot(plotlist[[i]])
}
