library(ggplot2)
library(dplyr)
library(Seurat)
packageVersion("Seurat")
options(stringsAsFactors=F)
#3  files in this folder
data <- Read10X("C:\\Users\\wilso\\Desktop\\scrna\\SRR9036396_filtered_feature_bc_matrix")
train <- CreateSeuratObject(counts = data,project = "train", min.cells = 3, min.features = 200)

expr_matrix <- train[["RNA"]]@counts
write.table(expr_matrix,file="C:\\Users\\wilso\\Desktop\\scrna\\train.UMI.counts.xls",
            col.names=T,row.names=T,quote=F,sep="\t")

train[["percent.mt"]] <- PercentageFeatureSet(train, pattern = "^MT-")
head(train@meta.data)

#QC violin plot
pdf("C:\\Users\\wilso\\Desktop\\scrna\\train.cellqc.gene.pdf")
VlnPlot(train, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#discard low-quality cells and QC
train
train <- subset(train, subset = percent.mt < 10 & nFeature_RNA >= 250 & nFeature_RNA  < 3000)
train
plot1 <- FeatureScatter(train, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(train, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("C:\\Users\\wilso\\Desktop\\scrna\\train.cellqc.scatter.pdf")
plot1 + plot2
dev.off()

#normalize
train <- NormalizeData(train, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
head(train[["RNA"]]@data[,1:5])

#plot high variable genes
train <- FindVariableFeatures(train, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(train), 10)
plot1 <- VariableFeaturePlot(train)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf("C:\\Users\\wilso\\Desktop\\scrna\\train.hvg.pdf")
plot2
dev.off()

#scale
train <- ScaleData(train,features=rownames(train))
head(train[["RNA"]]@scale.data[,1:5])

#PCA
train <- RunPCA(object=train,features=VariableFeatures(train),npcs=50)
print(train[["pca"]], dims = 1:5, nfeatures = 5)
pdf("C:\\Users\\wilso\\Desktop\\scrna\\train.pca.vizdim.pdf")
VizDimLoadings(train, dims = 1:5, nfeatures = 10, reduction = "pca")
dev.off()
pdf("C:\\Users\\wilso\\Desktop\\scrna\\train.pca.dimplot.pdf")
DimPlot(train, reduction = "pca",label = T)
dev.off()

#determine significant dimensions
pdf("C:\\Users\\wilso\\Desktop\\scrna\\train.pca.heatmap.pdf")
DimHeatmap(train, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()
train <- JackStraw(train, reduction="pca",num.replicate = 50,prop.freq=0.01)
train<- ScoreJackStraw(train, dims = 1:20)
pdf("C:\\Users\\wilso\\Desktop\\scrna\\train.JackStrawPlot.pdf")
JackStrawPlot(train, dims = 1:20)
dev.off()
pdf("C:\\Users\\wilso\\Desktop\\scrna\\train.ElbowPlot.pdf")
ElbowPlot(train)
dev.off()

train <- FindNeighbors(train, dims = 1:10)
train <- FindClusters(train, resolution = 0.5)
cellcluster <- train@meta.data
cellcluster$cellid <- rownames(cellcluster)
cellcluster <- subset(cellcluster,select=c("cellid","seurat_clusters"))
write.table(cellcluster,file="C:\\Users\\wilso\\Desktop\\scrna\\train.cell.cluster.xls",
            quote=F,col.names=T,row.names=F,sep="\t")

#TSNE AND UMAP
train <- RunTSNE(object=train, dims=1:30)
pdf("C:\\Users\\wilso\\Desktop\\scrna\\train.tsne.pdf")
DimPlot(object = train,reduction="tsne")
dev.off()
train <- RunUMAP(object = train, dims = 1:30)
pdf("C:\\Users\\wilso\\Desktop\\scrna\\train.umap.pdf")
DimPlot(object = train,reduction="umap")
dev.off()

#MARKER GENE heatmap
markers <- FindAllMarkers(train,logfc.threshold=0.5,test.use="wilcox",min.pct=0.25,only.pos=TRUE)
write.table(markers, file="C:\\Users\\wilso\\Desktop\\scrna\\train.cellmarker.xls",sep="\t",row.names=F,col.names=T,quote=F)
##this table contain all marker genes
#find most specific marker genes and save as marker_sub
markers_sub <- read.delim("C:/Users/wilso/Desktop/scrna/train.cellmarker.sub.xls")
top100 <- markers_sub %>% top_n(n = 100, wt=avg_log2FC)
pdf("C:\\Users\\wilso\\Desktop\\scrna\\train.marker_sub1.heatmap.pdf")
DoHeatmap(train, features = unique(top100$gene))
dev.off()

#violin plot
markers_sub1 <- read.delim("C:/Users/wilso/Desktop/scrna/train.cellmarker.sub1.xls")
top1 <- markers_sub1 %>% group_by(cluster) %>% top_n(n = 1, wt=avg_log2FC)
pdf("C:\\Users\\wilso\\Desktop\\scrna\\train.marker1.vlnplot.pdf")
VlnPlot(train, features = unique(top1$gene), pt.size=0)
dev.off()
ggsave("C:\\Users\\wilso\\Desktop\\scrna\\train.marker.vlnplot.png")
#
#
pdf("C:\\Users\\wilso\\Desktop\\scrna\\train.marker.featureplot.umap.pdf")
FeaturePlot(train, features = unique(top1$gene),reduction="umap")
dev.off()
ggsave("C:\\Users\\wilso\\Desktop\\scrna\\train.marker.featureplot.umap.png")
pdf("C:\\Users\\wilso\\Desktop\\scrna\\train.marker.dotplot.pdf")
DotPlot(object = train, features = unique(top1$gene)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
ggsave("marker.dotplot.png")
saveRDS(train, file = "C:\\Users\\wilso\\Desktop\\scrna\\train.rds")
#train  <- readRDS("C:\\Users\\wilso\\Desktop\\scrna\\train.rds")
#
#

