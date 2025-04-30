
library(GEOquery)
library(DESeq2)
library(gtools)
library(biobroom)
library(dplyr)
library(ggplot2)
library(pheatmap)


#### Loading datasets #### 

#load data
temp_mp <- read.table("Fcount_Vilcre_MP.txt", row.names =1 )
colnames(temp_mp) <- temp_mp[1,]
temp_mp <- temp_mp[2:nrow(temp_mp),]

#arrange dataset 
colnames(temp_mp)<- sub("\\-.*", "", colnames(temp_mp))
temp_mp <- temp_mp[,c("Chr", "Start", "End", "Strand", "Length",
                      "MP3", "MP4", "MP7", "MP8","MP10", "MP13a", "MP14", "MP16a")] 
# remove chr,start,end,strand,length
mp<- temp_mp[,6:ncol(temp_mp)]
mp_cond <- data.frame(colnames(mp),
                      mp = factor(c(rep("MP",2), rep("noMP", 2),
                                    rep("MP", 2), rep("noMP",2)), levels = c("noMP","MP")),
                      diet = factor(c(rep("LP",2), rep("LP", 2),
                                      rep("Ctrl", 2), rep("Ctrl",2)), levels = c("Ctrl","LP")),
                      conds = factor(c(rep("LP+MP",2), rep("LP", 2),
                                       rep("Ctrl+MP", 2), rep("Ctrl",2)), levels = c("Ctrl", "Ctrl+MP", "LP", "LP+MP")))


#convert numeric
mp1<- dplyr::mutate_all(mp, function(x) as.numeric(as.character(x)))


#### DESeq ####

model.matrix(~ conds, mp_cond)
dds_mp <- DESeqDataSetFromMatrix(countData = as.matrix(mp1), 
                                 colData = as.data.frame(mp_cond),
                                 design = ~ conds) 
dds_mp <- DESeq(dds_mp)
resultsNames(dds_mp)

#### Plot PCA ####

vds_mp <- vst(dds_mp, blind = F)
vsn::meanSdPlot(assay(vds_mp))
plotPCA(vds_mp, intgroup=c("conds"))

#### Volcano plot #### 

res.mp_filt<- results(dds_mp,alpha = 0.05, #lfcThreshold = log2(0), 
                      contrast = c("conds","Ctrl.MP", "Ctrl"))
summary(res.mp_filt)
DEG_mp_CtrlMPvCtrl<- res.mp_filt@rownames[res.mp_filt$padj < 0.05 & !is.na(res.mp_filt$padj)]

t.res.ctrl_mp <- tidy.DESeqResults(res.mp_filt)
t.res.ctrl_mp <- arrange(t.res.ctrl_mp, p.adjusted)
t.res.ctrl_mp <- inner_join (t.res.ctrl_mp, grcm38, by = c("gene" = "ensgene"))

ctrl_mp_volp <- data.frame(ens = t.res.ctrl_mp$gene,
                        gene = t.res.ctrl_mp$symbol,
                        log2FC = t.res.ctrl_mp$estimate,
                        p.adj = t.res.ctrl_mp$p.adjusted)
ctrl_mp_volp$deg <- "ns"

ctrl_mp_volp$deg[ctrl_mp_volp$log2FC < -log2(1.5)] <- "down"
ctrl_mp_volp$deg[ctrl_mp_volp$log2FC > log2(1.5)] <- "up"
ctrl_mp_volp$deg[ctrl_mp_volp$p.adj > 0.05] <- "ns"

table(ctrl_mp_volp$deg[!is.na(ctrl_mp_volp$p.adj)])


# For LP

res.mp_filt<- results(dds_mp,alpha = 0.05, #lfcThreshold = log2(0), 
                      contrast = c("conds","LP.MP", "LP"))
summary(res.mp_filt)
DEG_mp_LPMPvLP<- res.mp_filt@rownames[res.mp_filt$padj < 0.05 & !is.na(res.mp_filt$padj)]

t.res.LP_mp <- tidy.DESeqResults(res.mp_filt)
t.res.LP_mp <- arrange(t.res.LP_mp, p.adjusted)
t.res.LP_mp <- inner_join (t.res.LP_mp, grcm38, by = c("gene" = "ensgene"))

LP_mp_volp <- data.frame(ens = t.res.LP_mp$gene,
                           gene = t.res.LP_mp$symbol,
                           log2FC = t.res.LP_mp$estimate,
                           p.adj = t.res.LP_mp$p.adjusted)
LP_mp_volp$deg <- "ns"
LP_mp_volp$deg[LP_mp_volp$log2FC < -log2(1.5)] <- "down"
LP_mp_volp$deg[LP_mp_volp$log2FC > log2(1.5)] <- "up"
LP_mp_volp$deg[LP_mp_volp$p.adj > 0.05] <- "ns"


# plot adding up all layers we have seen so far
ggplot(data=ctrl_mp_volp[,2:5],aes(x=log2FC, y=-log10(p.adj), col=deg)) +
  geom_point() +
  theme_minimal() +
  #geom_text_repel() +
  scale_color_manual(values=c("blue", "grey", "tomato")) +
  geom_vline(xintercept=c(log2(2/3), log2(1.5)), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  #geom_text_repel() +
  #theme_minimal() #+
  xlim(-10,18)+
  ylim(0,18) +
  ggtitle("Ctrl+MP vs Ctrl")


ggplot(data=LP_mp_volp[,2:5],aes(x=log2FC, y=-log10(p.adj), col=deg)) +
  geom_point() +
  theme_minimal() +
  #geom_text_repel() +
  scale_color_manual(values=c("blue", "grey", "tomato")) +
  geom_vline(xintercept=c(log2(2/3), log2(1.5)), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  #geom_text_repel() +
  #theme_minimal() #+
  xlim(-10,18)+
  ylim(0,18) +
  ggtitle("LP+MP vs LP")


#### DESeq comparisons ####

res.mp <-results(dds_mp , contrast = c("conds","LP", "Ctrl")) 
t.res.mp <- tidy.DESeqResults(res.mp)
t.res.mp <- arrange(t.res.mp, p.adjusted)


res.mp_filt<- results(dds_mp,alpha = 0.05, contrast = c("conds","LP", "Ctrl"))
DEG_mp_LPvCtrl<- res.mp_filt@rownames[res.mp_filt$padj < 0.05 & !is.na(res.mp_filt$padj)]

res.mp_filt<- results(dds_mp,alpha = 0.05, contrast = c("conds","Ctrl+MP", "Ctrl"))
DEG_mp_CtrlMPvCtrl<- res.mp_filt@rownames[res.mp_filt$padj < 0.05 & !is.na(res.mp_filt$padj)]

res.mp_filt<- results(dds_mp,alpha = 0.05, contrast = c("conds","LP+MP", "LP"))
DEG_mp_LPMPvLP<- res.mp_filt@rownames[res.mp_filt$padj < 0.05 & !is.na(res.mp_filt$padj)]

res.mp_filt<- results(dds_mp,alpha = 0.05, contrast = c("conds","LP+MP", "Ctrl+MP"))
DEG_mp_LPMPvCtrlMP<- res.mp_filt@rownames[res.mp_filt$padj < 0.05 & !is.na(res.mp_filt$padj)]


# all mp
model.matrix(~ mp, mp_cond)
dds_mp <- DESeqDataSetFromMatrix(countData = as.matrix(mp1), 
                                 colData = as.data.frame(mp_cond),
                                 design = ~mp) 
dds_mp <- DESeq(dds_mp)
resultsNames(dds_mp)

res.mp <-results(dds_mp , contrast = c("mp","MP", "noMP")) 
t.res.mp <- tidy.DESeqResults(res.mp)
t.res.mp <- arrange(t.res.mp, p.adjusted)

res.mp_filt<- results(dds_mp,alpha = 0.05, contrast = c("mp","MP", "noMP"))
DEG_mp_all_MPvnoMP<- res.mp_filt@rownames[res.mp_filt$padj < 0.05 & !is.na(res.mp_filt$padj)]


# all diet
model.matrix(~ diet, mp_cond)
dds_mp <- DESeqDataSetFromMatrix(countData = as.matrix(mp1), 
                                 colData = as.data.frame(mp_cond),
                                 design = ~diet) 
dds_mp <- DESeq(dds_mp)

res.mp <-results(dds_mp , contrast = c("diet","LP", "Ctrl")) 
t.res.mp <- tidy.DESeqResults(res.mp)
t.res.mp <- arrange(t.res.mp, p.adjusted)

res.mp_filt<- results(dds_mp,alpha = 0.05, contrast = c("diet","LP", "Ctrl"))
DEG_mp_all_LPvCtrl<- res.mp_filt@rownames[res.mp_filt$padj < 0.05 & !is.na(res.mp_filt$padj)]


##### Venn diagram #####

ctrl_up <- ctrl_mp_volp[!is.na(ctrl_mp_volp$p.adj),]
ctrl_up <- ctrl_up$gene[which(ctrl_up$deg == "up")]

ctrl_down <- ctrl_mp_volp[!is.na(ctrl_mp_volp$p.adj),]
ctrl_down <- ctrl_down$gene[which(ctrl_down$deg == "down")]


LP_up <- LP_mp_volp[!is.na(LP_mp_volp$p.adj),]
LP_up <- LP_up$gene[which(LP_up$deg == "up")]

LP_down <- LP_mp_volp[!is.na(LP_mp_volp$p.adj),]
LP_down <- LP_down$gene[which(LP_down$deg == "down")]


AllMP_up <- AllMP_mp_volp[!is.na(AllMP_mp_volp$p.adj),]
AllMP_up <- AllMP_up$gene[which(AllMP_up$deg == "up")]

AllMP_down <- AllMP_mp_volp[!is.na(AllMP_mp_volp$p.adj),]
AllMP_down <- AllMP_down$gene[which(AllMP_down$deg == "down")]

library(VennDiagram)
venn.diagram(
  x= list(ctrl_up, LP_up, AllMP_up),
  category.names = c("ctrl_up", "LP_up", "AllMP_up"),
  filename = "up_deg_mp_venn.png",
  output = TRUE
)

venn.diagram(
  x= list(ctrl_down, LP_down, AllMP_down),
  category.names = c("ctrl_down", "LP_down", "AllMP_down"),
  filename = "down_deg_mp_venn.png",
  output = TRUE
)



####DEG_ALL#### 

DEG_all <- union(DEG_mp_LPvCtrl, DEG_mp_CtrlMPvCtrl) 
DEG_all <- union(DEG_all, DEG_mp_LPMPvLP) 
DEG_all <- union(DEG_all, DEG_mp_LPMPvCtrlMP) 
DEG_all <- union(DEG_all, DEG_mp_all_MPvnoMP)
DEG_all <- union(DEG_all, DEG_mp_all_LPvCtrl)

length(DEG_all)

pheatmap(assay(vds_mp)[DEG_all, c("MP14", "MP16a", "MP10", "MP13a", "MP7", "MP8", "MP3","MP4") ], 
         scale = "row", cluster_rows=T, show_rownames=FALSE, 
         cluster_cols=F, show_colnames = T,
         annotation_col = df_col, 
         annotation_colors = annot_color ) #confirm 


#### GO ####

library(clusterProfiler)
library(AnnotationHub)
library(dplyr)
library(annotables)
library(org.Mm.eg.db)
library(ggplot2)
library(DESeq2)
library(pheatmap)

set.seed(8888)


cl_mp<- pheatmap(assay(vds_mp)[DEG_all, c("MP14", "MP16a", "MP10", "MP13a", "MP7", "MP8", "MP3","MP4") ], 
                 scale = "row", cluster_rows=T, show_rownames=FALSE, 
                 cluster_cols=F, show_colnames = T, cutree_rows = 4,
                 annotation_col = df_col, 
                 annotation_colors = annot_color ) 

pheatmap(assay(vds_mp)[names(cl_mp)[cl_mp %in% c(3:10)], c("MP14", "MP16a", "MP10", "MP13a", "MP7", "MP8", "MP3","MP4") ], 
         scale = "row", cluster_rows=T, show_rownames=FALSE, 
         cluster_cols=F, show_colnames = T, cutree_rows = 6,
         annotation_col = df_col, 
         annotation_colors = annot_color,
         main = paste0("MP_cluster_1"))


cl_mp <-pheatmap(assay(vds_mp)[names(cl_mp)[cl_mp %in% c(3:10)], c("MP14", "MP16a", "MP10", "MP13a", "MP7", "MP8", "MP3","MP4") ], 
                 scale = "row", cluster_rows=T, show_rownames=FALSE, 
                 cluster_cols=F, show_colnames = T, cutree_rows = 6,
                 annotation_col = df_col, 
                 annotation_colors = annot_color,
                 main = paste0("MP_cluster_1"))


cl_mp <- cutree(cl_mp$tree_row, 6)
table(cl_mp)


##### only Ctrl+MP vs Ctrl #####


cl_mp <- pheatmap(assay(vds_mp)[DEG_all, c("MP14", "MP16a", "MP10", "MP13a", "MP7", "MP8", "MP3","MP4") ],
                  scale = "row", cluster_rows=T, show_rownames=F,
                  cluster_cols=F, show_colnames = T, cutree_rows = 10,
                  annotation_col = df_col, 
                  annotation_colors = annot_color ) 


cl_mp <- cutree(cl_mp$tree_row, 10)
table(cl_mp)


pheatmap(assay(vds_mp)[DEG_all, c("MP14", "MP16a", "MP10", "MP13a", "MP7", "MP8", "MP3","MP4") ],
         scale = "row", cluster_rows=T, show_rownames=F, 
         cluster_cols=F, show_colnames = T, cutree_rows = 10,
         annotation_col = df_col, 
         annotation_colors = annot_color )


pheatmap(assay(vds_mp)[names(cl_mp)[cl_mp == 9], c("MP14", "MP16a", "MP10", "MP13a", "MP7", "MP8", "MP3","MP4") ], 
         scale = "row", cluster_rows=F, show_rownames=FALSE, 
         cluster_cols=F, show_colnames = T, 
         annotation_col = df_col, 
         annotation_colors = annot_color,
         main = paste0("MP_cluster_1"))

pheatmap(assay(vds_mp)[names(cl_mp)[cl_mp %in% c(3:10)], c("MP14", "MP16a", "MP10", "MP13a", "MP7", "MP8", "MP3","MP4") ], 
         scale = "row", cluster_rows=T, show_rownames=FALSE, 
         cluster_cols=F, show_colnames = T, cutree_rows = 6,
         annotation_col = df_col, 
         annotation_colors = annot_color,
         main = paste0("MP_cluster_1"))


cl_mp <-pheatmap(assay(vds_mp)[names(cl_mp)[cl_mp %in% c(3:10)], c("MP14", "MP16a", "MP10", "MP13a") ], 
                 scale = "row", cluster_rows=T, show_rownames=FALSE, 
                 cluster_cols=F, show_colnames = T, cutree_rows = 3,
                 annotation_col = df_col, 
                 annotation_colors = annot_color,
                 main = paste0("MP_cluster_1"))


cl_mp <- cutree(cl_mp$tree_row, 6)
table(cl_mp)

pheatmap(assay(vds_mp)[names(cl_mp)[cl_mp == 1], c("MP14", "MP16a", "MP10", "MP13a", "MP7", "MP8", "MP3","MP4") ], 
         scale = "row", cluster_rows=T, show_rownames=FALSE, 
         cluster_cols=F, show_colnames = T, 
         annotation_col = df_col, 
         annotation_colors = annot_color,
         main = paste0("MP_cluster_1"))








##### All together #####


set.seed(8888)

cl_mp <-pheatmap(assay(vds_mp)[names(cl_mp)[cl_mp %in% c(3:10)], c("MP14", "MP16a", "MP10", "MP13a", "MP7", "MP8", "MP3","MP4") ], 
                 scale = "row", cluster_rows=T, show_rownames=FALSE, 
                 cluster_cols=F, show_colnames = T, cutree_rows = 6,
                 annotation_col = df_col, 
                 annotation_colors = annot_color,
                 main = paste0("MP_cluster_1"))


cl_mp <- cutree(cl_mp$tree_row, 6)
table(cl_mp)

pheatmap(assay(vds_mp)[names(cl_mp)[cl_mp == 1], c("MP14", "MP16a", "MP10", "MP13a", "MP7", "MP8", "MP3","MP4") ], 
         scale = "row", cluster_rows=T, show_rownames=FALSE, 
         cluster_cols=F, show_colnames = T, 
         annotation_col = df_col, 
         annotation_colors = annot_color,
         main = paste0("MP_cluster_1"))



## Cluster 1

HM_1<- pheatmap(assay(vds_mp)[names(cl_mp)[cl_mp == 1], c("MP14", "MP16a", "MP10", "MP13a", "MP7", "MP8", "MP3","MP4") ], 
                scale = "row", cluster_rows=F, show_rownames=FALSE, 
                cluster_cols=F, show_colnames = T, 
                annotation_col = df_col, 
                annotation_colors = annot_color,
                main = paste0("MP_cluster_1"))


GO_1 <-enrichGO(names(cl_mp)[cl_mp == 1],
                OrgDb         = org.Mm.eg.db,  
                ont           = "BP",
                pAdjustMethod = "BH",
                keyType       = 'ENSEMBL',
                pvalueCutoff  = 0.3,
                qvalueCutoff  = 0.3,
                readable      = TRUE)


write.csv(GO_1@result, paste0("MP_cluster_1.csv"))

dotplot(GO_1, showCategory = 20,
        title = "DEG_ALL MP GO cluster_1")



## Cluster 2

HM_2<- pheatmap(assay(vds_mp)[names(cl_mp)[cl_mp == 2], c("MP14", "MP16a", "MP10", "MP13a", "MP7", "MP8", "MP3","MP4") ], 
                scale = "row", cluster_rows=F, show_rownames=FALSE, 
                cluster_cols=F, show_colnames = T, 
                annotation_col = df_col, 
                annotation_colors = annot_color,
                main = paste0("MP_cluster_2"))


GO_2 <-enrichGO(names(cl_mp)[cl_mp == 2],
                OrgDb         = org.Mm.eg.db,  
                ont           = "BP",
                pAdjustMethod = "BH",
                keyType       = 'ENSEMBL',
                pvalueCutoff  = 0.3,
                qvalueCutoff  = 0.3,
                readable      = TRUE)


write.csv(GO_2@result, paste0("MP_cluster_2.csv"))

dotplot(GO_2, showCategory = 20,
        title = "DEG_ALL MP GO cluster_2")


## cluster 3

HM_3<- pheatmap(assay(vds_mp)[names(cl_mp)[cl_mp == 3], c("MP14", "MP16a", "MP10", "MP13a", "MP7", "MP8", "MP3","MP4") ], 
                scale = "row", cluster_rows=F, show_rownames=FALSE, 
                cluster_cols=F, show_colnames = T, 
                annotation_col = df_col, 
                annotation_colors = annot_color,
                main = paste0("MP_cluster_3"))


GO_3 <-enrichGO(names(cl_mp)[cl_mp == 3],
                OrgDb         = org.Mm.eg.db,  
                ont           = "BP",
                pAdjustMethod = "BH",
                keyType       = 'ENSEMBL',
                pvalueCutoff  = 0.3,
                qvalueCutoff  = 0.3,
                readable      = TRUE)


write.csv(GO_3@result, paste0("MP_cluster_3.csv"))

dotplot(GO_3, showCategory = 20,
        title = "DEG_ALL MP GO cluster_3")

## cluster 4

HM_4<- pheatmap(assay(vds_mp)[names(cl_mp)[cl_mp == 4], c("MP14", "MP16a", "MP10", "MP13a", "MP7", "MP8", "MP3","MP4") ], 
                scale = "row", cluster_rows=F, show_rownames=FALSE, 
                cluster_cols=F, show_colnames = T, 
                annotation_col = df_col, 
                annotation_colors = annot_color,
                main = paste0("MP_cluster_4"))


GO_4 <-enrichGO(names(cl_mp)[cl_mp == 4],
                OrgDb         = org.Mm.eg.db,  
                ont           = "BP",
                pAdjustMethod = "BH",
                keyType       = 'ENSEMBL',
                pvalueCutoff  = 0.3,
                qvalueCutoff  = 0.3,
                readable      = TRUE)


write.csv(GO_4@result, paste0("MP_cluster_4.csv"))

dotplot(GO_4, showCategory = 20,
        title = "DEG_ALL MP GO cluster_4")

## cluster 5

HM_5<- pheatmap(assay(vds_mp)[names(cl_mp)[cl_mp == 5], c("MP14", "MP16a", "MP10", "MP13a", "MP7", "MP8", "MP3","MP4") ], 
                scale = "row", cluster_rows=F, show_rownames=FALSE, 
                cluster_cols=F, show_colnames = T, 
                annotation_col = df_col, 
                annotation_colors = annot_color,
                main = paste0("MP_cluster_5"))


GO_5 <-enrichGO(names(cl_mp)[cl_mp == 5],
                OrgDb         = org.Mm.eg.db,  
                ont           = "BP",
                pAdjustMethod = "BH",
                keyType       = 'ENSEMBL',
                pvalueCutoff  = 0.3,
                qvalueCutoff  = 0.3,
                readable      = TRUE)


write.csv(GO_5@result, paste0("MP_cluster_5.csv"))

dotplot(GO_5, showCategory = 20,
        title = "DEG_ALL MP GO cluster_5")

## cluster 6

HM_6<- pheatmap(assay(vds_mp)[names(cl_mp)[cl_mp == 6], c("MP14", "MP16a", "MP10", "MP13a", "MP7", "MP8", "MP3","MP4") ], 
                scale = "row", cluster_rows=F, show_rownames=FALSE, 
                cluster_cols=F, show_colnames = T, 
                annotation_col = df_col, 
                annotation_colors = annot_color,
                main = paste0("MP_cluster_6"))

GO_6 <-enrichGO(names(cl_mp)[cl_mp == 6],
                OrgDb         = org.Mm.eg.db,  
                ont           = "BP",
                pAdjustMethod = "BH",
                keyType       = 'ENSEMBL',
                pvalueCutoff  = 0.3,
                qvalueCutoff  = 0.3,
                readable      = TRUE)

write.csv(GO_6@result, paste0("MP_cluster_6.csv"))

dotplot(GO_6, showCategory = 20,
        title = "DEG_ALL MP GO cluster_6")