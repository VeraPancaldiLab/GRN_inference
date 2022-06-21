#load files
p1_rep1 <- read.expr.matrix("Data/Gene_counts_Patient1_rep1_TPM.txt",form = "rows.are.genes")

p1_rep2 <- read.expr.matrix("Data/Gene_counts_Patient1_rep2_TPM.txt", form = "rows.are.genes")

p2_rep1 <- read.expr.matrix("Data/Gene_counts_Patient2_rep1_TPM.txt", form = "rows.are.genes")

p2_rep2 <- read.expr.matrix("Data/Gene_counts_Patient2_rep2_TPM.txt", form = "rows.are.genes")



for (i in 1:ncol(p2_rep1)){
  print(table(is.infinite(p2_rep1[,i])))
}

#corr
#patient1
library(corrplot)
corrplot(cor(p1_rep1[,2:ncol(p1_rep1)]),method="color",addCoef.col=T,cl.lim = c(0.94, 1),
         is.corr = F)

corrplot(cor(p1_rep2[,2:ncol(p1_rep2)]),method="color",addCoef.col=T,cl.lim = c(0.94, 1),
         is.corr = F)

#patient2
corrplot(cor(p2_rep1[,2:ncol(p2_rep1)]),method="color",addCoef.col=T,cl.lim = c(0.94, 1),
         is.corr = F)

corrplot(cor(p2_rep2[,2:ncol(p2_rep2)]),method="color",addCoef.col=T,cl.lim = c(0.94, 1),
         is.corr = F)



#PCA
library(ggplot2)
library(RColorBrewer)

#Determine batch effect between 2 replicates of patient 1
PCA <- prcomp(t(cbind(p1_rep1,p1_rep2)), scale = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG_p1 <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])
ggplot(dataGG_p1 %>% rownames_to_column("Patient"), aes(PC1, PC2,label=Patient)) +
  geom_point() +
  geom_text(aes(label=Patient),hjust=-.2, vjust=.4) 

#Determine batch effect between 2 replicates of patient 2
PCA <- prcomp(t(cbind(p2_rep1,p2_rep2)), scale = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG_p2 <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])
ggplot(dataGG_p2 %>% rownames_to_column("Patient"), aes(PC1, PC2,label=Patient)) +
  geom_point() +
  geom_text(aes(label=Patient),hjust=-.2, vjust=.4) 
#No batch effect, these results are expected because we are in the case of 2 technical replicates

#Determine batch effect between the 2 patients
all_data = lapply(list(p1_rep1,p1_rep2,p2_rep1,p2_rep2), cbind) %>% as.data.frame()
group = str_split_fixed(colnames(all_data),pattern = "_",n = 3)[,2]

PCA <- prcomp(t(all_data), scale = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])
ggplot(dataGG %>% rownames_to_column("Patient") %>% mutate(P=str_split_fixed(colnames(all_data),pattern = "_",n = 3)[,2],
                                                           J=str_split_fixed(colnames(all_data),pattern = "_",n = 3)[,1]), aes(PC1, PC2,label=Patient,col=P)) +
  geom_point() +
  geom_text(aes(label=Patient),hjust=-.2, vjust=.4) +
  ggtitle("PCA - Patient 1 and Patient 2") + 
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) + 
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
#Batch effect is present between patients

#Batch effect correction
batch_corrected_data = sva::ComBat_seq(counts = as.matrix(all_data),batch = group)

PCA <- prcomp(t(batch_corrected_data), scale = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])
ggplot(dataGG %>% rownames_to_column("Patient") %>% mutate(P=str_split_fixed(colnames(batch_corrected_data),pattern = "_",n = 3)[,2],
                                                           J=str_split_fixed(colnames(batch_corrected_data),pattern = "_",n = 3)[,1]), aes(PC1, PC2,label=Patient,col=P)) +
  geom_point() +
  geom_text(aes(label=Patient),hjust=-.2, vjust=.4) +
  ggtitle("PCA - Patient 1 and Patient 2 - Batch effect correction") + 
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) + 
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
#No more batch effect between patients

write.table(batch_corrected_data, file="Results/batch_corrected_data.txt", sep="\t", quote=F, col.names=NA)

#########
# TCseq #
#########
library(TCseq)
library(EnsDb.Hsapiens.v75)

###### 
#patient 1

##### create TCA object
### design table
design=data.frame(sampleid=c(colnames(p1_rep1)[2:12],colnames(p1_rep2)[2:12]),
                  timepoint=rep(c("1d","2d","3d","4d","5d","6d","7d","8d","9d","11d","13d"),2),
                  group=rep(1:11,2))

### Count matrix
cm_p1=merge(p1_rep1,p1_rep2,by="Gene")
rownames(cm_p1)=cm_p1$Gene;cm_p1=cm_p1[,-1]
cm_p1=as.matrix(cm_p1)
cm_p1=round(cm_p1)#if not, values are truncated

### genomic intervals
# Parse Ensembl genes
Ens=as.data.frame(genes(EnsDb.Hsapiens.v75)) #hg19
Ens=Ens[which(Ens$gene_biotype=="protein_coding"),]
Ens$seqnames=paste0("chr",Ens$seqnames)
rownames(Ens)=NULL
Ens=Ens[,c(1,2,3,7)]
colnames(Ens)=c("chr","start","end","id")
chr=c(paste0("chr",1:22),"chrX")
Ens=Ens[which(Ens$chr %in% chr),]

# Select genes
genes=data.frame()
for (i in rownames(cm_p1)){
  genes=rbind(genes,Ens[which(Ens$id==i),])
}

# keep first of duplicated genes
remove=vector()
for (i in unique(genes$id[which(isUnique(genes$id)==FALSE)])){
  remove=c(remove,rownames(genes[which(genes$id==i),])[2:nrow(genes[which(genes$id==i),])])
}
remove=as.numeric(remove)
genes=genes[which(!(rownames(genes) %in% remove)),]

cm_p1=cm_p1[which(rownames(cm_p1) %in% genes$id),]

# Verif
table(genes$id[which(isUnique(genes$id)==FALSE)]) #should be empty
dim(genes)
dim(cm_p1)
table(rownames(cm_p1)==genes$id) #should be only TRUE

### tca_p1 object
tca_p1 <- TCA(design = design, genomicFeature = genes, counts = cm_p1)

##### Diff Analysis
#tca_p1 <- DBanalysis(tca_p1)
tca_p1 <- DBanalysis(tca_p1, filter.type = "cpm", filter.value = 5, samplePassfilter = 2) #filter low read counts (at least 5 in 2 samples)

#Extract results
DBres <- DBresult(tca_p1, group1 = "1d", group2 = c("13d"))
res=as.data.frame(DBres$`13dvs1d`)

write.table(res[which(res$logFC>=2 & res$PValue<=0.05),],"/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Signif/TCseq_13dvs1d_UP.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(res[which(res$logFC<=-2 & res$PValue<=0.05),],"/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Signif/TCseq_13dvs1d_DOWN.txt",col.names=T,row.names=F,quote=F,sep="\t")

signif_2dvs1d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_1d/Patient1/TCseq_2dvs1d_DOWN.txt",header = T)
signif_3dvs1d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_1d/Patient1/TCseq_3dvs1d_DOWN.txt",header = T)
signif_4dvs1d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_1d/Patient1/TCseq_4dvs1d_DOWN.txt",header = T)
signif_5dvs1d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_1d/Patient1/TCseq_5dvs1d_DOWN.txt",header = T)
signif_6dvs1d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_1d/Patient1/TCseq_6dvs1d_DOWN.txt",header = T)
signif_7dvs1d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_1d/Patient1/TCseq_7dvs1d_DOWN.txt",header = T)
signif_8dvs1d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_1d/Patient1/TCseq_8dvs1d_DOWN.txt",header = T)
signif_9dvs1d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_1d/Patient1/TCseq_9dvs1d_DOWN.txt",header = T)
signif_11dvs1d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_1d/Patient1/TCseq_11dvs1d_DOWN.txt",header = T)
signif_13dvs1d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_1d/Patient1/TCseq_13dvs1d_DOWN.txt",header = T)

vec_UP=c(nrow(signif_2dvs1d_UP),
  nrow(signif_3dvs1d_UP),
  nrow(signif_4dvs1d_UP),
  nrow(signif_5dvs1d_UP),
  nrow(signif_6dvs1d_UP),
  nrow(signif_7dvs1d_UP),
  nrow(signif_8dvs1d_UP),
  nrow(signif_9dvs1d_UP),
  nrow(signif_11dvs1d_UP),
  nrow(signif_13dvs1d_UP))

vec_DOWN=c(nrow(signif_2dvs1d_DOWN),
      nrow(signif_3dvs1d_DOWN),
      nrow(signif_4dvs1d_DOWN),
      nrow(signif_5dvs1d_DOWN),
      nrow(signif_6dvs1d_DOWN),
      nrow(signif_7dvs1d_DOWN),
      nrow(signif_8dvs1d_DOWN),
      nrow(signif_9dvs1d_DOWN),
      nrow(signif_11dvs1d_DOWN),
      nrow(signif_13dvs1d_DOWN))

#UP
b=barplot(vec_UP,col="#77dd77",names=c("2d","3d","4d","5d","6d","7d","8d","9d","11d","13d"),ylim=c(0,250),main="UP deregulated genes - ref 1d")
text(b, vec_UP, labels=as.character(vec_UP), cex= .9,pos=3)

#DOWN
b=barplot(-vec_DOWN,col="#ff6961",names=c("2d","3d","4d","5d","6d","7d","8d","9d","11d","13d"),ylim=c(-250,0),main="DOWN deregulated genes - ref 1d")
text(b, -vec_DOWN, labels=as.character(vec_DOWN), cex= .9,pos=1)


DBres.sig <- DBresult(tca_p1, group1 = "1d", group2 = c("2d"), top.sig = TRUE) #only logFC > 2 and pval < 0.05



### Biomart ###
library("biomaRt")
#ensembl=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

UP=signif_13dvs1d_UP$id
DOWN=signif_13dvs1d_DOWN$id

#UP
genes.with.id=getBM(attributes=c("ensembl_gene_id","hgnc_symbol","ensembl_gene_id_version", "entrezgene_id"),values=UP, mart= ensembl)
id=genes.with.id[which(genes.with.id$hgnc_symbol %in% UP),]
myentrez_13dvs1d_UP=unique(id$entrezgene_id)

#DOWN
genes.with.id=getBM(attributes=c("ensembl_gene_id","hgnc_symbol","ensembl_gene_id_version", "entrezgene_id"),values=DOWN, mart= ensembl)
id=genes.with.id[which(genes.with.id$hgnc_symbol %in% DOWN),]
myentrez_13dvs1d_DOWN=unique(id$entrezgene_id)



myentrez_list_UP=list(UP_2dvs1d=myentrez_2dvs1d_UP,UP_3dvs1d=myentrez_3dvs1d_UP,UP_4dvs1d=myentrez_4dvs1d_UP,
                      UP_5dvs1d=myentrez_5dvs1d_UP,UP_6dvs1d=myentrez_6dvs1d_UP,UP_7dvs1d=myentrez_7dvs1d_UP,
                      UP_8dvs1d=myentrez_8dvs1d_UP,UP_9dvs1d=myentrez_9dvs1d_UP,UP_11dvs1d=myentrez_11dvs1d_UP,
                      UP_13dvs1d=myentrez_13dvs1d_UP)


myentrez_list_DOWN=list(DOWN_2dvs1d=myentrez_2dvs1d_DOWN,DOWN_3dvs1d=myentrez_3dvs1d_DOWN,DOWN_4dvs1d=myentrez_4dvs1d_DOWN,
                      DOWN_5dvs1d=myentrez_5dvs1d_DOWN,DOWN_6dvs1d=myentrez_6dvs1d_DOWN,DOWN_7dvs1d=myentrez_7dvs1d_DOWN,
                      DOWN_8dvs1d=myentrez_8dvs1d_DOWN,DOWN_9dvs1d=myentrez_9dvs1d_DOWN,DOWN_11dvs1d=myentrez_11dvs1d_DOWN,
                      DOWN_13dvs1d=myentrez_13dvs1d_DOWN)




myentrez_list_ALL=list(UP_2dvs1d=c(myentrez_2dvs1d_UP,myentrez_2dvs1d_DOWN),UP_3dvs1d=c(myentrez_3dvs1d_UP,myentrez_3dvs1d_DOWN),UP_4dvs1d=c(myentrez_4dvs1d_UP,myentrez_4dvs1d_DOWN),
                      UP_5dvs1d=c(myentrez_5dvs1d_UP,myentrez_5dvs1d_DOWN),UP_6dvs1d=c(myentrez_6dvs1d_UP,myentrez_6dvs1d_DOWN),UP_7dvs1d=c(myentrez_7dvs1d_UP,myentrez_7dvs1d_DOWN),
                      UP_8dvs1d=c(myentrez_8dvs1d_UP,myentrez_8dvs1d_DOWN),UP_9dvs1d=c(myentrez_9dvs1d_UP,myentrez_9dvs1d_DOWN),UP_11dvs1d=c(myentrez_11dvs1d_UP,myentrez_11dvs1d_DOWN),
                      UP_13dvs1d=c(myentrez_13dvs1d_UP,myentrez_13dvs1d_DOWN))


#GO
library(limma)
library(clusterProfiler)

#goana=goana(de=myentrez,geneid=myentrez,FDR=0.95,species="Hs")
#View(goana)

#GO Plot
compGO <- compareCluster(geneCluster     = myentrez_list_UP,
                         fun           = "enrichGO",
                         OrgDb         = 'org.Hs.eg.db',
                         ont           = 'BP',
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")
dotplot(compGO, showCategory = 15, title = "Gene Ontology Enrichment Analysis")


#Temporal pattern analysis
#tca_p1 <- timecourseTable(tca_p1, value = "FC", norm.method = "tpm", filter = TRUE)
tca_p1 <- timecourseTable(tca_p1, value = "expression", norm.method = "cpm", filter = TRUE)
t <- tcTable(tca_p1)

#clustering
tca_p1 <- timeclust(tca_p1, algo = "cm", k = 6, standardize = TRUE)

p <- timeclustplot(tca_p1, value = "z-score(PRKM)", cols = 3)

#plot cluster 1:
print(p[[1]])

#see which cluster each gene belongs to
tca_p1@clusterRes@cluster

table(tca_p1@clusterRes@cluster)

write.table(tca_p1@clusterRes@cluster,"/Users/fla./Desktop/CLL_RNA-Seq/TCseq/TCA_cmeans_cluster_patient1.txt",row.names=T,col.names = F,quote=F,sep="\t")


######
#patient 2
##### create TCA object
### design table
design=data.frame(sampleid=c(colnames(p2_rep1)[2:14],colnames(p2_rep2)[2:14]),
                  timepoint=rep(c("1d","2d","3d","4d","5d","6d","7d","8d","9d","10d","11d","12d","13d"),2),
                  group=rep(1:13,2))

### Count matrix
cm_p2=merge(p2_rep1,p2_rep2,by="Gene")
rownames(cm_p2)=cm_p2$Gene;cm_p2=cm_p2[,-1]
cm_p2=as.matrix(cm_p2)
cm_p2=round(cm_p2)#if not, values are truncated

### genomic intervals
# Parse Ensembl genes
Ens=as.data.frame(genes(EnsDb.Hsapiens.v75)) #hg19
Ens=Ens[which(Ens$gene_biotype=="protein_coding"),]
Ens$seqnames=paste0("chr",Ens$seqnames)
rownames(Ens)=NULL
Ens=Ens[,c(1,2,3,7)]
colnames(Ens)=c("chr","start","end","id")
chr=c(paste0("chr",1:22),"chrX")
Ens=Ens[which(Ens$chr %in% chr),]

# Select genes
genes=data.frame()
for (i in rownames(cm_p2)){
  genes=rbind(genes,Ens[which(Ens$id==i),])
}

# keep first of duplicated genes
remove=vector()
for (i in unique(genes$id[which(isUnique(genes$id)==FALSE)])){
  remove=c(remove,rownames(genes[which(genes$id==i),])[2:nrow(genes[which(genes$id==i),])])
}
remove=as.numeric(remove)
genes=genes[which(!(rownames(genes) %in% remove)),]

cm_p2=cm_p2[which(rownames(cm_p2) %in% genes$id),]

# Verif
table(genes$id[which(isUnique(genes$id)==FALSE)]) #should be empty
dim(genes)
dim(cm_p2)
table(rownames(cm_p2)==genes$id) #should be only TRUE

### TCA object
tca_p2 <- TCA(design = design, genomicFeature = genes, counts = cm_p2)

##### Diff Analysis
#tca_p1 <- DBanalysis(tca_p1)
tca_p2 <- DBanalysis(tca_p2, filter.type = "cpm", filter.value = 5, samplePassfilter = 2) #filter low read counts (at least 5 in 2 samples)

#Extract results
DBres <- DBresult(tca_p2, group1 = "12d", group2 = c("13d"))
res=as.data.frame(DBres$`13dvs12d`)

write.table(res[which(res$logFC>=2 & res$PValue<=0.05),],"/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_n-1/Patient2/TCseq_p2_13dvs12d_UP.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(res[which(res$logFC<=-2 & res$PValue<=0.05),],"/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_n-1/Patient2/TCseq_p2_13dvs12d_DOWN.txt",col.names=T,row.names=F,quote=F,sep="\t")


# vs 1st day
######

signif_2dvs1d_UP=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Signif/Patient2/TCseq_2dvs1d_UP.txt",header = T)
signif_3dvs1d_UP=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Signif/Patient2/TCseq_3dvs1d_UP.txt",header = T)
signif_4dvs1d_UP=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Signif/Patient2/TCseq_4dvs1d_UP.txt",header = T)
signif_5dvs1d_UP=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Signif/Patient2/TCseq_5dvs1d_UP.txt",header = T)
signif_6dvs1d_UP=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Signif/Patient2/TCseq_6dvs1d_UP.txt",header = T)
signif_7dvs1d_UP=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Signif/Patient2/TCseq_7dvs1d_UP.txt",header = T)
signif_8dvs1d_UP=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Signif/Patient2/TCseq_8dvs1d_UP.txt",header = T)
signif_9dvs1d_UP=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Signif/Patient2/TCseq_9dvs1d_UP.txt",header = T)
signif_10dvs1d_UP=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Signif/Patient2/TCseq_10dvs1d_UP.txt",header = T)
signif_11dvs1d_UP=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Signif/Patient2/TCseq_11dvs1d_UP.txt",header = T)
signif_12dvs1d_UP=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Signif/Patient2/TCseq_12dvs1d_UP.txt",header = T)
signif_13dvs1d_UP=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Signif/Patient2/TCseq_13dvs1d_UP.txt",header = T)

vec_UP=c(nrow(signif_2dvs1d_UP),
         nrow(signif_3dvs1d_UP),
         nrow(signif_4dvs1d_UP),
         nrow(signif_5dvs1d_UP),
         nrow(signif_6dvs1d_UP),
         nrow(signif_7dvs1d_UP),
         nrow(signif_8dvs1d_UP),
         nrow(signif_9dvs1d_UP),
         nrow(signif_10dvs1d_UP),
         nrow(signif_11dvs1d_UP),
         nrow(signif_12dvs1d_UP),
         nrow(signif_13dvs1d_UP))

vec_DOWN=c(nrow(signif_2dvs1d_DOWN),
           nrow(signif_3dvs1d_DOWN),
           nrow(signif_4dvs1d_DOWN),
           nrow(signif_5dvs1d_DOWN),
           nrow(signif_6dvs1d_DOWN),
           nrow(signif_7dvs1d_DOWN),
           nrow(signif_8dvs1d_DOWN),
           nrow(signif_9dvs1d_DOWN),
           nrow(signif_10dvs1d_DOWN),
           nrow(signif_11dvs1d_DOWN),
           nrow(signif_12dvs1d_DOWN),
           nrow(signif_13dvs1d_DOWN))
#####

# vs n-1
#####
signif_2dvs1d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_n-1/Patient2/TCseq_p2_2dvs1d_DOWN.txt",header = T)
signif_3dvs2d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_n-1/Patient2/TCseq_p2_3dvs2d_DOWN.txt",header = T)
signif_4dvs3d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_n-1/Patient2/TCseq_p2_4dvs3d_DOWN.txt",header = T)
signif_5dvs4d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_n-1/Patient2/TCseq_p2_5dvs4d_DOWN.txt",header = T)
signif_6dvs5d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_n-1/Patient2/TCseq_p2_6dvs5d_DOWN.txt",header = T)
signif_7dvs6d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_n-1/Patient2/TCseq_p2_7dvs6d_DOWN.txt",header = T)
signif_8dvs7d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_n-1/Patient2/TCseq_p2_8dvs7d_DOWN.txt",header = T)
signif_9dvs8d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_n-1/Patient2/TCseq_p2_9dvs8d_DOWN.txt",header = T)
signif_10dvs9d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_n-1/Patient2/TCseq_p2_10dvs9d_DOWN.txt",header = T)
signif_11dvs10d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_n-1/Patient2/TCseq_p2_11dvs10d_DOWN.txt",header = T)
signif_12dvs11d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_n-1/Patient2/TCseq_p2_12dvs11d_DOWN.txt",header = T)
signif_13dvs12d_DOWN=read.table("/Users/fla./Desktop/CLL_RNA-Seq/TCseq/Diff_Analysis_vs_n-1/Patient2/TCseq_p2_13dvs12d_DOWN.txt",header = T)

vec_UP=c(nrow(signif_2dvs1d_UP),
         nrow(signif_3dvs2d_UP),
         nrow(signif_4dvs3d_UP),
         nrow(signif_5dvs4d_UP),
         nrow(signif_6dvs5d_UP),
         nrow(signif_7dvs6d_UP),
         nrow(signif_8dvs7d_UP),
         nrow(signif_9dvs8d_UP),
         nrow(signif_10dvs9d_UP),
         nrow(signif_11dvs10d_UP),
         nrow(signif_12dvs11d_UP),
         nrow(signif_13dvs12d_UP))

vec_DOWN=c(nrow(signif_2dvs1d_DOWN),
           nrow(signif_3dvs2d_DOWN),
           nrow(signif_4dvs3d_DOWN),
           nrow(signif_5dvs4d_DOWN),
           nrow(signif_6dvs5d_DOWN),
           nrow(signif_7dvs6d_DOWN),
           nrow(signif_8dvs7d_DOWN),
           nrow(signif_9dvs8d_DOWN),
           nrow(signif_10dvs9d_DOWN),
           nrow(signif_11dvs10d_DOWN),
           nrow(signif_12dvs11d_DOWN),
           nrow(signif_13dvs12d_DOWN))

#####


#UP
b=barplot(vec_UP,col="#77dd77",names=c("2d","3d","4d","5d","6d","7d","8d","9d","10d","11d","12d","13d"),ylim=c(0,70),main="UP deregulated genes - Patient2 - ref 1d")
text(b, vec_UP, labels=as.character(vec_UP), cex= .9,pos=3)

#DOWN
b=barplot(-vec_DOWN,col="#ff6961",names=c("2d","3d","4d","5d","6d","7d","8d","9d","10d","11d","12d","13d"),ylim=c(-70,0),main="DOWN deregulated genes - Patient2 - ref 1d")
text(b, -vec_DOWN, labels=as.character(vec_DOWN), cex= .9,pos=1)


DBres.sig <- DBresult(tca_p1, group1 = "1d", group2 = c("2d"), top.sig = TRUE) #only logFC > 2 and pval < 0.05



### Biomart ###
library("biomaRt")
ensembl=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

UP=signif_13dvs1d_UP$id
DOWN=signif_13dvs1d_DOWN$id

#UP
genes.with.id=getBM(attributes=c("ensembl_gene_id","hgnc_symbol","ensembl_gene_id_version", "entrezgene_id"),values=UP, mart= ensembl)
id=genes.with.id[which(genes.with.id$hgnc_symbol %in% UP),]
myentrez_13dvs1d_UP=unique(id$entrezgene_id)

#DOWN
genes.with.id=getBM(attributes=c("ensembl_gene_id","hgnc_symbol","ensembl_gene_id_version", "entrezgene_id"),values=DOWN, mart= ensembl)
id=genes.with.id[which(genes.with.id$hgnc_symbol %in% DOWN),]
myentrez_13dvs1d_DOWN=unique(id$entrezgene_id)


# vs 1st day
#####
myentrez_list_UP=list(UP_2dvs1d=myentrez_2dvs1d_UP,UP_3dvs1d=myentrez_3dvs1d_UP,UP_4dvs1d=myentrez_4dvs1d_UP,
                      UP_5dvs1d=myentrez_5dvs1d_UP,UP_6dvs1d=myentrez_6dvs1d_UP,UP_7dvs1d=myentrez_7dvs1d_UP,
                      UP_8dvs1d=myentrez_8dvs1d_UP,UP_9dvs1d=myentrez_9dvs1d_UP,UP_10dvs1d=myentrez_10dvs1d_UP,
                      UP_11dvs1d=myentrez_11dvs1d_UP,UP_12dvs1d=myentrez_12dvs1d_UP,UP_13dvs1d=myentrez_13dvs1d_UP)


myentrez_list_DOWN=list(DOWN_2dvs1d=myentrez_2dvs1d_DOWN,DOWN_3dvs1d=myentrez_3dvs1d_DOWN,DOWN_4dvs1d=myentrez_4dvs1d_DOWN,
                        DOWN_5dvs1d=myentrez_5dvs1d_DOWN,DOWN_6dvs1d=myentrez_6dvs1d_DOWN,DOWN_7dvs1d=myentrez_7dvs1d_DOWN,
                        DOWN_8dvs1d=myentrez_8dvs1d_DOWN,DOWN_9dvs1d=myentrez_9dvs1d_DOWN,DOWN_10dvs1d=myentrez_10dvs1d_DOWN,
                        DOWN_11dvs1d=myentrez_11dvs1d_DOWN,DOWN_12dvs1d=myentrez_12dvs1d_DOWN,DOWN_13dvs1d=myentrez_13dvs1d_DOWN)




myentrez_list_ALL=list(UP_2dvs1d=c(myentrez_2dvs1d_UP,myentrez_2dvs1d_DOWN),UP_3dvs1d=c(myentrez_3dvs1d_UP,myentrez_3dvs1d_DOWN),UP_4dvs1d=c(myentrez_4dvs1d_UP,myentrez_4dvs1d_DOWN),
                       UP_5dvs1d=c(myentrez_5dvs1d_UP,myentrez_5dvs1d_DOWN),UP_6dvs1d=c(myentrez_6dvs1d_UP,myentrez_6dvs1d_DOWN),UP_7dvs1d=c(myentrez_7dvs1d_UP,myentrez_7dvs1d_DOWN),
                       UP_8dvs1d=c(myentrez_8dvs1d_UP,myentrez_8dvs1d_DOWN),UP_9dvs1d=c(myentrez_9dvs1d_UP,myentrez_9dvs1d_DOWN),UP_10dvs1d=c(myentrez_10dvs1d_UP,myentrez_10dvs1d_DOWN),
                       UP_11dvs1d=c(myentrez_11dvs1d_UP,myentrez_11dvs1d_DOWN),UP_12dvs1d=c(myentrez_12dvs1d_UP,myentrez_12dvs1d_DOWN),UP_13dvs1d=c(myentrez_13dvs1d_UP,myentrez_13dvs1d_DOWN))

#####

# vs n-1
#####
myentrez_list_UP=list(UP_2dvs1d=myentrez_2dvs1d_UP,UP_3dvs2d=myentrez_3dvs2d_UP,UP_4dvs3d=myentrez_4dvs3d_UP,
                      UP_5dvs4d=myentrez_5dvs4d_UP,UP_6dvs5d=myentrez_6dvs5d_UP,UP_7dvs6d=myentrez_7dvs6d_UP,
                      UP_8dvs7d=myentrez_8dvs7d_UP,UP_9dvs8d=myentrez_9dvs8d_UP,UP_10dvs9d=myentrez_10dvs9d_UP,
                      UP_11dvs10d=myentrez_11dvs10d_UP,UP_12dvs11d=myentrez_12dvs11d_UP,UP_13dvs12d=myentrez_13dvs12d_UP)

myentrez_list_DOWN=list(DOWN_2dvs1d=myentrez_2dvs1d_DOWN,DOWN_3dvs2d=myentrez_3dvs2d_DOWN,DOWN_4dvs3d=myentrez_4dvs3d_DOWN,
                      DOWN_5dvs4d=myentrez_5dvs4d_DOWN,DOWN_6dvs5d=myentrez_6dvs5d_DOWN,DOWN_7dvs6d=myentrez_7dvs6d_DOWN,
                      DOWN_8dvs7d=myentrez_8dvs7d_DOWN,DOWN_9dvs8d=myentrez_9dvs8d_DOWN,DOWN_10dvs9d=myentrez_10dvs9d_DOWN,
                      DOWN_11dvs10d=myentrez_11dvs10d_DOWN,DOWN_12dvs11d=myentrez_12dvs11d_DOWN,DOWN_13dvs12d=myentrez_13dvs12d_DOWN)


myentrez_list_ALL=list(UP_2dvs1d=c(myentrez_2dvs1d_UP,myentrez_2dvs1d_DOWN),UP_3dvs2d=c(myentrez_3dvs2d_UP,myentrez_3dvs2d_DOWN),UP_4dvs3d=c(myentrez_4dvs3d_UP,myentrez_4dvs3d_DOWN),
                       UP_5dvs4d=c(myentrez_5dvs4d_UP,myentrez_5dvs4d_DOWN),UP_6dvs5d=c(myentrez_6dvs5d_UP,myentrez_6dvs5d_DOWN),UP_7dvs6d=c(myentrez_7dvs6d_UP,myentrez_7dvs6d_DOWN),
                       UP_8dvs7d=c(myentrez_8dvs7d_UP,myentrez_8dvs7d_DOWN),UP_9dvs8d=c(myentrez_9dvs8d_UP,myentrez_9dvs8d_DOWN),UP_10dvs9d=c(myentrez_10dvs9d_UP,myentrez_10dvs9d_DOWN),
                       UP_11dvs10d=c(myentrez_11dvs10d_UP,myentrez_11dvs10d_DOWN),UP_12dvs11d=c(myentrez_12dvs11d_UP,myentrez_12dvs11d_DOWN),UP_13dvs12d=c(myentrez_13dvs12d_UP,myentrez_13dvs12d_DOWN))
#####

#GO
library(limma)
library(clusterProfiler)

#goana=goana(de=myentrez,geneid=myentrez,FDR=0.95,species="Hs")
#View(goana)

#GO Plot
compGO <- compareCluster(geneCluster     = myentrez_list_ALL,
                         fun           = "enrichGO",
                         OrgDb         = 'org.Hs.eg.db',
                         ont           = 'BP',
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")
dotplot(compGO, showCategory = 15, title = "Gene Ontology Enrichment Analysis")

#Temporal pattern analysis
#tca_p2 <- timecourseTable(tca_p2, value = "FC", norm.method = "tpm", filter = TRUE)
tca_p2 <- timecourseTable(tca_p2, value = "expression", norm.method = "cpm", filter = TRUE)
t <- tcTable(tca_p2)

#clustering
tca_p2 <- timeclust(tca_p2, algo = "cm", k = 6, standardize = TRUE)

p <- timeclustplot(tca_p2, value = "z-score(PRKM)", cols = 3)

#plot cluster 1:
print(p[[1]])

#see which cluster each gene belongs to
tca_p2@clusterRes@cluster

table(tca_p2@clusterRes@cluster)


saveRDS(tca_p2, file = "/Users/fla./Desktop/CLL_RNA-Seq/TCseq/TCA_patient2.rds")
write.table(tca_p2@clusterRes@cluster,"/Users/fla./Desktop/CLL_RNA-Seq/TCseq/TCA_cmeans_clusters_patient2.txt",row.names=T,col.names = F,quote=F,sep="\t")








######
#DESeq2
########
library(DESeq2)
######
#patient 1

#design
Time = rep(c("0h", "2h", "4h", "6h", "8h"), each = 2)
Treat = rep(c("Control", "Treat"), each = 10)
nameD <- paste(Treat, Time, c(rep(LETTERS[1:2], 5), rep(LETTERS[3:4], 5)), sep = "_")
coldata <- data.frame(row.names = nameD, Time = Time, Treat = Treat)


design=data.frame(row.names=c(colnames(p1_rep1)[2:12],colnames(p1_rep2)[2:12]),
                  timepoint=rep(c("1d","2d","3d","4d","5d","6d","7d","8d","9d","11d","13d"),2),
                  replicate=c(rep("rep1",11),rep("rep2",11)))
                  #treat=c("Control",rep("Treat",10)))




dds <- DESeqDataSetFromMatrix(countData=cm_p1, colData=design, design= ~ timepoint+replicate+timepoint:replicate)
dds <- DESeq(dds, test="LRT", reduced = ~ timepoint+replicate)
dds <- DESeq(dds, test="LRT")

res <- results(dds)












#vera request
library(TCseq)
#data=read.table("/Users/fla./Desktop/TCSeq/20180111 325 transcripto.xlsx - 20180111 325 transcripto.tsv",header = T)
subset=read.table("/Users/fla./Desktop/TCSeq/subset.txt")$V1
data=read.table("/Users/fla./Downloads/AraC_3PDXs_CinétiqueW0-1-2-4.xlsx - Feuil1.tsv",header=T,fill=T)
data=data[which(!(is.na(data$X53tx4wk))),]
View(data)

#replace "," by "."
library(stringr)
for (i in 3:7){
  data[,i]=as.double(str_replace(data[,i],",","."))
}

#average the genes duplicates
out <- aggregate(. ~ Gene_Symbol, data = data[2:14], mean)


##### create TCA object
### design table
design=data.frame(sampleid=colnames(data)[3:ncol(data)],
                  timepoint=rep(c("untx","tx1w","tx2w","tx3w"),each=3),
                  group=rep(1:4,each=3))


### Count matrix
cm=out[,2:13];rownames(cm)=out$Gene_Symbol
cm=as.matrix(cm)
cm=round(cm)


### genomic intervals
# Parse Ensembl genes
library(EnsDb.Hsapiens.v75)
Ens=as.data.frame(genes(EnsDb.Hsapiens.v75)) #hg19
#Ens=Ens[which(Ens$gene_biotype=="protein_coding"),]
Ens$seqnames=paste0("chr",Ens$seqnames)
rownames(Ens)=NULL
Ens=Ens[,c(1,2,3,7)]
colnames(Ens)=c("chr","start","end","id")
chr=c(paste0("chr",1:22),"chrX")
Ens=Ens[which(Ens$chr %in% chr),]

# Select genes
genes=data.frame()
for (i in rownames(cm)){
  genes=rbind(genes,Ens[which(Ens$id==i),])
}

# keep first of duplicated genes
remove=vector()
for (i in unique(genes$id[which(isUnique(genes$id)==FALSE)])){
  remove=c(remove,rownames(genes[which(genes$id==i),])[2:nrow(genes[which(genes$id==i),])])
}
remove=as.numeric(remove)
genes=genes[which(!(rownames(genes) %in% remove)),]

cm=cm[which(rownames(cm) %in% genes$id),]

# Verif
table(genes$id[which(isUnique(genes$id)==FALSE)]) #should be empty
dim(genes)
dim(cm)
table(rownames(cm)==genes$id) #should be only TRUE

### TCA object
tca <- TCA(design = design, genomicFeature = genes, counts = cm)

##### Diff Analysis
#tca_p1 <- DBanalysis(tca_p1)
tca <- DBanalysis(tca, filter.type = "cpm", filter.value = 5, samplePassfilter = 2) #filter low read counts (at least 5 in 2 samples)
tca <- DBanalysis(tca, filter.type = "raw", filter.value = 1, samplePassfilter = 1) #filter low read counts (at least 5 in 2 samples)
tca <- DBanalysis(tca, filter.type = "NULL", filter.value = 5, samplePassfilter = 2) #filter low read counts (at least 5 in 2 samples)

#Extract results
DBres <- DBresult(tca, group1 = "1d", group2 = c("13d"))
res=as.data.frame(DBres$`13dvs1d`)



#Temporal pattern analysis
tca <- timecourseTable(tca, value = "FC", norm.method = "tpm", filter = TRUE)
tca <- timecourseTable(tca, value = "expression", norm.method = "cpm", filter = F)
t <- tcTable(tca)

#clustering
tca_3c <- timeclust(tca, algo = "cm", k = 3, standardize = T)
p <- timeclustplot(tca_3c, value = "Z-Score", cols = 3)
table(tca_3c@clusterRes@cluster)
write.table(tca_3c@clusterRes@cluster,"/Users/fla./Desktop/TCSeq/TCA_cmeans_3custer.txt",row.names=T,col.names = F,quote=F,sep="\t")
tab3=tca_3c@clusterRes@cluster
tab3=tab3[which(names(tab3) %in% subset)]
write.table(tab3,"/Users/fla./Desktop/TCSeq/TCA_cmeans_3clusters_subset.txt",row.names=T,col.names = F,quote=F,sep="\t")

tca_6c <- timeclust(tca, algo = "cm", k = 6, standardize = T)
p <- timeclustplot(tca_6c, value = "Z-Score", cols = 3)
table(tca_6c@clusterRes@cluster)
write.table(tab6,"/Users/fla./Desktop/TCSeq/TCA_cmeans_6clusters_subset.txt",row.names=T,col.names = F,quote=F,sep="\t")
tab6=tca_6c@clusterRes@cluster
tab6=tab6[which(names(tab6) %in% subset)]

#plot cluster 1:
print(p[[1]])

#see which cluster each gene belongs to
tca@clusterRes@cluster

table(tca@clusterRes@cluster)

write.table(tca_p1@clusterRes@cluster,"/Users/fla./Desktop/CLL_RNA-Seq/TCseq/TCA_cmeans_cluster_patient1.txt",row.names=T,col.names = F,quote=F,sep="\t")




#PCA
#PCA
library(ggplot2)
library(RColorBrewer)
display.brewer.all()

J_rep1=c(paste0("W",0:4))
J_rep1=factor(J_rep1,levels=c("W0","W1","W2","W3","W4"))

#J_rep2=c("J1","J2","J3","J4","J5","J6","J7","J8","J9","J10","J11","J12","J13")
#J_rep2=factor(J_rep2,levels=c("J1","J2","J3","J4","J5","J6","J7","J8","J9","J10","J11","J12","J13"))

PCA <- prcomp(t(out[2:ncol(out)]), scale = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = J_rep1)) +
  geom_text(aes(label=J_rep1),hjust=-.2, vjust=.4) +
  ggtitle("PCA - Patient 2 rep 1") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = rainbow(length(percentVar)))

library(gplots)
heatmap.2(cm)

subset=read.table("/Users/fla./Desktop/TCSeq/subset.txt")$V1
cm_s=cm[which(rownames(cm) %in% subset),]
out_s=out[which(out$GENE_ID %in% subset),]

heatmap.2(cm_s)

PCA <- prcomp(t(out_s[2:ncol(out_s)]), scale = FALSE)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(colour = J_rep1)) +
  geom_text(aes(label=J_rep1),hjust=-.2, vjust=.4) +
  ggtitle("PCA - Patient 2 rep 1") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = rainbow(length(percentVar)))






### maSigPro
library(maSigPro)

Time=0:4
Replicates=rep(1,5)
Group=rep(1,5)
ss.edesign <- cbind(Time,Replicates,Group)
rownames(ss.edesign) <- paste("Array", c(1:5), sep = "")
ss.GENE <- function(n, r, var11 = 0.01, var12 = 0.02, var13 = 0.02,
                    var14 = 0.02, a1 = 0, a2 = 0, a3 = 0, a4 = 0) {
  tc.dat <- NULL
  for (i in 1:n) {
    gene <- c(rnorm(r, a1, var11), rnorm(r, a1, var12),
              rnorm(r, a3, var13), rnorm(r, a4, var14))
    tc.dat <- rbind(tc.dat, gene)
  }
  tc.dat }
flat <-ss.GENE(n = 85, r = 3) # flat
induc <- ss.GENE(n = 5, r = 3, a1 = 0, a2 = 0.2, a3 = 0.6, a4 = 1) # induction
sat <- ss.GENE(n = 5, r = 3, a1 = 0, a2 = 1, a3 = 1.1, a4 = 1.2) # saturation
ord <- ss.GENE(n = 5, r = 3, a1 = -0.8, a2 = -1, a3 = -1.3, a4 =-0.9) # intercept
ss.DATA <- rbind(flat, induc,sat,ord)
rownames(ss.DATA) <- paste("feature", c(1:100), sep = "")
colnames(ss.DATA) <- paste("Array", c(1:12), sep = "")
# run maSigPro
ss.example <- maSigPro(ss.DATA, ss.edesign, vars="each")


## make a single series edesign
Time <- rep(c(1,5,10,24), each = 3)
Replicates <- rep(c(1:4), each = 3)
Group <- rep(1,12)

ss.edesign <- cbind(Time,Replicates,Group)
rownames(ss.edesign) <- paste("Array", c(1:12), sep = "")
## Create data set
ss.GENE <- function(n, r, var11 = 0.01, var12 = 0.02, var13 = 0.02,
                      var14 = 0.02, a1 = 0, a2 = 0, a3 = 0, a4 = 0) {
  tc.dat <- NULL
  for (i in 1:n) {
    gene <- c(rnorm(r, a1, var11), rnorm(r, a1, var12),
              rnorm(r, a3, var13), rnorm(r, a4, var14))
    tc.dat <- rbind(tc.dat, gene)
  }
  tc.dat }
flat <-ss.GENE(n = 85, r = 3) # flat
induc <- ss.GENE(n = 5, r = 3, a1 = 0, a2 = 0.2, a3 = 0.6, a4 = 1) # induction
sat <- ss.GENE(n = 5, r = 3, a1 = 0, a2 = 1, a3 = 1.1, a4 = 1.2) # saturation
ord <- ss.GENE(n = 5, r = 3, a1 = -0.8, a2 = -1, a3 = -1.3, a4 =-0.9) # intercept
ss.DATA <- rbind(flat, induc,sat,ord)
rownames(ss.DATA) <- paste("feature", c(1:100), sep = "")
colnames(ss.DATA) <- paste("Array", c(1:12), sep = "")
# run maSigPro
ss.example <- maSigPro(ss.DATA, ss.edesign, vars="each")







#dysGENIE3


library(doRNG)
library(doParallel)
library(reshape2)


TS1 <- read.expr.matrix("example_data/time_series_1.txt",form="rows.are.samples")
TS2 <- read.expr.matrix("example_data/time_series_2.txt",form="rows.are.samples")
TS3 <- read.expr.matrix("example_data/time_series_3.txt",form="rows.are.samples")
time.points <- list(TS1[1,], TS2[1,], TS3[1,])
TS.data <- list(TS1[2:nrow(TS1),], TS2[2:nrow(TS2),], TS3[2:nrow(TS3),])



p1_rep1=read.table("/Volumes/CRCT21/Private/CLL_dec2020/RNAseqBulk_Eq21_Eq9/Gene_counts_tables/Gene_counts_Patient1_rep1_TPM.txt",header = T)

p1_rep1=read.expr.matrix("/Volumes/CRCT21/Private/CLL_dec2020/RNAseqBulk_Eq21_Eq9/Gene_counts_tables/Gene_counts_Patient1_rep1_TPM.txt",form="rows.are.genes")
p1_rep2=read.expr.matrix("/Volumes/CRCT21/Private/CLL_dec2020/RNAseqBulk_Eq21_Eq9/Gene_counts_tables/Gene_counts_Patient1_rep2_TPM.txt",form="rows.are.genes")

p1_rep1=p1_rep1[1:10,]
p1_rep2=p1_rep2[1:10,]

t1=c(1:11);names(t1)=colnames(p1_rep1)
t2=c(1:11);names(t2)=colnames(p1_rep1)
colnames(p1_rep2)=colnames(p1_rep1)

time.points <- list(t1,t2)
TS.data <- list(p1_rep1,p1_rep2)

res <- dynGENIE3(TS.data,time.points,regulators="")

write.table(p1_rep2,"/Volumes/CRCT21/Private/CLL_dec2020/RNAseqBulk_Eq21_Eq9/Gene_counts_tables/Gene_counts_Patient1_rep2_TPM.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(p2_rep1,"/Volumes/CRCT21/Private/CLL_dec2020/RNAseqBulk_Eq21_Eq9/Gene_counts_tables/Gene_counts_Patient2_rep1_TPM.txt",col.names=T,row.names=F,quote=F,sep="\t")
write.table(p2_rep2,"/Volumes/CRCT21/Private/CLL_dec2020/RNAseqBulk_Eq21_Eq9/Gene_counts_tables/Gene_counts_Patient2_rep2_TPM.txt",col.names=T,row.names=F,quote=F,sep="\t")



## TCSeq clusters / deregulated genes

c1=read.table("/Users/fla./Desktop/TCSeq/Macrophage_timecourse/TCA_cmeans_2cluster.txt")
c1=c1[which(c1$V2==1),1]
c2=c1[which(c1$V2==2),1]


# M1 vs Mono
m1_vs_mono=read.table("/Users/fla./Downloads/DE_All_gene_tab_M1_vs_Monocyte.txt",header = T)
m1_vs_mono_all=m1_vs_mono[which(m1_vs_mono$adj.P.Val<0.1),1]
m1_vs_mono_up=m1_vs_mono[which(m1_vs_mono$adj.P.Val<0.1 & m1_vs_mono$logFC>0),1]
m1_vs_mono_down=m1_vs_mono[which(m1_vs_mono$adj.P.Val<0.1 & m1_vs_mono$logFC<0),1]

# m2 vs Mono
m2_vs_mono=read.table("/Users/fla./Downloads/ALL_DE_tab_M2_vs_Monocyte.txt",header = T)
m2_vs_mono_all=m2_vs_mono[which(m2_vs_mono$adj.P.Val<0.1),1]
m2_vs_mono_up=m2_vs_mono[which(m2_vs_mono$adj.P.Val<0.1 & m2_vs_mono$logFC>0),1]
m2_vs_mono_down=m2_vs_mono[which(m2_vs_mono$adj.P.Val<0.1 & m2_vs_mono$logFC<0),1]

# Mono vs NLC
nlc_vs_mono=read.table("/Users/fla./Downloads/ALL_DE_tab_Monocyte_vs_NLC.txt",header = T)
nlc_vs_mono_all=nlc_vs_mono[which(nlc_vs_mono$adj.P.Val<0.1),1]
nlc_vs_mono_up=nlc_vs_mono[which(nlc_vs_mono$adj.P.Val<0.1 & nlc_vs_mono$logFC<0),1]
nlc_vs_mono_down=nlc_vs_mono[which(nlc_vs_mono$adj.P.Val<0.1 & nlc_vs_mono$logFC>0),1]

# Mono vs TAM
tam_vs_mono=read.table("/Users/fla./Downloads/ALL_DE_tab_Monocyte_vs_TAM.txt",header = T)
tam_vs_mono_all=tam_vs_mono[which(tam_vs_mono$adj.P.Val<0.1),1]
tam_vs_mono_up=tam_vs_mono[which(tam_vs_mono$adj.P.Val<0.1 & tam_vs_mono$logFC<0),1]
tam_vs_mono_down=tam_vs_mono[which(tam_vs_mono$adj.P.Val<0.1 & tam_vs_mono$logFC>0),1]

# NLC vs M1+M2
nlc_vs_mac=read.table("/Users/fla./Desktop/Vera_table/NLC_vs_M1+M2/Macrophage_vs_NLC_ALL.txt",header = T)
nlc_vs_mac_all=rownames(nlc_vs_mac[which(nlc_vs_mac$adj.P.Val<0.1),])
nlc_vs_mac_up=rownames(nlc_vs_mac[which(nlc_vs_mac$adj.P.Val<0.1 & nlc_vs_mac$logFC>0),])
nlc_vs_mac_down=rownames(nlc_vs_mac[which(nlc_vs_mac$adj.P.Val<0.1 & nlc_vs_mac$logFC<0),])

#piechart
ol_up=table(c1 %in% nlc_vs_mac_up)['TRUE']/length(c1)
ol_down=table(c1 %in% nlc_vs_mac_down)['TRUE']/length(c1)
vec=c(1-(ol_up+ol_down),ol_up,ol_down)
pie(vec,border="white",col=c("grey85","#598ec7","#b2332d"),labels=c("No Overlaps\n35.5%","\n\nUP\n7%","\n\nDOWN\n57.5%"),main="Overlap\nCluster1 / DE genes NLC vs M1+M2")


library(EnsDb.Hsapiens.v75)
genes=as.data.frame(genes(EnsDb.Hsapiens.v75))
nlc_vs_mac=read.table("/Users/fla./Desktop/Vera_table/NLC_vs_M1+M2/Macrophage_vs_NLC_ALL.txt",header = T)
nlc_vs_mac=merge(nlc_vs_mac,genes,by.x=0,by.y="gene_name")
