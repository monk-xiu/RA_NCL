options(stringsAsFactors = FALSE)

# 创建一个空的列表，用于存储所有的数据框
all_dfs <- list()

# 文件夹路径
folder_path <- "E:/QF/sh"

# 列出文件夹中的所有txt文件
txt_files <- list.files(folder_path, pattern = "\\.txt$", full.names = TRUE)

# 循环读取所有的txt文件
for (file_path in txt_files) {
  # 提取文件名作为数据框的名称
  df_name <- tools::file_path_sans_ext(basename(file_path))
  
  # 检查文件名是否以 .summary 结尾
  if (endsWith(df_name, ".summary")) {
    next  # 跳到下一个循环迭代
  }
  
  # 读取文件内容
  df <- read.table(file_path, header = TRUE)
  # 提取第一列和最后一列
  df_subset <- df[, c(1, ncol(df))]
  # 将数据框添加到列表中，使用文件名作为列表的名称
  all_dfs[[df_name]] <- df_subset
}


#在这个示例中，我们使用 Reduce() 函数和 merge() 函数按照每个子数据框的 "Geneid" 列进行合并，
#all = TRUE 用于保留所有匹配的行，
#即使在某些数据框中可能不存在匹配。最后，merged_df 包含了所有数据框按照 "Geneid" 列合并后的结果。
merged_df <- Reduce(function(x, y) merge(x, y, by = 1, all = TRUE), all_dfs)

write.csv(merged_df,file = 'allcount.csv')

merged_df <- read.csv('allcount.csv',row.names = 1)
# 注释
library(refGenome)
ens <- ensemblGenome()
read.gtf(ens, "../Homo_sapiens.GRCh38.105.chr.gtf")###导入gtf文件 比较耗时
my_gene <- getGenePositions(ens)
colnames(my_gene)
gene_id<- my_gene$gene_id
gene_name <- my_gene$gene_name
gene_biotype <- my_gene$gene_biotype

a <- cbind(gene_id,gene_name, gene_biotype)


write.table(a,file="gene_anno.txt",sep="\t",quote=F,row.names = T) ###保存一下

#a <- read.table("../../gene_anno.txt")
sample <- read.csv('allcount.csv',header = T, row.names = 1)           

alldata <- merge(sample,a, by.x = 1,by.y = 1)
write.table(alldata,file="all_exp_anno.txt",sep="\t",quote=F,row.names = T) ###保存一下

# library(tidyverse)
# alldata2 <- alldata %>% filter(OA_1.sort.bam>0)
# write.table(alldata2,file="OA2.txt",sep="\t",quote=F,row.names = F) ###保存一下



######################差异分析#################################
library(DESeq2)
dat <- read.table("all_exp_anno.txt",sep="\t")
#或者dat <- as.matrix(read.csv('allcount.csv',header = T,row.names = 'Geneid'))
dat <- dat[rowMeans(dat[,2:7])>1,]# 去除表达量过低的基因

library(tidyverse)
dat <- dat[!duplicated(dat$gene_name),] %>% na.omit()
rownames(dat) <- dat$gene_name
dat <- dat[,-c(1,8,9)]
colnames(dat) <- c('sh1','sh2','sh3','NC1','NC2','NC3')
write.table(alldata,file="all_exp_anno2.txt",sep="\t",quote=F,row.names = T) ###保存一下

dat <- as.matrix(dat)




colData <- data.frame(row.names=colnames(dat), condition = factor(rep(c("sh", "NC"), each = 3)))
#sampleTable <- data.frame(condition = factor(rep(c("OA", "RA"), each = 5)))
#sampleTable$condition <- factor(sampleTable$condition, levels = c("OA", "RA"))

dds <- DESeqDataSetFromMatrix(count=dat, colData=colData, design=~condition)
dds <- DESeq(dds)

#####################contrast后面三个参数，先是列名，再是分子和分母
res_shvsNC <- results(dds , contrast = c("condition","sh","NC"),tidy=TRUE)
#或者
deseq2.obj <- DESeq(dds)
res <- results(deseq2.obj)
res.df <- as.data.frame(res)

summary(res_shvsNC)
write.csv(res_shvsNC,file = 'alldf.csv')


### pca 作图
library(ggfortify)
## normalization
count.sf <-  estimateSizeFactors(dds)
normal.mtx <- counts(count.sf, normalized=T)
df_all <- cbind(t(normal.mtx), data.frame(group=c(rep('sh', 3),rep('NC', 3))))

pca_data <- prcomp(t(normal.mtx ), scale. = F)

library(ggsci)
pdf('pca.pdf',width = 8,height = 7)
autoplot(pca_data, data = df_all, colour = 'group', label = TRUE) + scale_color_aaas() +  theme_bw()
dev.off()


keep <- abs(res_shvsNC$log2FoldChange) > 1 & res_shvsNC$pvalue < 0.05                                                                                                                                                                        
deg <- res_shvsNC[keep, ]   %>% na.omit()                                                                                                                                                                                                        
dim(deg)                                                                                                                                                                                                                        
write.csv(deg,file = 'deg.csv')

######火山图######
volcano_plot=function(logfc,pvalue,symbol=NULL,cutFC=1,cutPvalue=0.05
                      ,showText=NULL
                      ,colors=c('red','grey','blue')
                      ,ylab='-log10(FDR)',
                      leg='Group',
                      xlab='log2(FoldChange)'){
  library(ggplot2)
  cange=rep('None',length(logfc))
  cange[which(logfc>cutFC&pvalue<cutPvalue)]='Up'
  cange[which(logfc< -cutFC&pvalue<cutPvalue)]='Down'
  if(is.null(symbol)){
    symbol=rep('',length(logfc))
    showText=NULL
  }
  dat.input=data.frame(logFC=logfc,FDR=pvalue,change=cange,SYMBOL=symbol)
  #print(head(dat.input))
  p1 <- ggplot(data = dat.input, 
               aes(x = logFC, 
                   y = -log10(FDR)))
  p1=p1+geom_point(alpha=0.4, size=3.5, aes(color=change))  
  p1=p1+scale_color_manual(values=colors,limits = c("Down",'None', "Up"),name=leg) 
  p1=p1+geom_vline(xintercept=c(-cutFC,cutFC),lty=2,col="black",lwd=0.8)  
  p1=p1+geom_hline(yintercept = -log10(cutPvalue),lty=2,col="black",lwd=0.8)  
  p1=p1+ylab(ylab)+xlab(xlab)
  p1=p1+theme_classic()
  p1=p1+theme(legend.background = element_rect(fill = NA, colour = NA))
  if(is.null(showText)|is.null(symbol)){
    showText=c()
  }
  if(length(showText)>0){
    for_label <-dat.input[match(intersect(showText,dat.input$SYMBOL),dat.input$SYMBOL),]
    p1=p1+geom_point(size = 3, shape = 1, data = for_label)+
      ggrepel::geom_label_repel(aes(label = SYMBOL),data = for_label,color="black" )
  }
  
  return(p1)
}

#先假设 dat 数据框包含了基因名字的一列，名为 "gene_name"
genes_to_label <- c("S100A9")  # 这里填你要标注的基因名

# 确保你要标注的基因名存在于数据中
if (!all(genes_to_label %in% dat$gene_name)) {
  stop("Not all genes to label are in the data!")
}

volcano <- volcano_plot(logfc = res_shvsNC$log2FoldChange,
                        pvalue = res_shvsNC$pvalue,
                        symbol = res_shvsNC$gene_name,  # 确保 "symbol" 参数指向包含基因名字的列
                        cutFC = log2(2), cutPvalue = 0.05,
                        colors = c('blue','grey','red'), 
                        leg = 'sh vs NC',
                        showText = genes_to_label)  # 把要标注的基因名传递给 "showText" 参数

ggsave('volcano.pdf',volcano,height = 4,width = 5)

#####heatmap#####
#devtools::install_github("stemangiola/tidyHeatmap")
library(tidyHeatmap)
library(tidyverse)

count <- read.table('all_exp_anno.txt')
count <- count[rowMeans(count[,2:15])>1,]
count <- count[!duplicated(count$gene_name),] %>% na.omit()
rownames(count) <- count$gene_name
count <- count[,-c(1,16,17)]
colnames(count) <- c('R1','R2','R3','R4','R5','R6','R7','NR1','NR2','NR3','NR4','NR5','NR6','NR7')
write.csv(count,'./count/count_anno.csv')

#leukocyte
leukocyte <- read.table('./heatmap/leukocyte chemotaxis.txt')
dat2 <- count[leukocyte$V1,]

#dat4['SNORD3',] <- NA
#dat4['SNORD3',1:6] <- colMeans(dat4[c('SNORD3B-2','SNORD3A','SNORD3C','SNORD3B-1','U3'),1:6]) #这里整合了一些基因

# 两种去除特定行的方法
# dat4 <- dat4[-c('SNORD3B-2','SNORD3A','SNORD3C','U3')] 这样删除行不可取
rownames_to_remove <- c('SNORD3B-2','SNORD3A','SNORD3C','SNORD3B-1','U3')
dat4 <- dat4[!rownames(dat4) %in% rownames_to_remove,]
# 或
# rows <- which(rownames(df) %in% c("b","d","f") )
# df <- df[-rows,]
dat4 <- dat4 %>% select(c('RA-1','RA-2','RA-3'),everything())

leukocyte <- leukocyte$V1
dat2 <- cbind(leukocyte,dat2)
colnames(dat2) <- c('leukocyte chemotaxis','R1','R2','R3','R4','R5','R6','R7','NR1','NR2','NR3','NR4','NR5','NR6','NR7')

dat_tidy <- 
  dat2 |> 
  # tidyfy
  pivot_longer(cols = c('R1','R2','R3','R4','R5','R6','R7','NR1','NR2','NR3','NR4','NR5','NR6','NR7'), names_to = "sample", values_to = "Value")
dat_tidy$group <- ifelse(grepl('^R.*',dat_tidy$sample),'R','NR')

dat_heatmap <- 
  dat_tidy |> 
  group_by(group) |> 
  heatmap(leukocyte, sample, Value,    scale = "row" ,palette_value = c( "blue", "white","red")  ) |>
  annotation_tile(group)


dat_heatmap <- 
  dat_tidy |> 
  heatmap(leukocyte, sample, Value,    scale = "row" )

dat_heatmap <- 
  dat_tidy |> 
  group_by(group) |> 
  heatmap(leukocyte, sample, Value,    scale = "row" ,palette_value = c( "blue", "white","red"),
          palette_grouping = list( # 可自定义注释条颜色
            # For second grouping (property_group)
            c("red", "blue")
          )  )

#添加分组框
pdf('./heatmap/in_heatmap.pdf',width = 6)
dat_heatmap <- 
  dat_tidy |> 
  group_by(group) |> 
  heatmap(`leukocyte chemotaxis`, sample, Value,    scale = "row" ,palette_value = c( "#84B1FB", "#fefada",'#E16867'),
          palette_grouping = list( # 可自定义注释条颜色
            # For second grouping (property_group)
            c('#E16867', "#84B1FB")
          )  )
dat_heatmap
dev.off()

#不添加分组框
pdf('in_heatmap.pdf',width = 6)
dat_heatmap <- 
  dat_tidy |> 
  heatmap(leukocyte, sample, Value,    scale = "row" ,palette_value = c("#fefada", "#21908CFF",  "#440154FF") )
dat_heatmap
dev.off()


###############富集分析#########
library("clusterProfiler")
library("org.Mm.eg.db")#小鼠
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(DOSE)
library(topGO)
#BiocManager::install("org.Mm.eg.db")

#com_gene=unique(c(rownames(tcga.diff.fit1),rownames(tcga.diff.fit2),rownames(tcga.diff.fit3)))

gene_ues <- dat %>% filter(pvalue<0.05)%>% filter(abs(log2FoldChange)>1)

gene_symbol <- gene_ues$gene_name
head(gene_symbol,2)
gene = bitr(gene_symbol, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库


#GO富集分析
ego <- enrichGO(gene = gene$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                pvalueCutoff =0.1, 
                #qvalueCutoff = 0.05,
                ont="all",
                readable =T)

write.table(ego@result,file="GO.txt",sep="\t",quote=F,row.names = F)
#KEGG富集
R.utils::setOption("clusterProfiler.download.method",'auto') 
kk <- enrichKEGG(gene = gene$ENTREZID,
                 organism = 'hsa', # 小鼠
                 pvalueCutoff = 0.05)
write.table(kk@result,file="kegg.txt",sep="\t",quote=F,row.names = F)

res<-rbind.data.frame(ego@result,
                      data.frame(ONTOLOGY='kegg',kk@result))
head(res)
table(res$ONTOLOGY)
res.fit=res[which(res$pvalue<0.05),]
table(res.fit$ONTOLOGY)
# BP kegg   MF 
# 2    2    6
write.table(res.fit,'enrich.txt',quote = F,row.names = F,sep='\t')

ggsci::pal_jama()(9)
ontology.col=ggsci::pal_lancet()(7)[1:4]

data=res.fit[order(res.fit$p.adjust),]
datasig=data[data$p.adjust<0.05,,drop=F]
BP = datasig[datasig$ONTOLOGY=="BP",,drop=F]
CC = datasig[datasig$ONTOLOGY=="CC",,drop=F]
MF = datasig[datasig$ONTOLOGY=="MF",,drop=F]
kegg = datasig[datasig$ONTOLOGY=="kegg",,drop=F]

BP = head(BP,5)
CC = head(CC,5)
MF = head(MF,5)
kegg=head(kegg,5)
data = rbind(BP,CC,MF,kegg)
main.col = ontology.col[c(rep(1,5),rep(2,5),rep(3,5),rep(4,5))]
#
BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
library(RColorBrewer)
library(circlize)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(GO=data$ID,start=1,end=max(BgGene))
rownames(df) = df$GO
bed2 = data.frame(GO=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(GO=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(GO=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5
library(ComplexHeatmap)
pdf('enrich.pdf',height = 9,width = 13)
{
  par(omi=c(0.1,0.1,0.1,1.5))
  circos.par(track.margin=c(0.01,0.01))
  circos.genomicInitialize(df,plotType="none")
  circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
  }, track.height = 0.08, bg.border = NA,bg.col = main.col)
  
  for(si in get.all.sector.index()) {
    circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
                major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
  }
  f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
  circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                      panel.fun = function(region, value, ...) {
                        i = getI(...)
                        circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                           border = NA, ...)
                        circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                      })
  circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                      panel.fun = function(region, value, ...) {
                        i = getI(...)
                        circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                           border = NA, ...)
                        circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                      })
  circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                      panel.fun = function(region, value, ...) {
                        cell.xlim = get.cell.meta.data("cell.xlim")
                        cell.ylim = get.cell.meta.data("cell.ylim")
                        for(j in 1:9) {
                          y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                          circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                        }
                        circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                           border = NA, ...)
                        #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                      })
  circos.clear()
  #绘制圈图中间的图例
  middle.legend = Legend(
    labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
    type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',ontology.col[1])),
    title="",nrow=3,size= unit(3, "mm")
  )
  circle_size = unit(1, "snpc")
  draw(middle.legend,x=circle_size*1.2,y=circle_size*0.3)
  #绘制GO分类的图例
  main.legend = Legend(
    labels = c("Biological Process","Cellular Component", "Molecular Function",'KEGG Pathway'),  type="points",pch=15,
    legend_gp = gpar(col=ontology.col), title_position = "topcenter",
    title = "ONTOLOGY", nrow = 4,size = unit(3, "mm"),grid_height = unit(5, "mm"),
    grid_width = unit(5, "mm")
  )
  #绘制pvalue的图例
  logp.legend = Legend(
    labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
    type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
    title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
    size = unit(3, "mm")
  )
  lgd = packLegend(main.legend,logp.legend)
  circle_size = unit(1, "snpc")
  #print(circle_size)
  draw(lgd, x = circle_size*1.2, y=circle_size*0.5,just = "left")
}
dev.off()


######自定义富集分析#########
library(clusterProfiler)
gmt.df=read.gmt("../c2.all.v2023.1.Hs.symbols.gmt") #gmt是一种存储基因集的格式，可以用文本编辑器直接打开
# UP
down_gene <- read.table('DOWN.txt')
SH_down_gmt=clusterProfiler::enricher(UP_gene$V1,TERM2GENE = gmt.df,minGSSize = 2,maxGSSize = 2000)
SH_down_res=SH_down_gmt@result
write.csv(OE_up_res,'SH_down_res.csv')



########GSEA#######
deg <- read.csv('alldf.csv',row.names = 1)
gsea_input <- arrange(deg,desc(log2FoldChange))
#准备目的gene列表；
gene_list <- gsea_input$log2FoldChange
names(gene_list) <- gsea_input$row
gene_list[1:6]

#读取从GSEA官网下载好的基因集文件；
h.all_gmt <- read.gmt("c2.all.v2023.1.Hs.symbols.gmt")

#预览转换格式后的gmt文件；
head(h.all_gmt)

#执行GSEA;
h.all_res <- GSEA(gene_list,TERM2GENE = h.all_gmt,pvalueCutoff = 1,verbose = F)
#预览分析结果；
head(h.all_res@result[1:10])
save(h.all_res,file = './gsea/all_gmt.Rdata')
load('./gsea/all_gmt.Rdata')


write.csv(h.all_res@result,file = 'all_gesa.csv')

#对感兴趣基因集的富集结果进行可视化；
gset="HALLMARK_OXIDATIVE_PHOSPHORYLATION"
gsetid<- h.all_res@result$ID
#获取位置编号；
index <- which(gsetid==gset)
#绘图；
gseaplot2(h.all_res, title = h.all_res$Description[index],
          geneSetID = index)


############go_gsea#############
#读取从GSEA官网下载好的基因集文件；
go_gsea <- read.gmt("c5.go.v2023.1.Hs.symbols.gmt")

#执行GSEA;
go_gsea_res <- GSEA(gene_list,TERM2GENE = go_gsea,pvalueCutoff = 1,verbose = F)
#预览分析结果；
head(go_gsea_res@result[1:10])
save(go_gsea_res,file = './gsea/go_gsea_res.Rdata')
load('./gsea/go_gsea_res.Rdata')
write.csv(go_gsea_res@result,file = 'go_gesa.csv')

#devtools::install_github("junjunlab/GseaVis")
library(GseaVis)
#自定义配色2：
#对感兴趣基因集的富集结果进行可视化；
kgset="GOBP_NEUTROPHIL_MEDIATED_KILLING_OF_SYMBIONT_CELL"
kgsetid<- go_gsea_res@result$ID
#获取位置编号；
index2 <- which(kgsetid==kgset)


gseaNb(object = go_gsea_res,
       geneSetID = go_gsea_res@result$ID[index2],
       subPlot = 3,
       addPval = T,
       pvalX = 0.95,
       pvalY = 0.8,
       curveCol = c("#99CC00","#c77cff"),
       htCol = c("#99CC00","#c77cff"),
       rankCol = c("#99CC00", "white", "#c77cff"))
ggsave('GOBP_NEUTROPHIL_MEDIATED_KILLING_OF_SYMBIONT_CELL.pdf',width = 6, height = 5)

gene_sets <- c('GOBP_NEUTROPHIL_MEDIATED_KILLING_OF_SYMBIONT_CELL',
               'GOBP_NEUTROPHIL_MEDIATED_CYTOTOXICITY',
               'GOBP_NEUTROPHIL_MEDIATED_IMMUNITY',
               'GOBP_NEUTROPHIL_CHEMOTAXIS',
               'GOBP_NEUTROPHIL_MIGRATION',
               'GOBP_NEUTROPHIL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
               'GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION',
               'GOBP_NEUTROPHIL_DEGRANULATION',
               'GOBP_REGULATION_OF_NEUTROPHIL_MIGRATION',
               'GOBP_NEUTROPHIL_EXTRAVASATION',
               'GOBP_REGULATION_OF_NEUTROPHIL_CHEMOTAXIS',
               'GOBP_NEUTROPHIL_HOMEOSTASIS',
               'GOBP_POSITIVE_REGULATION_OF_NEUTROPHIL_MIGRATION',
               'GOBP_NEUTROPHIL_DIFFERENTIATION')
for (i in gene_sets) {
  kgset= i
  kgsetid<- go_gsea_res@result$ID
  #获取位置编号；
  index2 <- which(kgsetid==kgset)
  gseaNb(object = go_gsea_res,
         geneSetID = go_gsea_res@result$ID[index2],
         subPlot = 3,
         addPval = T,
         pvalX = 0.95,
         pvalY = 0.8,
         curveCol = c("#99CC00","#c77cff"),
         htCol = c("#99CC00","#c77cff"),
         rankCol = c("#99CC00", "white", "#c77cff"))
  ggsave(paste0(i,'.pdf'),width = 6, height = 5)
  
}



############kegg_gsea#############
#读取从GSEA官网下载好的基因集文件；
kegg_gsea <- read.gmt("c2.cp.kegg.v2023.1.Hs.symbols.gmt")

#执行GSEA;
kegg_gsea_res <- GSEA(gene_list,TERM2GENE = kegg_gsea,pvalueCutoff = 1,verbose = F)
#预览分析结果；
head(kegg_gsea_res@result[1:10])
save(kegg_gsea_res,file = './gsea/kegg_gsea_res.Rdata')
#load('./gsea/go_gsea_res.Rdata')
write.csv(kegg_gsea_res@result,file = './gsea/kegg_gesa.csv')

#devtools::install_github("junjunlab/GseaVis")
library(GseaVis)
#自定义配色2：
#对感兴趣基因集的富集结果进行可视化；
kgset="GOBP_NEUTROPHIL_MEDIATED_KILLING_OF_SYMBIONT_CELL"
kgsetid<- go_gsea_res@result$ID
#获取位置编号；
index2 <- which(kgsetid==kgset)


gseaNb(object = go_gsea_res,
       geneSetID = go_gsea_res@result$ID[index2],
       subPlot = 3,
       addPval = T,
       pvalX = 0.95,
       pvalY = 0.8,
       curveCol = c("#99CC00","#c77cff"),
       htCol = c("#99CC00","#c77cff"),
       rankCol = c("#99CC00", "white", "#c77cff"))
ggsave('GOBP_NEUTROPHIL_MEDIATED_KILLING_OF_SYMBIONT_CELL.pdf',width = 6, height = 5)

gene_sets <- c('GOBP_NEUTROPHIL_MEDIATED_KILLING_OF_SYMBIONT_CELL',
               'GOBP_NEUTROPHIL_MEDIATED_CYTOTOXICITY',
               'GOBP_NEUTROPHIL_MEDIATED_IMMUNITY',
               'GOBP_NEUTROPHIL_CHEMOTAXIS',
               'GOBP_NEUTROPHIL_MIGRATION',
               'GOBP_NEUTROPHIL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
               'GOBP_REGULATION_OF_NEUTROPHIL_ACTIVATION',
               'GOBP_NEUTROPHIL_DEGRANULATION',
               'GOBP_REGULATION_OF_NEUTROPHIL_MIGRATION',
               'GOBP_NEUTROPHIL_EXTRAVASATION',
               'GOBP_REGULATION_OF_NEUTROPHIL_CHEMOTAXIS',
               'GOBP_NEUTROPHIL_HOMEOSTASIS',
               'GOBP_POSITIVE_REGULATION_OF_NEUTROPHIL_MIGRATION',
               'GOBP_NEUTROPHIL_DIFFERENTIATION')
for (i in gene_sets) {
  kgset= i
  kgsetid<- go_gsea_res@result$ID
  #获取位置编号；
  index2 <- which(kgsetid==kgset)
  gseaNb(object = go_gsea_res,
         geneSetID = go_gsea_res@result$ID[index2],
         subPlot = 3,
         addPval = T,
         pvalX = 0.95,
         pvalY = 0.8,
         curveCol = c("#99CC00","#c77cff"),
         htCol = c("#99CC00","#c77cff"),
         rankCol = c("#99CC00", "white", "#c77cff"))
  ggsave(paste0(i,'.pdf'),width = 6, height = 5)
  
}


###############gene_boxplot#####################
#加载数据
df <- read.table("all_exp_anno.txt",header = 1)
df <- df[, -c(1,17)]
# 转换成TPM


s100A9 <- df[which(df$gene_name == 'S100A9'),]
rownames(s100A9) <- s100A9$gene_name
s100A9 <- s100A9[,-15]
s100A9 <- t(s100A9)
group <- c(rep('R',7),rep('NR',7))
s100A9 <- as.data.frame(cbind(s100A9,group))
str(s100A9)
s100A9$S100A9 <- as.numeric(s100A9$S100A9)
s100A9$S100A9 <- log2(s100A9$S100A9)


library(ggplot2)
library(gghalves)
library(ggsignif)
##绘图
# ggplot(s100A9,aes(group,S100A9,fill=group))+
#   geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.8,color=NA)+
#   geom_boxplot(width=0.4,size=1.2,outlier.color =NA)+
#   geom_jitter(aes(fill=group),shape=21,size=2.5,width=0.2)+
#   geom_hline(yintercept = 0, linetype = 2, color = "red",linewidth=1)+
#   geom_hline(yintercept = 80, linetype = 2, color = "red",linewidth=1)+
#   geom_signif(comparisons = list(c("A","B"),  c("A","C"),  c("C","D")),
#               map_signif_level = T, test = t.test, y_position = c(100,120,130),
#               tip_length = c(0,0,0,0,0,0),size=1,color="black",textsize = 7)+
#   scale_y_continuous(limits = c(-20,140),breaks = c(0,40,80,120))+
#   theme_bw()+
#   theme(panel.grid = element_blank(), panel.border = element_rect(size = 1),
#         axis.text.x = element_text(color = "black", size = 13),
#         axis.text.y = element_text(color = "black",size = 13),legend.position = "none",
#         axis.ticks = element_line(color="black",linewidth = 1))+
#   labs(x=NULL,y=NULL)+
#   scale_fill_manual(values = c("#5cc3e8","#ffdb00","#79ceb8","#e95f5c"))


ggplot(s100A9,aes(group,S100A9,fill=group))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.8,color=NA)+
  geom_boxplot(width=0.4,size=1.2,outlier.color =NA)+
  geom_jitter(aes(fill=group),shape=21,size=2.5,width=0.2)+
  geom_hline(yintercept = 14, linetype = 2, color = "#79ceb8",linewidth=1)+
  geom_hline(yintercept = 18, linetype = 2, color = "#79ceb8",linewidth=1)+
  geom_signif(comparisons = list(c("NR","R")),
              map_signif_level = T, test = t.test, y_position = c(19),
              tip_length = c(0,0,0,0,0,0),size=1,color="black",textsize = 7)+
  scale_y_continuous(limits = c(13,20),breaks = c(13,15,17,19))+
  theme_bw()+
  theme(panel.grid = element_blank(), panel.border = element_rect(size = 1),
        axis.text.x = element_text(color = "black", size = 13),
        axis.text.y = element_text(color = "black",size = 13),legend.position = "none",
        axis.ticks = element_line(color="black",linewidth = 1))+
  labs(x=NULL,y=NULL)+
  ggtitle('S100A9',)+
  theme(plot.title = element_text(hjust = 0.5, size = 16))+
  scale_fill_manual(values = c("#5cc3e8","#ffdb00","#79ceb8","#e95f5c"))

ggsave('S100A9_boxplot.pdf',width = 6, height = 6)


###############循环gene_boxplot#####################
#加载数据
df <- read.table("all_exp_anno.txt",header = 1)
df <- df[, -c(1,17)]
# 转换成TPM
genes <- c('S100A8','S100A9','S100A12')
gene_log <- data.frame()

for (i in genes) {
  gene_exp <- df[which(df$gene_name == i),]
  gene_log <- rbind(gene_log,gene_exp)
}


rownames(gene_log) <- gene_log$gene_name
gene_log <- gene_log[,-15]
gene_log <- t(gene_log)
group <- c(rep('R',7),rep('NR',7))
gene_log <- as.data.frame(cbind(gene_log,group))
gene_log[,1:3] <- apply(gene_log[,1:3],2,as.numeric)
gene_log[,1:3] <- log2(gene_log[,1:3])

library(ggplot2)
library(gghalves)
library(ggsignif)
##绘图
# ggplot(s100A9,aes(group,S100A9,fill=group))+
#   geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.8,color=NA)+
#   geom_boxplot(width=0.4,size=1.2,outlier.color =NA)+
#   geom_jitter(aes(fill=group),shape=21,size=2.5,width=0.2)+
#   geom_hline(yintercept = 0, linetype = 2, color = "red",linewidth=1)+
#   geom_hline(yintercept = 80, linetype = 2, color = "red",linewidth=1)+
#   geom_signif(comparisons = list(c("A","B"),  c("A","C"),  c("C","D")),
#               map_signif_level = T, test = t.test, y_position = c(100,120,130),
#               tip_length = c(0,0,0,0,0,0),size=1,color="black",textsize = 7)+
#   scale_y_continuous(limits = c(-20,140),breaks = c(0,40,80,120))+
#   theme_bw()+
#   theme(panel.grid = element_blank(), panel.border = element_rect(size = 1),
#         axis.text.x = element_text(color = "black", size = 13),
#         axis.text.y = element_text(color = "black",size = 13),legend.position = "none",
#         axis.ticks = element_line(color="black",linewidth = 1))+
#   labs(x=NULL,y=NULL)+
#   scale_fill_manual(values = c("#5cc3e8","#ffdb00","#79ceb8","#e95f5c"))


# .data[[i]] #指代前面的数据
# .data[[i]] 是一种用于在ggplot2中处理变量的方式。在您的代码中，i 是您在genes列表中迭代的每个基因的名称。
# 
# .data[[i]] 的作用是让ggplot2在aes()函数中获取gene_log数据框中的特定变量（基因），而不需要引用数据框的名称。这在循环中特别有用，因为您可以根据i的值轻松选择不同的基因，而无需显式指定数据框的名称。
# 
# 这里是解释：
# .data：这是ggplot2的一种特殊方式来访问数据框中的变量。它允许您通过变量的名称来引用数据框中的列。
# [[i]]：这部分是从.data中选择特定变量的语法。[[i]] 表示从数据框中选择名称为i的列。
# 所以，.data[[i]] 的作用是选择gene_log数据框中的i列，其中i是您正在循环中迭代的每个基因的名称。这使您能够动态地选择不同的基因进行绘图，而无需手动更改数据框的名称。

for (i in genes) {
  p <- ggplot(gene_log, aes(group, .data[[i]], fill = group)) +
    geom_half_violin(position = position_nudge(x = 0.25), side = "r", width = 0.8, color = NA) +
    geom_boxplot(width = 0.4, size = 1.2, outlier.color = NA) +
    geom_jitter(aes(fill = group), shape = 21, size = 2.5, width = 0.2) +
    geom_hline(yintercept = 10, linetype = 2, color = "#79ceb8", linewidth = 1) +
    geom_hline(yintercept = 18, linetype = 2, color = "#79ceb8", linewidth = 1) +
    geom_signif(comparisons = list(c("NR", "R")),
                map_signif_level = T, test = t.test, y_position = c(19),
                tip_length = c(0, 0, 0, 0, 0, 0), size = 1, color = "black", textsize = 7) +
    scale_y_continuous(limits = c(10, 20), breaks = c(11, 14, 17, 20)) +
    theme_bw() +
    theme(panel.grid = element_blank(), panel.border = element_rect(size = 1),
          axis.text.x = element_text(color = "black", size = 13),
          axis.text.y = element_text(color = "black", size = 13), legend.position = "none",
          axis.ticks = element_line(color = "black", linewidth = 1)) +
    labs(x = NULL, y = NULL) +
    ggtitle(i) +  # 设置标题为基因名称
    theme(plot.title = element_text(hjust = 0.5, size = 16)) +
    scale_fill_manual(values = c("#5cc3e8","#ffdb00","#79ceb8","#e95f5c"))
  
  ggsave(p, file = paste0(i, '_boxplot.pdf'), width = 6, height = 6)
}


###################GSVA######################
library(GSVA)
library(clusterProfiler)
library(limma)
library(stringr)
library(ggplot2)

#gsva条形图
count <- read.csv('表达谱.csv',row.names = 1)
count <- avereps(count[,-c(1)], ID=count$anno)
gs <-read.gmt("c2.cp.kegg.v7.5.1.symbols.gmt")

#gs$term <- gsub('KEGG_','',gs$term)

gs.list <- gs %>% split(.$term) %>% lapply( "[[", 2)

gsva_es <- gsva(as.matrix(count), 
                gs.list,
                method="ssgsea",
                abs.ranking=T)

#增加样本来源属性

group <- c(rep('normal',10),rep('RA',13))




group_list <- data.frame(sample = colnames(count), 
                         group = group)
head(group_list)

design <- model.matrix(~ 0 + factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva_es)
design


contrast.matrix <- makeContrasts(RA-normal, levels = design)

fit <- lmFit(gsva_es, design)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)

x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")

head(x)

write.csv(x, "gsva_limma.csv", quote = F)

df <- data.frame(ID = rownames(x), score = x$t)

mean(df$score)
#按照score的值分组
library(tidyverse)
df_1 <- df %>%filter(abs(score)>1.5) 
median(df_1$score)
cutoff <- 2
df_1$group <- cut(df_1$score, 
                  breaks = c(-Inf, -cutoff, cutoff, Inf),
                  labels = c(1,2,3))

#按照score排序
sortdf <- df_1[order(df_1$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
sortdf$group <- factor(sortdf$group)

head(sortdf)
write.csv(sortdf,file = 'gsva_result.csv')

gsva_input <- data.table::fread('gsva_result.csv',)
gsva_input <- gsva_input[order(gsva_input$score),]
gsva_input$ID <- factor(gsva_input$ID, levels = gsva_input$ID)
gsva_input$group <- factor(gsva_input$group)

ggplot(sortdf, aes(ID, 
                   score, 
                   fill = group)) + 
  geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c('palegreen3', 
                               'snow3', 
                               'dodgerblue4'), 
                    guide = FALSE) + 
  geom_hline(yintercept = c(-cutoff,cutoff), 
             color="white",
             linetype = 2, #画虚线
             size = 0.3) + 
  geom_text(data = subset(sortdf, score < 0),
            aes(x=ID, 
                y= 0.1, 
                label= ID, 
                color = group),#bar跟坐标轴间留出间隙
            size = 3, #字的大小
            hjust = 0 ) +  #字的对齐方式
  geom_text(data = subset(sortdf, score > 0),
            aes(x=ID, 
                y= -0.1, 
                label=ID, 
                color = group),
            size = 3, 
            hjust = 1) +  
  scale_colour_manual(values = c("black","snow3","black"), guide = FALSE) +
  xlab("") +ylab("t value of GSVA score, NEC \n versus health")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 0.6)) + #边框粗细
  theme(axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank()) #去除y轴

ggsave(filename = "gsva.pdf",height = 8,width = 11)


###样本gsva
sortdf_1 <- sortdf %>%filter(abs(score)>2) 

pdf("Fig3f.pdf",height = 10,width = 10)
pheatmap::pheatmap(gsva_es[sortdf_1$ID,] ,
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                   cluster_cols = F,
                   cluster_rows = F,
                   scale = "row")
dev.off()

# 绘制差异分析热图
load(file = "data/GSVA_DEG_limma.Rda")
load(file = "data/GSVA_MsigDb_result.Rda")

DEG_sig <- x[x$P.Value<0.05 & abs(x$logFC) > 0.3,]

dat <- gsva_es[match(rownames(DEG_sig),rownames(gsva_es)),]
annotation_col <- data.frame(group_list)
rownames(annotation_col) <- colnames(count)
pheatmap::pheatmap(dat, 
                   width = 20, 
                   height = 11,
                   annotation_col = annotation_col,
                   show_colnames = F,
                   cluster_cols = FALSE,
                   filename = 'GSVA_heatmap.pdf')


##kegg_Metabolism
df <- read.csv('kegg_Metabolism.csv',row.names = 1)

library(clusterProfiler)

#######自制kegg代谢gmt#################
df
head(df)
library(tidyverse)
df_expanded <- separate_rows(df, genes, sep = ";")
Metabolism <- df_expanded
Metabolism$nm <-  gsub("- Homo sapiens \\(human\\)", "", Metabolism$nm)
Metabolism <- Metabolism[,-1]
colnames(Metabolism) <- c('term','gene')

##################GSVA#######################
library(GSVA)
library(clusterProfiler)
library(limma)
library(stringr)
library(ggplot2)

#gsva条形图
count <- read.table('all_exp_anno.txt')
count <- count[rowMeans(count[,2:15])>1,]
count <- count[!duplicated(count$gene_name),] %>% na.omit()
rownames(count) <- count$gene_name
count <- count[,-c(1,16,17)]

gs <-Metabolism

#gs$term <- gsub('KEGG_','',gs$term)

gs.list <- gs %>% split(.$term) %>% lapply( "[[", 2)

gsva_es <- gsva(as.matrix(count), 
                gs.list,
                method="ssgsea",
                abs.ranking=T)

#增加样本来源属性
# paste0('control',1:11),paste0('OA',1:11)
group <- c(rep('R',7),rep('NR',7))

group_list <- data.frame(sample = colnames(count), 
                         group = group)
head(group_list)

design <- model.matrix(~ 0 + factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva_es)
design


contrast.matrix <- makeContrasts(NR-R, levels = design)

fit <- lmFit(gsva_es, design)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)

x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")

head(x)
rownames(x) <- gsub(',','_',rownames(x))

write.csv(x, "gsva_limma.csv", quote = F,row.names = T)
#x<- x[x$P.Value<0.1,]

df <- data.frame(ID = rownames(x), score = x$t)

mean(df$score)
#按照score的值分组
library(tidyverse)
#df_1 <- df %>%filter(abs(score)>1.5) 
median(df_1$score)
cutoff <- 2
#cutoff <- 1.5
df_1$group <- cut(df_1$score, 
                  breaks = c(-Inf, -cutoff, cutoff, Inf),
                  labels = c(1,2,3))
table(df_1$group)
#按照score排序
sortdf <- df_1[order(df_1$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
sortdf$group <- factor(sortdf$group)

head(sortdf)
write.csv(sortdf,file = 'gsva_result.csv')

gsva_input <- data.table::fread('gsva_result.csv',)
gsva_input <- gsva_input[order(gsva_input$score),]
gsva_input$ID <- factor(gsva_input$ID, levels = gsva_input$ID)
gsva_input$group <- factor(gsva_input$group)

ggplot(sortdf, aes(ID, 
                   score, 
                   fill = group)) + 
  geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c("navy", 'snow3', "firebrick3"), 
                    guide = FALSE) + 
  geom_hline(yintercept = c(-cutoff,cutoff), 
             color="white",
             linetype = 2, #画虚线
             size = 0.3) + 
  geom_text(data = subset(sortdf, score < 0),
            aes(x=ID, 
                y= 0.1, 
                label= ID, 
                color = group),#bar跟坐标轴间留出间隙
            size = 3, #字的大小
            hjust = 0 ) +  #字的对齐方式
  geom_text(data = subset(sortdf, score > 0),
            aes(x=ID, 
                y= -0.1, 
                label=ID, 
                color = group),
            size = 3, 
            hjust = 1) +  
  scale_colour_manual(values = c("black","snow3","black"), guide = FALSE) +
  xlab("") +ylab("t value of GSVA score, \n NR versus R")+
  theme_classic() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  #theme(panel.border = element_rect(size = 0.6)) + #边框粗细
  theme(axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank()) #去除y轴

ggsave(filename = "gsva.pdf",height = 9.5,width = 10)


###样本gsva
sortdf_1 <- sortdf %>%filter(abs(score)>2) 

pdf("sample_gsva.pdf",height = 10,width = 10)
pheatmap::pheatmap(gsva_es[sortdf$ID,] ,
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                   cluster_cols = F,
                   cluster_rows = F,
                   scale = "row")
dev.off()

# 绘制差异分析热图
load(file = "data/GSVA_DEG_limma.Rda")
load(file = "data/GSVA_MsigDb_result.Rda")

DEG_sig <- x[x$P.Value<0.05 & abs(x$logFC) > 0.3,]

dat <- gsva_es[match(rownames(DEG_sig),rownames(gsva_es)),]
annotation_col <- data.frame(group_list)
rownames(annotation_col) <- colnames(count)
pheatmap::pheatmap(dat, 
                   width = 20, 
                   height = 11,
                   annotation_col = annotation_col,
                   show_colnames = F,
                   cluster_cols = FALSE,
                   filename = 'GSVA_heatmap.pdf')


####################免疫浸润####################
count <- read.table('all_exp_anno.txt')
count <- count %>% na.omit()
count <- count[rowMeans(count[,2:15] ) >1,]
count <- count[!duplicated(count$gene_name),]
rownames(count) <- count$gene_name
count <- count[,2:15]

library(IOBR)
library(tidyr)
tpm=count2tpm(count,idType = "SYMBOL")
tpm[1:4,1:4]
cibersort<-deconvo_tme(eset = tpm, method = "cibersort", arrays = FALSE, perm = 200 )

# cibersort<-deconvo_tme(mat, method = "cibersort",perm = 1000 )
write.table(cibersort,"CIBERSORT.txt",sep="\t",quote = F,row.names = F)
cibersort$ID=seq(1,ncol(tpm),1)
res<-cell_bar_plot(input = cibersort, title = "CIBERSORT Cell Fraction",coord_filp = F)

pdf(file="CIBERSORT.pdf",  width=10, height=7)
print(res)
dev.off()

library(ggsci)
library(randomcoloR)
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
immFile="CIBERSORT.txt"   
#cluFile="sample.txt"             
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value_CIBERSORT"]<0.05,]
data=as.matrix(immune[,1:(ncol(immune)-3)])
colnames(immune)=gsub("_CIBERSORT"," ",colnames(immune))
colnames(data)=gsub("_CIBERSORT"," ",colnames(data))
colnames(immune)=gsub("_"," ",colnames(immune))
colnames(data)=gsub("_"," ",colnames(data))
cluster=read.table('sample.txt', header=T, sep="\t", check.names=F, row.names=1)
colnames(cluster)="Type"
sameSample=intersect(row.names(data), row.names(cluster))
data=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])
data=data[order(data$Type),]
gaps=c(1, as.vector(cumsum(table(data$Type))))
xlabels=levels(factor(data$Type))

data=melt(data,id.vars=c("Type"))
colnames(data)=c("Type", "Immune", "Expression")

group=levels(factor(data$Type))
data$Type=factor(data$Type, levels=group)

# 调色板
# bioCol=pal_jco()(6)
# bioCol=bioCol[1:length(group)]
bioCol=c("#80B1D3", "#F9918C")




boxplot=ggboxplot(data, x="Immune", y="Expression",  fill="Type",
                  xlab="",
                  ylab="CIBERSORT Fraction",
                  legend.title="Typer", 
                  width=0.8,
                  palette=bioCol,add.params = list(size=0.1))
boxplot=boxplot+stat_compare_means(aes(group=Type),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")+theme_bw()+rotate_x_text(50)

pdf(file="immune.diff.pdf", width=9, height=4.5)
print(boxplot)
dev.off()


pdf(file="class_immune.diff.pdf", width=9, height=4.5)

print(boxplot+stat_compare_means(aes(group=Type),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")+theme_classic()+rotate_x_text(50))
dev.off()


#-------------------------#
# 复现 相关性分析 #
#-------------------------#

# 设置工作路径
workdir <- "G:/ExternalProject/SNOW/Figure9B相关性分析"; setwd(workdir)

# 加载R包
# install.packages("ggplot2")
# install.packages("ggpubr")
library(ggpubr)
library(ggplot2)

# 数据源自Data下载目录里的tpms_mrna.txt表达谱以及Figure8Abarplot目录下的cibersort_output.txt文件

# 加载表达谱
#expr <- read.table("aaaa.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

# 选取目标基因
target.gene <- "S100A9"
if(!is.element(target.gene,rownames(dat2))) { # 如果表达谱里不存在这个基因则报错
  stop("The gene you choose is not included in the expression profile!")
}

# 取出目标基因表达（表达谱对数化）
target.gene.expr <- as.numeric(log2(dat2[target.gene,] + 1))
names(target.gene.expr) <- colnames(dat2)

# 加载CIBERSORT免疫丰度数据
ciber_output <- read.table("CIBERSORT.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ciber.res <- ciber_output[,1:22] # 取前22列为细胞丰度数据
ciber.res <- ciber.res[,colSums(ciber.res) > 0] # 去除丰度全为0的细胞

# 循环生成相关性散点图
p.cutoff <- 0.05 # 相关性检验p值的阈值
outTab <- NULL
for (i in colnames(ciber.res)) {
  message(paste0("--analysis of ",i," done..."))
  # 构建数据库，包括基因表达以及对应细胞的丰度
  dat <- data.frame(gene.dat2 = as.numeric(target.gene.expr),
                    cell.score = as.numeric(ciber.res[names(target.gene.expr),i]),
                    stringsAsFactors = F)
  
  cor.res <- cor.test(dat$gene.dat2,dat$cell.score, method = "pearson") # 相关性分析
  outTab <- rbind.data.frame(outTab,
                             data.frame(gene = target.gene,
                                        cell = i,
                                        rho = cor.res$estimate,
                                        pval = cor.res$p.value,
                                        stringsAsFactors = F),
                             stringsAsFactors = F)
  
  if(cor.res$p.value < p.cutoff) { # 若检验p值小于阈值则出图
    sp <- ggscatter(dat, 
                    x = "gene.dat2", # x为基因表达
                    y = "cell.score", # y为细胞丰度
                    xlab = paste0(target.gene," expression level"),
                    ylab = i,
                    size = 0.8,
                    color = "blue", # 散点颜色
                    add = "reg.line",  # 添加回归线
                    add.params = list(color = "red", fill = "grey40"), # 设置回归线以及回归置信区间的颜色
                    conf.int = TRUE) + 
      stat_cor(method = "pearson", # 相关性分析方法
               label.x = min(dat$gene.dat2), # 结果标注的x位置
               label.y = max(dat$cell.score)) # 结果标注的y位置
    ggsave(file = paste0("correlation scatter plot between expression of ", target.gene," and ", i, ".pdf"), width = 4, height = 4)
  }
}
write.table(outTab,paste0("correlation results between ",target.gene," and cibersort cells.txt"),sep = "\t",row.names = F,col.names = T,quote = F)


##################2.ssGSEA############################
ssGSEAScore_by_genes<-function(gene.exp,genes){
  gs=GSEABase::GeneSet(setName='GeneSet', setIdentifier=paste0("101"),geneIds=unique(genes),GSEABase::SymbolIdentifier()) 
  
  gsc <- GSEABase::GeneSetCollection(list(gs))
  fl <- tempfile()
  GSEABase::toGmt(gsc, fl)
  cgeneset=GSEABase::getGmt(fl)
  ssGSEA.geneset <- GSVA::gsva(as.matrix(gene.exp), cgeneset,method='ssgsea',
                               min.sz=1, max.sz=5000, verbose=TRUE)
  return(ssGSEA.geneset)
}
pathway.score<-function(exp,gene){
  pathway_score<-data.frame()
  for (i in unique(gene[,2])){
    gene_set=gene[gene[,2]==i,1]
    score=ssGSEAScore_by_genes(exp,gene_set)
    rownames(score)=i
    pathway_score=rbind.data.frame(pathway_score,score)
  }
  return(t(pathway_score))
}
ssgase.imm<-read.delim('./immune/ssgsea/28_immune_pmid_28052254.txt',sep='\t',header = T)
ssgase.imm=ssgase.imm[,c(1,2)]
colnames(ssgase.imm)[1]='gene'
ssgsea<-pathway.score(exp = exp,gene =ssgase.imm )

write.table(ssgsea,'./immune/ssgsea/ssgsea.txt',quote = F,sep = '\t',row.names = T)
ssgsea1=reshape2::melt(ssgsea)
subtype=read.table('sample.txt')
colnames(subtype) <- c("ID","group")
ssgsea1=merge(subtype[,c("ID","group")],
              ssgsea1,by.x='ID',by.y='Var1')
library(ggpubr)

table(ssgsea1$Var2)
figb<-ggboxplot(ssgsea1, x='Var2', y='value', 
                fill = "group", color = "black",
                # palette = ggsci::pal_aaas()(9)[4:6], 
                ylab="Score",xlab='',
                add = "boxplot")+ 
  stat_compare_means(aes(group=group),method = 'kruskal.test',
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),
                     label = "p.signif", size=4) +  # 调整标签的大小)
  theme_classic2() +  # 使用最小化主题
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 7),
        legend.title = element_blank(),  # 隐藏图例标题
        legend.position = "top",
        axis.title.x = element_text(size = 12),  # 调整x轴标签的大小
        axis.title.y = element_text(size = 12),  # 调整y轴标签的大小
        plot.title = element_text(size = 16))+     # 调整图形标题的大小
   # 添加图形标题
  scale_fill_manual(values = c("#223D6C","#D20A13"),  # 设置填充色
                    labels = c("NR", "R"))  # 设置图例标

figb
pdf(file="class_ssgsea.diff.pdf", width=9, height=4.5)
# 显示图形
print(figb)
dev.off()

# 调整图形的大小
options(repr.plot.width=8, repr.plot.height=5)  # 调整图形的宽度和高度


