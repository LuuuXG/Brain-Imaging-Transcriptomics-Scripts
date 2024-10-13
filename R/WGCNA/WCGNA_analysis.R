#site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
#install.packages(c("WGCNA", "stringr", "reshape2"), repos=site)

#BiocManager::install(c("AnnotationDbi", "impute", "preprocessCore"))

library(WGCNA)
library(reshape2)
library(stringr)

# 避免转换为因子
options(stringsAsFactors = FALSE)
# 打开多线程
enableWGCNAThreads()

setwd('E:/Codes/Imaging-Transcriptomics/R')

dataExpr <- gene_expression_data_raw
dataTraits <- imaging_data_raw
suffix <- 'DK_all'

# TOMType: 是否考虑相关性的方向
TOMType = "signed" # one of "none", "unsigned", "signed", "signed Nowick", 
                   # "unsigned 2", "signed 2" and "signed Nowick 2"
# corType：相关性计算方法（连续数据选择pearson，离散数据选择bicor）
corType = "pearson" # "pearson" and "bicor"

# 与corType选择相关的选项
# 简单来说，pearson相关对应cor、maxPOutliers无意义、robustY为TRUE
corFnc = ifelse(corType=="pearson", cor, bicor) # ?

# 对二元变量，如样本性状信息计算相关性时，
# 或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05) # 这个参数与排除离群值有关

robustY = ifelse(corType=="pearson",T,F)

## 筛选中位绝对偏差前75%的基因，至少MAD大于0.01
## 筛选后会降低运算量，也会失去部分信息
## 也可不做筛选，使MAD大于0即可
m.mad <- apply(dataExpr, 2, mad)
dataExpr <- dataExpr[, which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2], 0.01))]

## 检测缺失值
gsg = goodSamplesGenes(dataExpr, verbose = 3)

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average")
sampleTree_plot_path <- file.path(WGCNA_output_folder, '01_sampleTree.jpeg')
#plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
jpeg(filename = sampleTree_plot_path, width = 1600, height = 1200)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.main = 2, cex.lab = 2, cex.axis = 2, cex = 1.2)
dev.off()

powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=TOMType, verbose=5)

soft_threshold_plot_path <- file.path(WGCNA_output_folder, '02_soft_threshold.jpeg')
jpeg(filename = soft_threshold_plot_path, width = 1600, height = 800)
par(mfrow = c(1,2), mar = c(5, 7, 4, 2) + 0.1)

# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"),
     cex.axis = 2, cex.lab = 2, cex.main = 2)
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=2,col="red")
# 筛选标准。R-square=0.85
abline(h=0.85,col="red")

# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"),
     cex.axis = 2, cex.lab = 2, cex.main = 2)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=2, col="red")
dev.off()

power = sft$powerEstimate

# 无向网络在power小于15或有向网络power小于30内，没有一个power值可以使
# 无标度网络图谱结构R^2达到0.8，平均连接度较高如在100以上，可能是由于
# 部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对
# 表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。
# 如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))       
                 )
  )
}

##一步法网络构建：One-step network construction and module detection##
# power: 上一步计算的软阈值
# maxBlockSize: 计算机能处理的最大模块的基因数量 (默认5000)；
#  4G内存电脑可处理8000-10000个，16G内存电脑可以处理2万个，32G内存电脑可
#  以处理3万个
#  计算资源允许的情况下最好放在一个block里面。
# corType: pearson or bicor
# numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
# saveTOMs：最耗费时间的计算，存储起来，供后续使用
# mergeCutHeight: 合并模块的阈值，越大模块越少
net = blockwiseModules(dataExpr, power = 20, maxBlockSize = nGenes,
                       TOMType = TOMType, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOM=T,
                       saveTOMFileBase = paste0(WGCNA_output_folder, '/', suffix, "_TOM"),
                       verbose = 3)

# 根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
# **0 (grey)**表示**未**分入任何模块的基因。 
table(net$colors)

## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)

MEs = net$MEs
geneTree = net$dendrograms[[1]]

save(MEs, moduleLabels, moduleColors, geneTree,
     file = paste0(WGCNA_output_folder, '/', suffix, "_auto.RData"))

load(paste0(WGCNA_output_folder, '/', suffix, "_auto.RData"))

# Plot the dendrogram and the module colors underneath
# 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间
cluster_plot_path <- file.path(WGCNA_output_folder, '03_cluster_plot.jpeg')
jpeg(filename = cluster_plot_path, width = 1200, height = 800)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    cex.colorLabels = 1, cex.main = 2, cex.lab = 2)
dev.off()

### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
MEs_path <- file.path(WGCNA_output_folder, 'MEs.csv')
write.csv(MEs_col, MEs_path)

# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
eigengene_adjacency_plot_path <- file.path(WGCNA_output_folder, '04_eigengene_adjacency_plot.jpeg')
jpeg(filename = eigengene_adjacency_plot_path, width = 1200, height = 800)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()

#load(net$TOMFiles[1], verbose=T)
#TOM <- as.matrix(TOM)

# Module-trait associations
# 前面已经计算了nGenes和nSamples
MEs0 = moduleEigengenes(dataExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, dataTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

module_trait_association_plot_path <- file.path(WGCNA_output_folder, 
                                                '05_module-trait_association_plot2.jpeg')
jpeg(filename = module_trait_association_plot_path, width = 2400, height = 3000, res = 450)
#sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(4, 8, 4, 4))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(dataTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# GS and MM
# GS：单个基因与临床数据的相关性（绝对值）
# MM：单个基因与所属模块ME值的相关性

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(dataExpr, dataTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(dataTraits), sep="")
names(GSPvalue) = paste("p.GS.", names(dataTraits), sep="")

# identifying genes with high GS and MM
module = 'red'
column = match(module, modNames)
moduleGenes = moduleColors==module

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for imaging data",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# summary output
#names(dataExpr)
#names(dataExpr)[moduleColors=="red"]
geneInfo0 = data.frame(geneSymbol = names(dataExpr),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

modOrder = order(-abs(cor(MEs, dataTraits, use = "p")))

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Measures));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "geneInfo.csv")
