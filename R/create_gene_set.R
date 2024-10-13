library(jsonlite)
library(clusterProfiler)
library(org.Hs.eg.db)

# 设置文件夹路径
folder_path <- "../data/02_Gene_Sets/keyword_iron"

# 读取文件夹中所有JSON文件
file_paths <- list.files(folder_path, pattern="*.json", full.names=TRUE)

# 初始化一个空的数据框用于存储所有结果
all_data <- data.frame(gs_name=character(), ENTREZID=character(), stringsAsFactors=FALSE)

# 遍历每个文件
for (file_path in file_paths) {
  # 读取JSON文件
  json_data <- fromJSON(file_path)
  geneset_name <- names(json_data)
  gene_symbols <- json_data[[1]]$geneSymbols
  
  # 转换基因符号到Entrez ID
  gene_info <- bitr(gene_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  # 生成部分结果数据框
  partial_data <- data.frame(gs_name = rep(geneset_name, nrow(gene_info)), ENTREZID = gene_info$ENTREZID, stringsAsFactors = FALSE)
  
  # 合并到全局数据框
  all_data <- rbind(all_data, partial_data)
}