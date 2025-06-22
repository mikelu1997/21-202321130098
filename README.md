# 加载软件包
library(clusterProfiler)
library(org.Hs.eg.db)
library(DO.db)
BiocManager::install("clusterProfiler")
# 读取基因列表
input <- read.table("data/clusterProfiler_input.txt",sep="\t",header=T,check.names=F)#以tab（制表符）分隔，需要表头，不检查名字
input <- input[is.na(input[,"ENTREZID"])==F,]#is.na判断是否为NA，是返回T，否返回F，则这里返回的是非NA的，也就是把NA的去掉
gene <- input$ENTREZID# input$ENTREZID叫gene

# 如果你的基因列表是gene symbol，可以使用如下函数转换，将SYMBOL转换为基因id
genee1 <- bitr(GSE89632_DEG_phua2$symbol,
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)# 物种是人类
head(genee)
load( "D:/Program Files/基因集泛癌加预后模型分析/genee.Rda")
########## GO富集分析 ########## 
ego <- enrichGO(gene = genee$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont="All",
                pvalueCutoff =0.01, 
                qvalueCutoff =0.05,
                readable = TRUE)
head(ego)
# 将GO富集结果写出到文件
ego1<-as.data.frame(ego)
kk1<-as.data.frame(kk)
write.table(ego,file="D:/Program Files/clusterProfiler_GO_result.txt",sep="\t",
            quote=F,row.names = F)

########## KEGG富集分析 ########## 
install.packages("R.utils")
remotes::install_github("YuLab-SMU/createKEGGdb")
#创建自己的物种的包create_kegg_db，会自动创建名称为KEGG.db_1.0.tar,gz的包。物种名称的简写，在
createKEGGdb::create_kegg_db('zma')
source("http://bioconductor.org/biocLite.R")
BiocManager::install("KEGG.db")
biocLite("KEGG.db")
R.utils::setOption("clusterProfiler.download.method",'auto') 
R.utils::setOption("clusterProfiler.download.method",'auto') 
options(clusterProfiler.download.method = "wininet")
options(clusterProfiler.download.method = "auto")
install.packages("R.utils")

kk <- enrichKEGG(gene = genee$ENTREZID[1:1000],
                 organism = "hsa",
                 keyType  = 'kegg',
                 pvalueCutoff = 0.2,
                 qvalueCutoff =0.2,
                 use_internal_data = F)
