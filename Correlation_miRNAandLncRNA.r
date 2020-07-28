library(corrgram)
data2 <- read.table('miRNA_gene_togather.txt', row.names = 1, sep = "\t", header = T,check.names = F,stringsAsFactors = F)
correl <- matrix(0,25*25,3)
colnames(correl) = c("gene_miRNA","correlation","p")

#spearman correlation
for(i in 1:25){
  for(j in 1:25){
    x <- cor.test(data2[,i],data2[,(j+25)], conf.level = 0.95, method = "spearman")
    correl[(25*(i-1)+j),3] = x$p.value
    correl[(25*(i-1)+j),2] = x$estimate
    correl[(25*(i-1)+j),1] = paste(colnames(data2)[i],colnames(data2)[j+25], sep='_')
    }
}
write.table(correl, "miRNA_gene_correlation(25).txt", quote = F, sep = "\t", row.names = F)

data3 <- read.table('lncRNA_gene_togather.txt', row.names = 1, sep = "\t", header = T,check.names = F,stringsAsFactors = F)
correl <- matrix(0,25*2484,3)
colnames(correl) = c("gene_lncRNA","correlation","p")

#spearman correlation
for(i in 1:25){
  for(j in 1:2484){
    x <- cor.test(data3[,i],data3[,(j+25)], conf.level = 0.95, method = "spearman")
    correl[(25*(i-1)+j),3] = x$p.value
    correl[(25*(i-1)+j),2] = x$estimate
    correl[(25*(i-1)+j),1] = paste(colnames(data3)[i],colnames(data3)[j+25], sep='_')
  }
}
write.table(correl, "lncRNA_gene_correlation(NCBI).txt", quote = F, sep = "\t", row.names = F)

library(corrgram)
data4 <- read.table('miRNA_gene_togather(2miRNA2gene).txt', row.names = 1, sep = "\t", header = T,check.names = F,stringsAsFactors = F)
correl <- matrix(0,2*2,3)
colnames(correl) = c("gene_miRNA","correlation","p")

#spearman correlation系数
for(i in 1:2){
  for(j in 1:2){
    x <- cor.test(data4[,i],data4[,(j+2)], conf.level = 0.95, method = "spearman")
    correl[(2*(i-1)+j),3] = x$p.value
    correl[(2*(i-1)+j),2] = x$estimate
    correl[(2*(i-1)+j),1] = paste(colnames(data4)[i],colnames(data4)[j+2], sep='_')
  }
}
write.table(correl, "miRNA_gene_correlation(2miRNA2gene).txt", quote = F, sep = "\t", row.names = F)


library(corrgram)
data5 <- read.table('miRNA_gene_togather(inflammation factors).txt', row.names = 1, sep = "\t", header = T,check.names = F,stringsAsFactors = F)
correl <- matrix(0,54*240,3)
colnames(correl) = c("gene_miRNA","correlation","p")

#spearman correlation系数
for(i in 1:54){
  for(j in 1:240){
    x <- cor.test(data5[,i],data5[,(j+54)], conf.level = 0.95, method = "spearman")
    correl[(54*(i-1)+j),3] = x$p.value
    correl[(54*(i-1)+j),2] = x$estimate
    correl[(54*(i-1)+j),1] = paste(colnames(data5)[i],colnames(data5)[j+54], sep='_')
  }
}
write.table(correl, "miRNA_gene_correlation(inflammation factors).txt", quote = F, sep = "\t", row.names = F)
