mRNA<-read.csv("All_mRNA_FPKM.csv",header=T,row.names=1)
gc<-c("CD28","CD3D","CD8A","LCK","ZAP70",
                     "GATA3","EOMES","IL23A",
                     "CXCL8","CXCR1","CXCL1","CXCR2",
                     "PTGS2","IL1R2","IL1R1","MMP8","MMP9","S100A12","S100A8",
                     "TLR4","TLR6","IL1B","TNF","NEDD4L","CDC34","UBE2E3","FOXO3","GABARAPL2","PINK1")

exp<-log2(mRNA+1)
bar_mat<-t(exp)
bar_mat<-bar_mat[,gc]
anno<-read.table("sample_index.txt",header=T,row.names=1)
anno$type2<-anno$type
bar_mat<-bar_mat[rownames(anno),]
bar_mat<-as.data.frame(bar_mat)
bar_mat$sam=anno$type

bar_mat$sam<-factor(bar_mat$sam,levels=c("Asymptomatic","Mild","Severe","Critical"))
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
c=getPalette(9)
library(ggpubr)
library(ggplot2)
plist2<-list()
co<-c("#5CB85C","#337AB7","#F0AD4E","#D9534F")
for (i in 1:length(gc)){
  bar_tmp<-bar_mat[,c(gc[i],"sam")]
  colnames(bar_tmp)<-c("Expression","sam")
  bar_tmp$color<-c(rep(x="1",times=64),rep(x="2",times=64),rep(x="3",times=34),rep(x="4",times=16))
  my_comparisons1 <- list(c("Asymptomatic", "Mild"))
  my_comparisons2 <- list(c("Asymptomatic", "Severe"))
  my_comparisons3 <- list(c("Asymptomatic", "Critical"))
  my_comparisons4 <- list(c("Mild", "Severe"))
  my_comparisons5 <- list(c("Mild", "Critical"))
  my_comparisons6 <- list(c("Severe", "Critical"))
  pb1<-ggboxplot(bar_tmp,x="sam",y="Expression",color="color",fill=NULL,add = "jitter",bxp.errorbar.width = 0.6,width = 0.4,size=0.01,font.label = list(size=30), palette = c(co[1],co[2],co[3],co[4]))+theme(panel.background =element_blank())
  pb1<-pb1+theme(axis.line=element_line(colour="black"))+theme(axis.title.x = element_blank())
  #pb1<-pb1+scale_fill_manual(name="",labels=c("Ctrl","COV"),values=c(c[2],c[1]))
  pb1<-pb1+theme(axis.title.y = element_blank())+theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))
  pb1<-pb1+theme(axis.text.y = element_text(size = 15))+ggtitle(gc[i])+theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))
  pb1<-pb1+theme(legend.position = "NA")
  pb1<-pb1+stat_compare_means(method="t.test",hide.ns = F,comparisons =c(my_comparisons1,my_comparisons2,my_comparisons3,my_comparisons4,my_comparisons5,my_comparisons6),label="p.signif")
  plist2[[i]]<-pb1
} 
pall<-plot_grid(plist2[[1]],plist2[[2]],plist2[[3]],
                plist2[[4]],plist2[[5]],plist2[[6]],
                plist2[[7]],plist2[[8]],plist2[[9]],
                plist2[[10]],plist2[[11]],plist2[[12]],plist2[[13]],plist2[[14]],
                plist2[[15]],plist2[[16]],plist2[[17]],plist2[[18]],
                plist2[[19]],plist2[[20]],plist2[[21]],
                plist2[[22]],plist2[[23]],plist2[[24]],
                plist2[[25]],plist2[[26]],plist2[[27]],
                plist2[[28]],plist2[[29]],ncol=6)
pall