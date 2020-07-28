DX<-read.table("metabolite_tryptophan_update_clean.txt",sep="\t",header=T,row.names=1)
A<-colnames(DX)[grepl(colnames(DX),pattern = "Asymptomatic")]
B<-colnames(DX)[grepl(colnames(DX),pattern = "Mild")]
C<-colnames(DX)[grepl(colnames(DX),pattern = "Severe")]
D<-colnames(DX)[grepl(colnames(DX),pattern = "Critical")]
E<-colnames(DX)[grepl(colnames(DX),pattern = "Death")]
D<-c(D,E)
order<-c(A,B,C,D)
order_anno<-c(rep(x="Asymptomatic",times=length(A)),rep(x="Mild",times=length(B)),rep(x="Severe",times=length(C)),rep(x="Critical",times=length(D)))

gene<-rownames(DX)

exp<-log2(DX+1)
exp<-exp[,order]
bar_mat<-as.data.frame(t(exp))
bar_mat$sam=order_anno


bar_mat$sam<-factor(bar_mat$sam,levels=c("Asymptomatic","Mild","Severe","Critical"))
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
c=getPalette(9)
library(ggpubr)
library(ggplot2)
plist2<-list()
for (i in 1:length(gene)){
  bar_tmp<-bar_mat[,c(gene[i],"sam")]
  colnames(bar_tmp)<-c("Expression","sam")
  bar_tmp$color<-c(rep(x="1",times=53),rep(x="2",times=54),rep(x="3",times=33),rep(x="4",times=21))
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
  pb1<-pb1+theme(axis.text.y = element_text(size = 15))+ggtitle(gene[i])+theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))
  pb1<-pb1+theme(legend.position = "NA")
  pb1<-pb1+stat_compare_means(method="t.test",hide.ns = F,comparisons =c(my_comparisons1,my_comparisons2,my_comparisons3,my_comparisons4,my_comparisons5,my_comparisons6),label="p.signif")
  plist2[[i]]<-pb1
} 
library(cowplot)

pdf("metabolites_new.pdf",width=25,height=16)
pall<-plot_grid(plist2[[1]],plist2[[2]],plist2[[3]],
                plist2[[4]],plist2[[5]],plist2[[6]],
                plist2[[7]],plist2[[8]],plist2[[9]],
                plist2[[10]],plist2[[11]],plist2[[12]],
                plist2[[13]],plist2[[14]],plist2[[15]],
                plist2[[16]],ncol=6)
pall
dev.off()