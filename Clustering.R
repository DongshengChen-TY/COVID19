c <- 7 #cluster number
m <- 0 #fuzzifier prevent clustering of random data
x <- 4 #row number
y <- 3 #column number
library("reshape2")
library("data.table")
library("Mfuzz")
library("dplyr")
library("readr")
library("magrittr")


raw_sub <- read_delim("clean_protein_meta_lipid.txt",delim="\t")

output_prefix <- "Asymptomatic.Mild.Severe.Critical"
groups_name <- c("Asymptomatic","Mild","Severe","Critical") 

rawColName <- raw_sub[1,]
colnames(raw_sub) <- raw_sub[1,]
raw_sub <- raw_sub[-1,]
timeIndex <- c()
for (i in groups_name) { 
  index <- grep(paste("^",i,sep=""),rawColName)
  timeIndex <- c(timeIndex,index)
}

raw_sub <- raw_sub[,c(1,timeIndex)]

raw_melt <- melt(raw_sub,id.vars = "Lable", variable.name="class", value.name="Intensity")
raw_melt$Intensity <- as.numeric(raw_melt$Intensity)
setDT(raw_melt)
raw_melt[,Condition:=unlist(strsplit(as.character(class),"_"))[1],.(Lable,class)]
raw_melt$Condition <- factor(raw_melt$Condition,levels=groups_name)
data <- dcast(raw_melt,Lable~Condition,fun=mean,na.rm=T,value.var="Intensity")
write.table(data,file=paste(output_prefix,".average.xls",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
data <- table2eset(paste(output_prefix,".average.xls",sep=""))  #input matrix
data.r <- filter.NA(data,thres=0.5)
data.m <- fill.NA(data.r,mode="knn",k=10)
data.m <- filter.NA(data.m,thres=0)
data.f <- filter.std(data.m,min.std=0,visu=F)
data.s <- standardise(data.f)
if (m==0) {
  m=mestimate(data.s)
}
m
cl <- mfuzz(data.s,c=c,m=m)
pdf(paste(output_prefix,".mfuzz.pdf",sep=""),width=10,height=10)
mfuzz.plot2(data.s,cl=cl,mfrow=c(x,y),time.labels=groups_name,x11=FALSE)
dev.off()
pdf(paste(output_prefix,".mfuzz.split.pdf",sep=""),width=10,height=10)
mfuzz.plot2(data.s,cl=cl,mfrow=c(1,1),time.labels=groups_name,x11=FALSE)
dev.off()
cluster <- as.data.frame(cl$cluster)
colnames(cluster) <- c("cluster")
membership <- as.data.frame(cl$membership)
membership_name <- c()
for (i in 1:c) {
  tmp <- paste("cluster",i,"_membership",sep="")
  membership_name <- c(membership_name,tmp)
}
colnames(membership) <- membership_name
#expStandard <- exprs(data.s)
write.table(cluster,file=paste(output_prefix,".cluster.xls",sep=""),col.names=NA,quote=F,sep="\t")
write.table(membership,file=paste(output_prefix,".membership.xls",sep=""),col.names=NA,quote=F,sep="\t")
system(paste("convert -density 300 -resize 40% ",output_prefix,".mfuzz.pdf ",output_prefix,".mfuzz.png",sep=""))
system(paste("convert -density 300 -resize 40% ",output_prefix,".mfuzz.split.pdf ",output_prefix,".mfuzz.split.png",sep=""))

