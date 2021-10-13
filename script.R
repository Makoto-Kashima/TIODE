# Load sample attribute
library(openxlsx)
at.mouse = read.xlsx("sample_info.xlsx",sheet = 1)
at.mouse = at.mouse[!is.na(at.mouse$name),]
tmp = unlist(lapply(strsplit(at.mouse$name, " "),FUN = function(x){return(x[1])}))
at.mouse = cbind(at.mouse, data.frame("stage" = tmp))
at.mouse$stage = factor(at.mouse$stage,levels = c("E7.5","E8.5","E9.5","E10.5","E11.5","E12.5","E13.5"))

# Load count data
list = list.files("salmon")
list.mouse = list[grep("Kuro", list)]
rawcnt.mouse = NULL
for(file in list.mouse){
  file = sprintf("salmon/%s/quant.sf",file)
  tmp = read.table(file, header = T, stringsAsFactors = F)
  rawcnt.mouse = cbind(rawcnt.mouse, tmp$NumReads)
}
rownames(rawcnt.mouse) = tmp$Name
colnames(rawcnt.mouse) = 1:96
rawcnt.mouse[,21:22] = rawcnt.mouse[,95:96] # mis-applied sample
rawcnt.mouse = rawcnt.mouse[,at.mouse$index]
save(rawcnt.mouse, file = "rawcnt-original")

# Marge count data of each gene
transcript2gene.mouse = read.csv("/mnt/hdd/hdd/data/species/M.musculus/mm9/ensembl/Mus_musculus.GRCm38.101.transcript2gene",row.names = 1)
transcript2gene.mouse = transcript2gene.mouse[is.element(rownames(transcript2gene.mouse), rownames(rawcnt.mouse)),]
all.gene = unique(transcript2gene.mouse$GeneName)
rawcnt.gn = NULL
for(gn in all.gene){
  tmp = rawcnt.mouse[rownames(transcript2gene.mouse)[transcript2gene.mouse$GeneName==gn],]
  if(is.matrix(tmp)){
    tmp = apply(tmp, 2, sum)
  }
  rawcnt.gn = rbind(rawcnt.gn, tmp)
}
colnames(rawcnt.gn) = colnames(rawcnt.mouse)
rownames(rawcnt.gn) = all.gene
save(rawcnt.gn, file = "rawcnt-gn")

# PCA with Seurat
library(Seurat)
seurat.mouse <- CreateSeuratObject(rawcnt.gn, project = "mouse", min.cells = 1)
seurat.mouse@meta.data$stage <- at.mouse$stage
seurat.mouse <- NormalizeData(seurat.mouse)
seurat.mouse <- ScaleData(seurat.mouse, do.scale = F, do.center = T)
seurat.mouse <- FindVariableFeatures(seurat.mouse, verbose = F, selection.method = "mvp")
seurat.mouse = RunPCA(seurat.mouse,verbose = F)
seurat.mouse = FindNeighbors(seurat.mouse,verbose = F)
seurat.mouse = FindClusters(seurat.mouse,verbose = T,resolution = 2.5)
DimPlot(seurat.mouse, reduction = "pca", group.by = c("stage","ident"))
markers = FindAllMarkers(seurat.mouse,logfc.threshold = 0.5)
pca = seurat.mouse@reductions$pca
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / sum(eigValues) # PC1:64.7%, PC2:21.2%

DoHeatmap(seurat.mouse, group.by = "stage",features = markers$gene)
DoHeatmap(seurat.mouse,features = markers$gene)


# Pseudotime analysis
sce = as.SingleCellExperiment(seurat.mouse)
library(slingshot)
library(SingleCellExperiment)
sce = slingshot(sce,clusterLabels = "stage",reducedDim = "PCA",start.clus = "E7.5", end.clus = "E13.5")
library(ggplot2)
tmp = SlingshotDataSet(sce)
pdf("Fig1ab.pdf")
g1 = DimPlot(seurat.mouse, reduction = "pca", group.by = c("stage"),pt.size = 2)+
  ylim(c(8,-7))+
  NoLegend()+
  geom_line(data = data.frame("x"= tmp@curves$curve1$s[,1], "y" = tmp@curves$curve1$s[,2]), aes(x = x, y = y))
g2 = ggplot(data.frame("stage" = sce$stage, "Pseudotime" = sce$slingPseudotime_1), aes(y = stage, x  = Pseudotime, color = stage))+
  geom_jitter(height = 0, height = 0.1,size = 2)+
  theme_bw()+
  NoGrid()+
  theme(legend.position = "none")
print(g1)
print(g2)
dev.off()

scode.1 =  as.matrix(seurat.mouse@assays$RNA@data)
scode.2 = data.frame("id" = 1:ncol(scode.1), "pseudotime" = sce$slingPseudotime_1)
save(scode.1, scode.2, file = "scode_input")

# Trendy for Fig.1c,d

SSR = NULL
for(gn in rownames(seurat.mouse@assays$RNA@data)){
  data = data.frame("ex" = scode.1[gn,], "pseudotime" = scode.2$pseudotime, "stage" = seurat.mouse$stage)
  stage = smooth.spline(as.numeric(data$stage),data$ex,all.knots = T, lambda = 0.001)
  pt = smooth.spline(data$pseudotime,data$ex,all.knots = T,lambda = 0.001)
  data.fit.stage = data.frame(x = as.numeric(data$stage), y = predict(stage,as.numeric(data$stage))$y)
  data.fit.pt = data.frame(x = data$pseudotime, y = predict(pt,data$pseudotime)$y)
  stage = sum(sqrt((data$ex - data.fit.stage$y )^2))
  pt = sum(sqrt((data$ex - data.fit.pt$y)^2))
  SSR = rbind(SSR,data.frame("stage" = stage,"pseudotime" = pt))
  # data.fit.stage = data.frame(x = seq(1,7,length.out = 1000), y = predict(stage,seq(1,7,length.out = 1000))$y)
  # data.fit.pt = data.frame(x = seq(range(data$pseudotime)[1],range(data$pseudotime)[2],length.out = 1000), y = predict(pt,seq(range(data$pseudotime)[1],range(data$pseudotime)[2],length.out = 1000))$y)
  # g1 = ggplot(data, aes(x = as.numeric(stage), y = ex, color = stage))+
  #   geom_line(data = data.fit.stage, aes(x = as.numeric(x), y = y), color = "black")+
  #   geom_point()+
  #   theme_bw()+
  #   theme(legend.position = "none")
  # g2 = ggplot(data, aes(x = pseudotime, y = ex, color = stage))+
  #   geom_line(data = data.fit.pt, aes(x = x, y = y), color = "black")+
  #   geom_point()+
  #   theme_bw()+
  #   theme(legend.position = "none")
  # g = grid.arrange(g1,g2,ncol =2, top = gn)
  # print(g)
}
pdf("Fig1d.pdf")
ggplot(SSR, aes(x = stage, y = pseudotime))+
  geom_point(alpha = 1)+
  geom_abline(intercept = 0,slope = 1,color = "red")+
  theme_bw()
dev.off()
rownames(SSR) = rownames(seurat.mouse@assays$RNA@data)
100-mean(SSR$pseudotime/SSR$stage*100)
100-mean(SSR[seurat.mouse@assays$RNA@var.features,]$pseudotime/SSR[seurat.mouse@assays$RNA@var.features,]$stage*100)
pdf("fig1c.pdf")
gn = "Fabp7"
  data = data.frame("ex" = scode.1[gn,], "pseudotime" = scode.2$pseudotime, "stage" = seurat.mouse$stage)
  stage = smooth.spline(as.numeric(data$stage),data$ex,all.knots = T, lambda = 0.001)
  pt = smooth.spline(data$pseudotime,data$ex,all.knots = T,lambda = 0.001)
  data.fit.stage = data.frame(x = seq(1,7,length.out = 1000), y = predict(stage,seq(1,7,length.out = 1000))$y)
  data.fit.pt = data.frame(x = seq(range(data$pseudotime)[1],range(data$pseudotime)[2],length.out = 1000), y = predict(pt,seq(range(data$pseudotime)[1],range(data$pseudotime)[2],length.out = 1000))$y)
  g1 = ggplot(data, aes(x = as.numeric(stage), y = ex, color = stage))+
    geom_line(data = data.fit.stage, aes(x = as.numeric(x), y = y), color = "black")+
    geom_point()+
    theme_bw()+
    theme(legend.position = "none")
  g2 = ggplot(data, aes(x = pseudotime, y = ex, color = stage))+
    geom_line(data = data.fit.pt, aes(x = x, y = y), color = "black")+
    geom_point()+
    theme_bw()+
    theme(legend.position = "none")
  g = grid.arrange(g1,g2,ncol =1, top = gn)
  print(g)
dev.off()
SSR[gn,]

#optimize D
library(MASS)
pnums = seq(1,20,0.5)
tfnum <- nrow(scode.1) #gene_number
cnum <- ncol(scode.1)# sample
maxite <- 100 # iter

maxB <- 2.0
minB <- -10.0

cal.RSS = function(x){
  RSSs = NULL
  pnum = pnums[x]
  X <- as.matrix(scode.1)
  W <- matrix(rep(0,tfnum*pnum), nrow=tfnum, ncol=pnum)
  Z <- matrix(rep(0,pnum*cnum), nrow=pnum, ncol=cnum)
  WZ <- matrix(nrow=tfnum, ncol=cnum)
  
  #read pseudo-time and normalize pseudo-time
  pseudotime <- scode.2$pseudotime
  pseudotime <- pseudotime/max(pseudotime)
  
  new_B <- rep(0, pnum)
  old_B <- rep(0, pnum)
  
  #initialization
  RSS <- Inf
  for(i in 1:pnum){
    new_B[i] <- runif(1, min=minB, max=maxB)
    old_B[i] <- new_B[i]
  }
  
  #function to sample Z
  sample_Z <- function(){
    for(i in 1:pnum){
      for(j in 1:cnum){
        Z[i,j] <<- exp(new_B[i]*pseudotime[j]) + runif(1, min=-0.001, max=0.001)
      }
    }
  }
  
  #optimize W and B iteratively
  for(ite in 1:maxite){
    cat(as.character(Sys.time()), "\t",RSS, "\t", ite,"\t",pnum,"\n")
    #sampling B
    target <- floor(runif(1, min=1, max=pnum+1))
    new_B[target] <- runif(1, min=minB, max=maxB)
    
    #for last calculation
    if(ite == maxite){
      for(i in 1:pnum){
        new_B[i] <- old_B[i]
      }
    }
    
    #sample Z from new B
    sample_Z()
    
    #regression
    for(i in 1:tfnum){
      X.lm <- lm(X[i,] ~ t(Z)-1)
      for(j in 1:pnum){
        W[i,j] <- X.lm$coefficients[j]
      }
      WZ[i,] <- W[i,] %*% Z
    }
    
    #RSS
    tmp_RSS <- sum((X-WZ)**2)
    if(tmp_RSS < RSS){
      RSS <- tmp_RSS
      old_B[target] <- new_B[target]
    }
    else{
      new_B[target] <- old_B[target]
    }
    RSSs = c(RSSs, RSS)
  }
  write.table(RSSs, file = sprintf("D_optimization/RSS_iter_%s",pnums[x]),row.names = F,col.names = F)
  write.table(RSS, file = sprintf("D_optimization/RSS_%s",pnums[x]),row.names = F,col.names = F)
}
dir.create("D_optimization")
library(doParallel)
registerDoParallel(20)
foreach (i=1:length(pnums)) %dopar% cal.RSS(i)

files = list.files("D_optimization/")
files = files[grep("iter", files, invert = T)]
RSSs = NULL
for (file in files) {
  RSSs = c(RSSs, read.table(sprintf("D_optimization/%s",file))[1,1,])
}
data = data.frame("D" = pnums, "RSS" = RSSs)
data = data[data$D<=10,]
pdf("optimizationD.pdf",height = 3.5)
ggplot(data, aes(x=D, y = RSS))+
  geom_line()+
  scale_x_continuous(breaks=seq(0,10,1))+
  theme_bw()
dev.off()

# Infering gene regulatory network with SCODE with D = 5
scode = function(i,pnum){
  RSSs = NULL
  dir = sprintf("SCODEmouse_out_all/out%s",i)
  dir.create(dir)
  X <- as.matrix(scode.1)
  W <- matrix(rep(0,tfnum*pnum), nrow=tfnum, ncol=pnum)
  Z <- matrix(rep(0,pnum*cnum), nrow=pnum, ncol=cnum)
  WZ <- matrix(nrow=tfnum, ncol=cnum)
  
  #read pseudo-time and normalize pseudo-time
  pseudotime <- scode.2$pseudotime
  pseudotime <- pseudotime/max(pseudotime)
  
  new_B <- rep(0, pnum)
  old_B <- rep(0, pnum)
  
  #initialization
  RSS <- Inf
  for(i in 1:pnum){
    new_B[i] <- runif(1, min=minB, max=maxB)
    old_B[i] <- new_B[i]
  }
  
  #function to sample Z
  sample_Z <- function(){
    for(i in 1:pnum){
      for(j in 1:cnum){
        Z[i,j] <<- exp(new_B[i]*pseudotime[j]) + runif(1, min=-0.001, max=0.001)
      }
    }
  }
  
  #optimize W and B iteratively
  for(ite in 1:maxite){
    cat(as.character(Sys.time()), "\t",RSS, "\t", ite,"\t",pnum,"\n")
    #sampling B
    target <- floor(runif(1, min=1, max=pnum+1))
    new_B[target] <- runif(1, min=minB, max=maxB)
    
    #for last calculation
    if(ite == maxite){
      for(i in 1:pnum){
        new_B[i] <- old_B[i]
      }
    }
    
    #sample Z from new B
    sample_Z()
    
    #regression
    for(i in 1:tfnum){
      X.lm <- lm(X[i,] ~ t(Z)-1)
      for(j in 1:pnum){
        W[i,j] <- X.lm$coefficients[j]
      }
      WZ[i,] <- W[i,] %*% Z
    }
    
    #RSS
    tmp_RSS <- sum((X-WZ)**2)
    if(tmp_RSS < RSS){
      RSS <- tmp_RSS
      old_B[target] <- new_B[target]
    }
    else{
      new_B[target] <- old_B[target]
    }
    RSSs = c(RSSs, RSS)
  }
  
  #output RSS
  write.table(RSS, paste(dir,"/RSS.txt",sep=""), row.names=F, col.names=F, sep="\t")
  
  #output W
  write.table(W, paste(dir,"/W.txt",sep=""), row.names=F, col.names=F, sep="\t")
  
  #infer A
  B <- matrix(rep(0,pnum*pnum), nrow=pnum, ncol=pnum)
  for(i in 1:pnum){
    B[i,i] <- new_B[i]
  }
  invW <- ginv(W)
  A <- W %*% B %*% invW
  
  #output A and B
  save(A,file = paste(dir,"/A",sep=""))
  #write.table(A, paste(dir,"/A.txt",sep=""), row.names=F, col.names=F, sep="\t")
  save(B, RSS, file = paste(dir,"/B-RSS",sep=""))
}
dir.create("SCODEmouse_out_all")
repnum = 20
registerDoParallel(repnum)
foreach (i=1:repnum) %dopar% scode(i,4,dir)

dir = "SCODEmouse_out_all"
load(paste(dir,"/out1/A",sep=""))
meanA <- A
for(i in 2:repnum){
  cat(i,"\n")
  load(paste(dir,"/out",i,"/A",sep=""))
  tmp <- A
  meanA <- meanA + tmp
  gc(gc())
}
meanA <- meanA / repnum

# Check reproducibility of A and RSS
dir = "SCODEmouse_out_all"
cors = NULL
for(a in 1:repnum){
  cat(a,"\n")
  load(paste(dir,"/out",a,"/A",sep=""))
  cors = c(cors,cor(matrix(A,ncol = 1),matrix(meanA,ncol = 1)))
  gc(gc(gc()))
}
pdf("reproducibilityA.pdf", height = 3.5)
top10 = order(cors, decreasing = T)[1:10]
ggplot(data.frame("cor" = cors),aes(x = 1:20,y = cor))+
  geom_point()+
  theme_bw()
dev.off()
gc(gc())
dir = "SCODEmouse_out_all"
RSSs = NULL
for(a in 1:repnum){
  cat(a,"\n")
  tmp = read.table(paste(dir,"/out",a,"/RSS.txt", sep = ""))[1,1]
  RSSs = c(RSSs, tmp)
}
plot(cors, RSSs)


# Calculate average A
dir = "SCODEmouse_out_all"
meanA = NULL
for(i in top10){
  cat(i,"\n")
  if(i == top10[1]){
    load(paste(dir,"/out1/A",sep=""))
    meanA <- A
  } else {
    load(paste(dir,"/out",i,"/A",sep=""))
    tmp <- A
    meanA <- meanA + tmp
  }
  gc(gc())
}
meanA <- meanA / repnum
res.scode = as.matrix(meanA)
colnames(res.scode) = rownames(scode.1)
rownames(res.scode) = rownames(scode.1)
save(res.scode, file = "res-scode")

# Check correlation of A and expression levels
mean = apply(seurat.mouse@assays$RNA@data, 1, mean)
cor.col = apply(res.scode, 2, FUN = function(x){
  return(cor(mean,abs(x)))
})
cor.row = apply(res.scode, 1, FUN = function(x){
  return(cor(mean,abs(x)))
})
data = data.frame("Regulator" = cor.col, "Target" = cor.row )
library(reshape)
data = melt(data)

pdf("A-expression.pdf",height = 3.5)
ggplot(data, aes(y = value, x = variable))+
  geom_violin(fill = "gray")+
  xlab(NULL)+
  ylab("Peason's correlation")+
  theme_bw()
dev.off()

# Define threathhold for each gene for GRNI
gn = "Sox8"
totals = abs(res.scode[,gn])
o <- order(totals, decreasing = TRUE)
stuff <- rle(totals[o])
run.rank <- cumsum(stuff$lengths) - (stuff$lengths - 1)/2
run.totals <- stuff$values
y <- (run.totals)
x <- (run.rank)
# fit.sp = smooth.spline(x, y, df =20)
# fit.lo = loess(y ~x)
fit <- nls(y ~ (a * x + b), start = c(a=0, b=0))
plot(y)
lines(predict(fit, x =x), col = "red")
# lines(predict(fit.lo), col = "blue")
# lines(predict(fit.sp), col = "green")
top = which((y - predict(fit,x = x))<0)[1]-1
bottom = which((rev(y - predict(fit,x = x)))>0)[1]-1
# top = which((y - predict(fit)$y)<0)[1]-1
# bottom = which((rev(y - predict(fit)$y))>0)[1]-1
bottom = length(y)-bottom+1
abline(v = c(top, bottom), col = "blue")
mean = apply(seurat.mouse@assays$RNA@data, 1,mean)
library(gridExtra)
library(ggplot2)
pdf("A-mean-sox8.pdf",height = 5)
g1 = ggplot(data.frame("A" = abs(res.scode[,gn]), "Expression" = mean), aes(x=Expression, y = A))+
  geom_point(color = "gray")+
  xlab("Average expression of targets")+
  theme_bw()
g2 = ggplot(data.frame("A" = abs(res.scode[gn,]), "Expression" = mean), aes(x=Expression, y = A))+
  geom_point(color = "gray")+
  xlab("Average expression of regulators")+
  theme_bw()
grid.arrange(g2,g1, nrow = 2)
dev.off()
pdf("A-threshold-sox8.pdf", width = 9)
ggplot(data.frame("A" = y, "Index" = 1:length(y)), aes(x=Index, y = A))+
  geom_point(color = "gray")+
  geom_line(data = data.frame("A"=predict(fit, x =x),"Index" = 1:length(y)), color = "black")+
  # geom_vline(xintercept = c(top, bottom), linetype="dashed")+
  geom_hline(yintercept = fitted(fit)[c(top, bottom)], linetype="dashed")+
  ggtitle(gn)+
  theme_bw()
dev.off()

cal.knee = function(totals){
  totals = abs(totals)
  o <- order(totals, decreasing = TRUE)
  stuff <- rle(totals[o])
  run.rank <- cumsum(stuff$lengths) - (stuff$lengths - 1)/2
  run.totals <- stuff$values
  y <- (run.totals)
  x <- (run.rank)
  fit <- nls(y ~ (a * x + b), start = c(a=0, b=0))
  top = which((y - predict(fit,x = x))<0)[1]-1
  knee <- predict(fit)[top]
  return(knee)
}

knees = NULL
for(gn in colnames(res.scode)){
  cat(gn, "\n")
  knees = c(knees, cal.knee(res.scode[,gn]))
}
names(knees) = rownames(res.scode)  
save(knees, file = "knees")

regulation.mat = matrix(0,ncol = ncol(res.scode), nrow = nrow(res.scode))
colnames(regulation.mat) = colnames(res.scode)
rownames(regulation.mat) = rownames(res.scode)
for(gn in colnames(res.scode)){
  pass = abs(res.scode[,gn])>knees[gn] 
  target = colnames(res.scode)[pass]
  order = order(abs(res.scode[,gn])[target],decreasing = T)
  target = target[order]
  regulation.mat[target,gn] = 1:length(target)
}
regulation.mat.knee = regulation.mat

# TF2DNA
answer = NULL
list = list.dirs("TF2DNAdb")
list = list[grep("Mus",list)]
for (dir in list) {
  files = list.files(dir)
  for (file in files) {
    cat(file, "\n")
    tmp = read.table(sprintf("%s/%s",dir,file),header = T,sep = "\t")[,1:2]
    answer = rbind(answer,tmp)
  }
}
library(dplyr)
answer = answer %>% distinct(tf_name,target_name)
TF.list = colnames(res.scode)[is.element(colnames(res.scode),unique(answer$tf_name))]
save(answer, TF.list, file = "TF2DNA")
correct.rate = NULL
all.gene = colnames(regulation.mat)
for (gn in TF.list) {
  tmp3 = colnames(regulation.mat)[regulation.mat[,gn]>0]
  tmp3 = tmp3[is.element(tmp3, all.gene)]
  tmp3 = unique(tmp3)
  tmp2 = unique(answer[answer$tf_name==gn,]$target_name)
  tmp2 = tmp2[is.element(tmp2, all.gene)]
  mat = matrix(c( (length(all.gene) - length(tmp2)) ,(length(tmp3)-sum(is.element(tmp3,tmp2))),length(tmp2),sum(is.element(tmp3,tmp2))),byrow = T,ncol = 2)
  tmp = data.frame("gene.name" = gn,"Num.of.predicted.target" = length(tmp3),"correct.rate" = sum(is.element(tmp3,tmp2))/length(tmp3), "target.rate" = (length(tmp2)/length(all.gene)),"p.value" = fisher.test(mat)$p.value, "Num.of.target" = length(tmp2))
  correct.rate = rbind(correct.rate,tmp)
}
correct.rate = cbind(correct.rate,data.frame("q.value" = p.adjust(correct.rate$p.value,method = "fdr")))
mean(correct.rate$Num.of.predicted.target)
mean(correct.rate$correct.rate)
mean(correct.rate$target.rate)
sum(correct.rate$q.value<0.05)


data = correct.rate
library(ggplot2)
library(gridExtra)
g1 = ggplot(data, aes(x = target.rate))+
  geom_histogram(fill = "black")+
  scale_x_continuous(breaks=seq(0,1,by=0.1),limits=c(0,1))+
  theme_bw()
g2 = ggplot(data, aes(x = correct.rate))+
  geom_histogram(fill = "black")+
  scale_x_continuous(breaks=seq(0,1,by=0.1),limits=c(0,1))+
  coord_flip()+  
  theme_bw()
g3 =ggplot(data, aes(x = target.rate, y = correct.rate))+
  scale_x_continuous(breaks=seq(0,1,by=0.1),limits=c(0,1))+
  scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,1))+
  geom_point(size = 2, color = "black")+
  geom_abline(slope = 1,intercept = 0)+
  theme_bw()
g4 =ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())

pdf("Fig1d.pdf")
grid.arrange(g1, g4, g3, g2, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()
data = data.frame("TF" = correct.rate$gene.name,"Total" = correct.rate$Num.of.predicted.target,"Positive" = correct.rate$Num.of.predicted.target*correct.rate$correct.rate, "Negative" = correct.rate$Num.of.predicted.target-(correct.rate$Num.of.predicted.target*correct.rate$correct.rate))
data = data[base::order(data$Total, decreasing = T),]
library(reshape2)
data = melt(data, id.vars = c("TF","Total"))
data$TF = factor(data$TF, levels = as.character(unique(data$TF)))
data$variable = factor(data$variable, levels = c("Negative","Positive"))
pdf("Top30TF.pdf")
ggplot(data[is.element(data$TF,data$TF[1:30]),], aes(x = TF,y=value, fill = variable))+
  geom_bar(stat = "identity", color = "black")+
  scale_fill_manual(values=c("black","white")) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "none")
dev.off()

ex = data.frame("Tcf4" = seurat.mouse@assays$RNA@data["Tcf4",],"Tcf21" = seurat.mouse@assays$RNA@data["Tcf21",],"Pseudotime" = rep(scode.2$pseudotime))
g1 = ggplot(ex, aes(x = Pseudotime, y = Tcf4))+
  geom_smooth(se = F,method = "loess",color = "gray")+
  geom_point()+
  ggtitle("Tcf4")+
  theme_bw()
g2 = ggplot(ex, aes(x = Pseudotime, y = Tcf21))+
  geom_smooth(se = F,method = "loess",color = "gray")+
  geom_point()+
  ggtitle("Tcf21")+
  theme_bw()
pdf("Tcf4-21-explot.pdf")
grid.arrange(g1,g2, ncol = 1)
dev.off()

regulation.mat.rank = apply(res.scode,2,FUN = function(x){
  x2 = order(abs(x), decreasing = T)
  x2[x2] = 1:length(x2)
  return(x2)
})
rownames(regulation.mat.rank) = rownames(res.scode)

PR.diff = NULL
pdf("TF_th_check.pdf",height = 3.5)
for(gn in c("E2f3","Arid5b")){
  cat(gn,"\n")
  out = NULL
  for (th in seq(100,nrow(regulation.mat),100)) {
    tmp3 = colnames(regulation.mat)[regulation.mat.rank[,gn]<=th]
    tmp2 = unique(answer[answer$tf_name==gn,]$target_name)
    tmp2 = tmp2[is.element(tmp2, colnames(regulation.mat))]
    out = rbind(out, data.frame("th" = th,"PR" = sum(is.element(tmp3,tmp2))/length(tmp3)))
  }
  data = data.frame("th" = out$th, "PR" = out$PR)
  ss = out$th[out$PR==max(out$PR)]
  th = max(regulation.mat[,gn])
  PR.diff = rbind(PR.diff, data.frame("TF" = gn,"ideal" = ss, "used.th" = th, "PR" = max(out$PR)-correct.rate$correct.rate[correct.rate$gene.name==gn]))
  g = ggplot(data, aes(x=th, y = PR))+
    geom_point()+
    ylim(c(0.3,0.95))+
    geom_vline(xintercept = th)+
    geom_vline(xintercept = ss, linetype="dashed")+
    ggtitle(gn)+
    theme_bw()
  # abs(th-ss)
  print(g)
}
dev.off()

# ROC
TF.target =  matrix(0,ncol = ncol(res.scode), nrow = nrow(res.scode))
colnames(TF.target) = colnames(res.scode)
rownames(TF.target) = rownames(res.scode)
for(gn in TF.list){
    TF.target[unique(answer$target_name[answer$tf_name==gn])[is.element(unique(answer$target_name[answer$tf_name==gn]),rownames(res.scode))],gn] = 1
}
library(ROCR)

aucs = NULL
for(gn in TF.list){
  cat(gn,"\n")
  scode = abs(res.scode[,gn])
  answer.TF = TF.target[,gn]
  pred <- prediction(scode, answer.TF)
  perf <- performance(pred, "tpr", "fpr")
  plot(perf)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  aucs = c(aucs, auc)
}

gn = "E2f3"

scode = abs(res.scode[,gn])
answer.TF = TF.target[,gn]
pred <- prediction(scode, answer.TF)
perf <- performance(pred, "tpr", "fpr")
data = data.frame("FP" = perf@x.values[[1]], "TP" = perf@y.values[[1]])
g1 = ggplot(data, aes(x = FP, y = TP))+
  geom_point()+
  geom_abline(slope = 1,intercept = 0)+
  ggtitle(gn)+
  theme_bw()

gn = "Arid5b"
scode = abs(res.scode[,gn])
answer.TF = TF.target[,gn]
pred <- prediction(scode, answer.TF)
perf <- performance(pred, "tpr", "fpr")
data = data.frame("FP" = perf@x.values[[1]], "TP" = perf@y.values[[1]])
g2 = ggplot(data, aes(x = FP, y = TP))+
  geom_point()+
  geom_abline(slope = 1,intercept = 0)+
  ggtitle(gn)+
  theme_bw()

grid.arrange(g1,g2, ncol = 1)
names(aucs) = TF.list
aucs.scode = aucs
pdf("AUC.pdf",width = 3.5)
ggplot(data.frame("AUC" = aucs), aes(y = AUC, x = 1))+
  geom_beeswarm()+
  theme_bw()
dev.off()
which(aucs == max(aucs))
which(aucs == min(aucs))

ligand.recepter.num = function(a, mat){
  pair = as.character(a[1])
  l = as.character(a[2])
  r = as.character(a[3])
  pair = unlist(strsplit(pair,"-"))
  ligand = unlist(strsplit(l,", "))
  recepter = unlist(strsplit(r,", "))
  ligand = ligand[is.element(ligand, rownames(mat))]
  recepter = recepter[is.element(recepter, rownames(mat))]
  if(length(ligand)>1){
    ligand.target = rownames(mat)[rowSums(mat[,ligand]>0)>0]
  } else {
    ligand.target = rownames(mat)[mat[,ligand]>0]
  }
  if(length(recepter)>1){
    recepter.target = rownames(mat)[rowSums(mat[,recepter]>0)>0]
  } else {
    recepter.target = rownames(mat)[mat[,recepter]>0]
  }
  ligand.target.positive = toy = list("Ligand" = ligand.target, "Recepter" = recepter.target)
  names(toy) = pair
  toy = Venn(toy)
  p =  ggvenn(toy,fill = c("white","white"))+theme_void()
  common = sum(is.element(ligand.target, recepter.target))
  toy = list("All" = all.gene, "Ligand" = ligand.target, "Recepter" = recepter.target)
  names(toy)[2:3] = pair
  toy = Venn(toy)
  en = enrichment_test(toy,2,3)
  out = data.frame("Ligand.name" = pair[1],"Recepter.name" = pair[2],"Num.ligand.target" = length(ligand.target), "Num.recepter.target" = length(recepter.target), "Num.overlap" = common, "p.value" =  en$Significance)
  return(list("Plot" = p, "Table" = out,"ligand" = pair[1],"receptor" = pair[2],"ligand.target" = ligand.target,"receptor.target" = recepter.target))
}
mat = regulation.mat
ligand.target.table = NULL
list = read.csv("TableS2.csv")
pdf("ligand-recepter.pdf")
for (i in 1:nrow(list)) {
  a = ligand.recepter.num(list[i,],mat)
  ligand.target.table = rbind(ligand.target.table,a$Table)
}
dev.off()

write.csv(correct.rate, file = "TableS1.csv")
write.csv(ligand.target.table, file = "ligand-target.csv")

data = data.frame("Unique recepter targets" = ligand.target.table$Num.recepter.target-ligand.target.table$Num.overlap,"Common targets" = ligand.target.table$Num.overlap,"Unique ligand targets" = ligand.target.table$Num.ligand.target-ligand.target.table$Num.overlap, "pair" = sprintf("%s-%s",ligand.target.table$Ligand.name,ligand.target.table$Recepter.name))
library(reshape2)
data = melt(data)
pdf("ligand-recepter.pdf")
ggplot(data, aes(x = value, y = pair, fill = variable))+
  geom_bar(stat = "identity", color = "black")+
  scale_fill_manual(values = c("gray40","white","gray"))+
  theme_bw()
dev.off()

# Top ligand and receotor
library(UpSetR)
pdf("ligand-receptor-top.pdf", height = 3.5)
for(i in 1:nrow(list)){
  gnl = list[i,]
  ligand  = unlist(strsplit(gnl$Ligand,", "))
  if(length(ligand)>1){
    tmp = colSums(regulation.mat[,ligand]>0)
    ligand = names(which(tmp==max(tmp))[1])
  }
  receptor  = unlist(strsplit(gnl$Receptor,", "))
  if(length(receptor)>1){
    tmp = colSums(regulation.mat[,receptor]>0)
    receptor = names(which(tmp==max(tmp)))
  }
  ligand.positive = names(res.scode[regulation.mat[,ligand]>0, ligand])[res.scode[regulation.mat[,ligand]>0, ligand]>0]
  ligand.negative = names(res.scode[regulation.mat[,ligand]>0, ligand])[res.scode[regulation.mat[,ligand]>0, ligand]<0]
  receptor.positive = names(res.scode[regulation.mat[,receptor]>0, receptor])[res.scode[regulation.mat[,receptor]>0, receptor]>0]
  receptor.negative = names(res.scode[regulation.mat[,receptor]>0, receptor])[res.scode[regulation.mat[,receptor]>0, receptor]<0]
  toy = list(ligand.positive, ligand.negative, receptor.positive, receptor.negative)
  names(toy) = c(sprintf("%s_positive", ligand),sprintf("%s_negative", ligand),sprintf("%s_positive", receptor),sprintf("%s_negative", receptor))
  p = upset(fromList(toy), order.by = "freq", mb.ratio = c(0.2,0.8),point.size = 10,set_size.show = T,keep.order = T, sets = c(sprintf("%s_positive", ligand),sprintf("%s_positive", receptor),sprintf("%s_negative", ligand),sprintf("%s_negative", receptor)), intersections = list(list(sprintf("%s_positive", ligand),sprintf("%s_positive", receptor)),list(sprintf("%s_negative", ligand),sprintf("%s_negative", receptor)),list(sprintf("%s_positive", ligand),sprintf("%s_negative", receptor)),list(sprintf("%s_negative", ligand),sprintf("%s_positive", receptor))))
  print(p)
  # toy1 = list(ligand.positive, ligand.negative, receptor.positive)
  # names(toy1) = c(sprintf("%s_positive", ligand),sprintf("%s_negative", ligand),sprintf("%s_positive", receptor))
  # toy2 = list(ligand.positive, ligand.negative, receptor.negative)
  # names(toy2) = c(sprintf("%s_positive", ligand),sprintf("%s_negative", ligand),sprintf("%s_negative", receptor))
  # toy1= Venn(toy1)
  # toy2= Venn(toy2)
  # p1 = ggvenn(toy1)+theme_void()+NoLegend()
  # p2 = ggvenn(toy2)+theme_void()+NoLegend()
  # grid.arrange(p1,p2,ncol = 2, top = gnl$Signaling.pathway)
}
dev.off()

pdf("wnt-Ctnnb1-top.pdf", height = 3.5)
ligand = "Wnt7a"
receptor = "ctnnb1"
ligand.positive = names(res.scode[regulation.mat[,ligand]>0, ligand])[res.scode[regulation.mat[,ligand]>0, ligand]>0]
ligand.negative = names(res.scode[regulation.mat[,ligand]>0, ligand])[res.scode[regulation.mat[,ligand]>0, ligand]<0]
receptor.positive = names(res.scode[regulation.mat[,receptor]>0, receptor])[res.scode[regulation.mat[,receptor]>0, receptor]>0]
receptor.negative = names(res.scode[regulation.mat[,receptor]>0, receptor])[res.scode[regulation.mat[,receptor]>0, receptor]<0]
toy = list(ligand.positive, ligand.negative, receptor.positive, receptor.negative)
names(toy) = c(sprintf("%s_positive", ligand),sprintf("%s_negative", ligand),sprintf("%s_positive", receptor),sprintf("%s_negative", receptor))
p = upset(fromList(toy), order.by = "freq", mb.ratio = c(0.2,0.8),point.size = 10,set_size.show = T,keep.order = T, sets = c(sprintf("%s_positive", ligand),sprintf("%s_positive", receptor),sprintf("%s_negative", ligand),sprintf("%s_negative", receptor)), intersections = list(list(sprintf("%s_positive", ligand),sprintf("%s_positive", receptor)),list(sprintf("%s_negative", ligand),sprintf("%s_negative", receptor)),list(sprintf("%s_positive", ligand),sprintf("%s_negative", receptor)),list(sprintf("%s_negative", ligand),sprintf("%s_positive", receptor))))
print(p)
dev.off()


pdf("wnt-Fzd-top.pdf", height = 3.5)
ligand = "Wnt7a"
receptor = "Fzd6"
ligand.positive = names(res.scode[regulation.mat[,ligand]>0, ligand])[res.scode[regulation.mat[,ligand]>0, ligand]>0]
ligand.negative = names(res.scode[regulation.mat[,ligand]>0, ligand])[res.scode[regulation.mat[,ligand]>0, ligand]<0]
receptor.positive = names(res.scode[regulation.mat[,receptor]>0, receptor])[res.scode[regulation.mat[,receptor]>0, receptor]>0]
receptor.negative = names(res.scode[regulation.mat[,receptor]>0, receptor])[res.scode[regulation.mat[,receptor]>0, receptor]<0]
toy = list(ligand.positive, ligand.negative, receptor.positive, receptor.negative)
names(toy) = c(sprintf("%s_positive", ligand),sprintf("%s_negative", ligand),sprintf("%s_positive", receptor),sprintf("%s_negative", receptor))
p = upset(fromList(toy), order.by = "freq", mb.ratio = c(0.2,0.8),point.size = 10,set_size.show = T,keep.order = T, sets = c(sprintf("%s_positive", ligand),sprintf("%s_positive", receptor),sprintf("%s_negative", ligand),sprintf("%s_negative", receptor)), intersections = list(list(sprintf("%s_positive", ligand),sprintf("%s_positive", receptor)),list(sprintf("%s_negative", ligand),sprintf("%s_negative", receptor)),list(sprintf("%s_positive", ligand),sprintf("%s_negative", receptor)),list(sprintf("%s_negative", ligand),sprintf("%s_positive", receptor))))
print(p)
ligand = "Wnt8a"
receptor = "Fzd5"
ligand.positive = names(res.scode[regulation.mat[,ligand]>0, ligand])[res.scode[regulation.mat[,ligand]>0, ligand]>0]
ligand.negative = names(res.scode[regulation.mat[,ligand]>0, ligand])[res.scode[regulation.mat[,ligand]>0, ligand]<0]
receptor.positive = names(res.scode[regulation.mat[,receptor]>0, receptor])[res.scode[regulation.mat[,receptor]>0, receptor]>0]
receptor.negative = names(res.scode[regulation.mat[,receptor]>0, receptor])[res.scode[regulation.mat[,receptor]>0, receptor]<0]
toy = list(ligand.positive, ligand.negative, receptor.positive, receptor.negative)
names(toy) = c(sprintf("%s_positive", ligand),sprintf("%s_negative", ligand),sprintf("%s_positive", receptor),sprintf("%s_negative", receptor))
p = upset(fromList(toy), order.by = "freq", mb.ratio = c(0.2,0.8),point.size = 10,set_size.show = T,keep.order = T, sets = c(sprintf("%s_positive", ligand),sprintf("%s_positive", receptor),sprintf("%s_negative", ligand),sprintf("%s_negative", receptor)), intersections = list(list(sprintf("%s_positive", ligand),sprintf("%s_positive", receptor)),list(sprintf("%s_negative", ligand),sprintf("%s_negative", receptor)),list(sprintf("%s_positive", ligand),sprintf("%s_negative", receptor)),list(sprintf("%s_negative", ligand),sprintf("%s_positive", receptor))))
print(p)
dev.off()

pdf("Wnt-explot.pdf", width = 14, height = 21)
VlnPlot(seurat.mouse, features = c("Wnt1", "Wnt2" , "Wnt2b" , "Wnt3"  , "Wnt3a", "Wnt4" , "Wnt5a" , "Wnt5b", "Wnt6" ,   "Wnt7a",   "Wnt7b" , "Wnt8a","Wnt8b", "Wnt9a" , "Wnt9b" ,  "Wnt10a","Wnt10b",      "Wnt11"  ,  "Wnt16" ),group.by = "stage", ncol = 2)
dev.off()

pdf("Fzd-explot.pdf", width = 14, height = 10.5)
VlnPlot(seurat.mouse, features = c( "Fzd1", "Fzd2","Fzd3" ,"Fzd4","Fzd5", "Fzd6",  "Fzd7"  ,  "Fzd8","Fzd9","Fzd10"      ),group.by = "stage", ncol = 2)
dev.off()

pdf("headihog-top.pdf", height = 3.5)
ligand = "Ssh2"
receptor = "Gli1"
ligand.positive = names(res.scode[regulation.mat[,ligand]>0, ligand])[res.scode[regulation.mat[,ligand]>0, ligand]>0]
ligand.negative = names(res.scode[regulation.mat[,ligand]>0, ligand])[res.scode[regulation.mat[,ligand]>0, ligand]<0]
receptor.positive = names(res.scode[regulation.mat[,receptor]>0, receptor])[res.scode[regulation.mat[,receptor]>0, receptor]>0]
receptor.negative = names(res.scode[regulation.mat[,receptor]>0, receptor])[res.scode[regulation.mat[,receptor]>0, receptor]<0]
toy = list(ligand.positive, ligand.negative, receptor.positive, receptor.negative)
names(toy) = c(sprintf("%s_positive", ligand),sprintf("%s_negative", ligand),sprintf("%s_positive", receptor),sprintf("%s_negative", receptor))
p = upset(fromList(toy), order.by = "freq", mb.ratio = c(0.2,0.8),point.size = 10,set_size.show = T,keep.order = T, sets = c(sprintf("%s_positive", ligand),sprintf("%s_positive", receptor),sprintf("%s_negative", ligand),sprintf("%s_negative", receptor)), intersections = list(list(sprintf("%s_positive", ligand),sprintf("%s_positive", receptor)),list(sprintf("%s_negative", ligand),sprintf("%s_negative", receptor)),list(sprintf("%s_positive", ligand),sprintf("%s_negative", receptor)),list(sprintf("%s_negative", ligand),sprintf("%s_positive", receptor))))
print(p)
dev.off()

# dynGENIE3
source("dynGENIE3/dynGENIE3_R_C_wrapper/dynGENIE3.R")
# TS.data: a list of matrices containing time series expression data. Each
# matrix (genes x time points) corresponds to a time series experiment. Each
# row of a matrix is a gene, each column is a time point.
# time.points: a list of vectors containing the time points. The k-th vector
# must correspond to the k-th time series of TS.data
# TS1 <- read.expr.matrix("dynGENIE3/dynGENIE3_R_C_wrapper/example_data/time_series_1.txt", form="rows.are.samples")
# tmp = TS1[2:nrow(TS1),]
# TS.data <- list(tmp)
# time.points <- list(TS1[1,])
# TS1 = read.expr.matrix("dynGENIE3/dynGENIE3_R_C_wrapper/example_data/time_series_1.txt", form="rows.are.samples")
# TS2 = read.expr.matrix("dynGENIE3/dynGENIE3_R_C_wrapper/example_data/time_series_2.txt", form="rows.are.samples")
# TS3 = read.expr.matrix("dynGENIE3/dynGENIE3_R_C_wrapper/example_data/time_series_3.txt", form="rows.are.samples")
# time.points = list(TS1[1,],TS2[1,],TS3[1,])
# TS.data = list(TS1[2:nrow(TS1),],TS2[2:nrow(TS2),],TS3[2:nrow(TS3),])
tmp2 = as.integer(seurat.mouse@meta.data$stage)
names(tmp2) = colnames(tmp)
tmp = as.array(seurat.mouse@assays$RNA@data)
rownames(tmp) = rownames(seurat.mouse@assays$RNA@data)
colnames(tmp) = sprintf("sample_%s",colnames(seurat.mouse@assays$RNA@data))
TS.data = NULL
time.points = NULL
for(a in 1:11){
  sample.list = NULL
  time.list = NULL
  for(i in 1:7){
    if((a!=11) | (i !=1)){
      sample.list = c(sample.list,which(tmp2==i)[a])
      time.list = c(time.list,i)
    }
  }
  ex = tmp[,sample.list]
  names(time.list) = colnames(ex)
  TS.data = c(TS.data,list(ex))
  time.points = c(time.points,list(time.list-1))
}
pass = NULL
for(a in 1:length(TS.data)){
  if(a==1){
    pass = rowSums(TS.data[[a]]>0)==0
  } else {
    pass = pass + (rowSums(TS.data[[a]]>0)==0)
  }
}
pass = names(which(pass==0))
for(a in 1:length(TS.data)){
  TS.data[[a]] = TS.data[[a]][pass,]
}
  res.dynGENIE3 <- dynGENIE3(TS.data,time.points,ncores = 32)
res.dynGENIE3 = t(res.dynGENIE3$weight.matrix)
save(res.dynGENIE3,file = "res-dynGENIE3-real-time")
plot(sort(res.dynGENIE3[,4]))
load("TF2DNA")

#BETS
tmp2 = as.integer(seurat.mouse@meta.data$stage)
names(tmp2) = colnames(tmp)
tmp = as.array(seurat.mouse@assays$RNA@scale.data)
rownames(tmp) = rownames(seurat.mouse@assays$RNA@scale.data)
colnames(tmp) = sprintf("sample_%s",colnames(seurat.mouse@assays$RNA@scale.data))
TS.data = NULL
time.points = NULL
for(a in 1:11){
  sample.list = NULL
  time.list = NULL
  for(i in 1:7){
    if((a!=11) | (i !=1)){
      sample.list = c(sample.list,which(tmp2==i)[a])
      time.list = c(time.list,i)
    }
  }
  ex = tmp[,sample.list]
  names(time.list) = colnames(ex)
  TS.data = c(TS.data,list(ex))
  time.points = c(time.points,list(time.list-1))
}
pass = NULL
for(a in 1:length(TS.data)){
  if(a==1){
    pass = rowSums(TS.data[[a]]>0)==0
  } else {
    pass = pass + (rowSums(TS.data[[a]]>0)==0)
  }
}
pass = names(which(pass==0))
for(a in 1:length(TS.data)){
  TS.data[[a]] = TS.data[[a]][pass,]
}
for(i in 1:10){
  write.table(TS.data[[i]],quote = F,file = sprintf("BETS/data/DREAM/insilico_size100_1/0mean/mouse-%s.txt",i),sep = "\t")
}

for(i in 1:10){
  write.table(t(apply(TS.data[[i]],1,sample)),quote = F,file = sprintf("BETS/data/DREAM/insilico_size100_1/0mean/mouse-%s-rand.txt",i),sep = "\t")
}

write.table(rownames(TS.data[[i]],1,sample),row.names = F,col.names = F,quote = F,file = "BETS/data/DREAM/insilico_size100_1/0mean/genes.txt",sep = "\t")

# ROC
TF.target =  matrix(0,ncol = ncol(res.dynGENIE3), nrow = nrow(res.dynGENIE3))
colnames(TF.target) = colnames(res.dynGENIE3)
rownames(TF.target) = rownames(res.dynGENIE3)
TF.list = TF.list[is.element(TF.list,colnames(TF.target))]
for(gn in TF.list){
  TF.target[unique(answer$target_name[answer$tf_name==gn])[is.element(unique(answer$target_name[answer$tf_name==gn]),rownames(res.dynGENIE3))],gn] = 1
}
library(ROCR)

aucs = NULL
for(gn in TF.list){
  cat(gn,"\n")
  scode = abs(res.dynGENIE3[,gn])
  answer.TF = TF.target[,gn]
  pred <- prediction(scode, answer.TF)
  perf <- performance(pred, "tpr", "fpr")
  plot(perf)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  aucs = c(aucs, auc)
}

gn = "E2f3"

scode = abs(res.dynGENIE3[,gn])
answer.TF = TF.target[,gn]
pred <- prediction(scode, answer.TF)
perf <- performance(pred, "tpr", "fpr")
data = data.frame("FP" = perf@x.values[[1]], "TP" = perf@y.values[[1]])
g1 = ggplot(data, aes(x = FP, y = TP))+
  geom_point()+
  geom_abline(slope = 1,intercept = 0)+
  ggtitle(gn)+
  theme_bw()

gn = "Arid5b"
scode = abs(res.dynGENIE3[,gn])
answer.TF = TF.target[,gn]
pred <- prediction(scode, answer.TF)
perf <- performance(pred, "tpr", "fpr")
data = data.frame("FP" = perf@x.values[[1]], "TP" = perf@y.values[[1]])
g2 = ggplot(data, aes(x = FP, y = TP))+
  geom_point()+
  geom_abline(slope = 1,intercept = 0)+
  ggtitle(gn)+
  theme_bw()

grid.arrange(g1,g2, ncol = 1)

names(aucs) = TF.list
aucs.dyGENIE3 = aucs
save(aucs.dyGENIE3, file = "aucs-dyGENIE3")
pdf("AUC_dynGENIE3.pdf",width = 3.5)
ggplot(data.frame("AUC" = aucs), aes(y = AUC, x = 1))+
  geom_beeswarm()+
  theme_bw()
dev.off()
which(aucs == max(aucs))
which(aucs == min(aucs))

data = rbind(data.frame("Method" = "scode","AUC" = aucs.scode), data.frame("Method"= "dyGENIE3","AUC" = aucs.dyGENIE3))
data$Method = factor(data$Method, levels = c("scode","dyGENIE3"))
pdf("fig3a.pdf")
ggplot(data, aes(x = Method,y = AUC))+
  geom_beeswarm()+theme_bw()
dev.off()

############ SCODE with stage
scode.1 =  as.matrix(seurat.mouse@assays$RNA@data)
scode.2 = data.frame("id" = 1:ncol(scode.1), "pseudotime" = as.numeric(seurat.mouse@meta.data$stage))

###SCODE
# Infering gene regulatory network with SCODE with D = 5
library(MASS)
tfnum <- nrow(scode.1) #gene_number
cnum <- ncol(scode.1)# sample
maxite <- 100 # iter

dir = sprintf("SCODE_stage")
dir.create(dir)
repnum = 20
library(doParallel)
registerDoParallel(repnum)
foreach (i=1:repnum) %dopar% scode(i,4,dir)

load(paste(dir,"/out1/A",sep=""))
meanA <- A
for(i in 2:repnum){
  cat(i,"\n")
  load(paste(dir,"/out",i,"/A",sep=""))
  tmp <- A
  meanA <- meanA + tmp
  gc(gc())
}
meanA <- meanA / repnum

# Check reproducibility of A and RSS
cors = NULL
for(a in 1:repnum){
  cat(a,"\n")
  load(paste(dir,"/out",a,"/A",sep=""))
  cors = c(cors,cor(matrix(A,ncol = 1),matrix(meanA,ncol = 1)))
  gc(gc(gc()))
}
top10 = order(cors, decreasing = T)[1:10]
ggplot(data.frame("cor" = cors),aes(x = 1:20,y = cor))+
  geom_point()+
  theme_bw()
gc(gc())


# Calculate average A
meanA = NULL
for(i in top10){
  cat(i,"\n")
  if(i == top10[1]){
    load(paste(dir,"/out1/A",sep=""))
    meanA <- A
  } else {
    load(paste(dir,"/out",i,"/A",sep=""))
    tmp <- A
    meanA <- meanA + tmp
  }
  gc(gc())
}
meanA <- meanA / 10
res.scode = as.matrix(meanA)
# gnl = NULL
# for(gn in rownames(scode.1)){
#   tmp = unique(transcript2gene.mouse$GeneName[transcript2gene.mouse$GeneID==gn])
#   if(length(tmp)>0){
#     names(tmp) = gn
#     gnl = c(gnl,tmp)
#   }
# }
colnames(res.scode) = rownames(scode.1)
rownames(res.scode) = rownames(scode.1)
# res.scode = res.scode[names(gnl),names(gnl)]
save(res.scode, file = sprintf("res-scode-stage"))

# ROC
load("TF2DNA")
TF.target =  matrix(0,ncol = ncol(res.scode), nrow = nrow(res.scode))
colnames(TF.target) = colnames(res.scode)
rownames(TF.target) = rownames(res.scode)
for(gn in TF.list[is.element(TF.list, rownames(res.scode))]){
  TF.target[unique(answer$target_name[answer$tf_name==gn])[is.element(unique(answer$target_name[answer$tf_name==gn]),rownames(res.scode))],gn] = 1
}
library(ROCR)

aucs = NULL
for(gn in TF.list[is.element(TF.list, rownames(res.scode))]){
  cat(gn,"\n")
  scode = abs(res.scode[,gn])
  answer.TF = TF.target[,gn]
  pred <- prediction(scode, answer.TF)
  perf <- performance(pred, "tpr", "fpr")
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  aucs = c(aucs, auc)
}
save(aucs,file = sprintf("AUC-SCODE-stage"))
library(ggplot2)
aucs.stage = aucs
data = rbind(data.frame("Method" = "scode","AUC" = aucs.scode), data.frame("Method"= "scode.stage","AUC" = aucs.stage), data.frame("Method"= "dyGENIE3","AUC" = aucs.dyGENIE3))
data$Method = factor(data$Method, levels = c("scode","scode.stage","dyGENIE3"))
median(aucs.scode)
median(aucs.stage)
t.test(aucs.scode,aucs.stage)
pdf("fig3a.pdf")
ggplot(data, aes(x = Method,y = AUC))+
  geom_beeswarm()+theme_bw()
dev.off()
data = data.frame("pseudotime"=aucs.scode,"stage" = aucs.stage)
ggplot(data, aes(x = stage,y = pseudotime))+
  geom_abline(slope = 1,intercept = 0)+
  geom_point()+
  theme_bw()


# Cardoso et al, 2019
library(Seurat)
library(doParallel)
library(ggplot2)
cpm.cardoso = read.table("Mouse.CPM.txt")
id = colnames(cpm.cardoso)
tissue = unlist(lapply(id, FUN = function(x){
  x = strsplit(x,"\\.")
  return(x[[1]][1])
}))
stage = unlist(lapply(id, FUN = function(x){
  x = strsplit(x,"\\.")
  return(x[[1]][2])
}))
at.cardoso = data.frame("stage" = stage, "tissue" = tissue)
rownames(at.cardoso) = id
library(Seurat)
seurat.carsodo = CreateSeuratObject(cpm.cardoso)
seurat.carsodo@meta.data$stage = at.cardoso$stage
seurat.carsodo@meta.data$tissue = at.cardoso$tissue

# Pseudotime analysis
rm(meanA,A, answer, TF.list,TF.target)
gc()
gc()
gc()
for(tissue in unique(at.cardoso$tissue)){
  sample = rownames(at.cardoso)[at.cardoso$tissue==tissue]
  seurat.carsodo.Heart = subset(seurat.carsodo, cell = sample)
  seurat.carsodo.Heart = NormalizeData(seurat.carsodo.Heart)
  seurat.carsodo.Heart = ScaleData(seurat.carsodo.Heart,do.scale = F)
  seurat.carsodo.Heart = FindVariableFeatures(seurat.carsodo.Heart,selection.method = "mvp")
  seurat.carsodo.Heart = RunPCA(seurat.carsodo.Heart,verbose = F,npcs = length(sample)-1)
  g0 = DimPlot(seurat.carsodo.Heart,group.by = "stage")
  sce = as.SingleCellExperiment(seurat.carsodo.Heart)
  library(slingshot)
  library(SingleCellExperiment)
  sce = slingshot(sce,reducedDim = "PCA", allow.breaks = F)
  library(ggplot2)
  tmp = SlingshotDataSet(sce)
  data = data.frame("Stage" = seurat.carsodo.Heart@meta.data$stage,"Pseudotime" = sce$slingPseudotime_1)
  data$Stage = factor(data$Stage, levels = unique(at.cardoso$stage))
  g1 = ggplot(data,aes(y = Stage,x = Pseudotime, color = Stage))+
    geom_point()+theme_bw()
  pdf(sprintf("pseudotime-%s.pdf",tissue))
  print(g0)
  print(g1)
  dev.off()
  if(cor(1:length(sce$slingPseudotime_1),sce$slingPseudotime_1)<0){
    sce$slingPseudotime_1 = -(sce$slingPseudotime_1-max(sce$slingPseudotime_1) )
  }
  scode.1 =  as.matrix(seurat.carsodo.Heart@assays$RNA@data)
  scode.2 = data.frame("id" = 1:ncol(scode.1), "pseudotime" = sce$slingPseudotime_1)
  save(scode.1, scode.2, file = sprintf("scode_input_%s",tissue))
  

  # Infering gene regulatory network with SCODE with D = 5
  library(MASS)
  tfnum <- nrow(scode.1) #gene_number
  cnum <- ncol(scode.1)# sample
  maxite <- 100 # iter
  
  maxB <- 2.0
  minB <- -10.0
  dir = sprintf("SCODE_%s",tissue)
  scode = function(i,pnum,dir=dir){
    RSSs = NULL
    dir2 = sprintf("%s/out%s",dir,i)
    dir.create(dir2)
    X <- as.matrix(scode.1)
    W <- matrix(rep(0,tfnum*pnum), nrow=tfnum, ncol=pnum)
    Z <- matrix(rep(0,pnum*cnum), nrow=pnum, ncol=cnum)
    WZ <- matrix(nrow=tfnum, ncol=cnum)
    
    #read pseudo-time and normalize pseudo-time
    pseudotime <- scode.2$pseudotime
    pseudotime <- pseudotime/max(pseudotime)
    
    new_B <- rep(0, pnum)
    old_B <- rep(0, pnum)
    
    #initialization
    RSS <- Inf
    for(i in 1:pnum){
      new_B[i] <- runif(1, min=minB, max=maxB)
      old_B[i] <- new_B[i]
    }
    
    #function to sample Z
    sample_Z <- function(){
      for(i in 1:pnum){
        for(j in 1:cnum){
          Z[i,j] <<- exp(new_B[i]*pseudotime[j]) + runif(1, min=-0.001, max=0.001)
        }
      }
    }
    
    #optimize W and B iteratively
    for(ite in 1:maxite){
      cat(as.character(Sys.time()), "\t",RSS, "\t", ite,"\t",pnum,"\n")
      #sampling B
      target <- floor(runif(1, min=1, max=pnum+1))
      new_B[target] <- runif(1, min=minB, max=maxB)
      
      #for last calculation
      if(ite == maxite){
        for(i in 1:pnum){
          new_B[i] <- old_B[i]
        }
      }
      
      #sample Z from new B
      sample_Z()
      
      #regression
      for(i in 1:tfnum){
        X.lm <- lm(X[i,] ~ t(Z)-1)
        for(j in 1:pnum){
          W[i,j] <- X.lm$coefficients[j]
        }
        WZ[i,] <- W[i,] %*% Z
      }
      
      #RSS
      tmp_RSS <- sum((X-WZ)**2)
      if(tmp_RSS < RSS){
        RSS <- tmp_RSS
        old_B[target] <- new_B[target]
      }
      else{
        new_B[target] <- old_B[target]
      }
      RSSs = c(RSSs, RSS)
    }
    
    #output RSS
    write.table(RSS, paste(dir2,"/RSS.txt",sep=""), row.names=F, col.names=F, sep="\t")
    
    #output W
    write.table(W, paste(dir2,"/W.txt",sep=""), row.names=F, col.names=F, sep="\t")
    
    #infer A
    B <- matrix(rep(0,pnum*pnum), nrow=pnum, ncol=pnum)
    for(i in 1:pnum){
      B[i,i] <- new_B[i]
    }
    invW <- ginv(W)
    A <- W %*% B %*% invW
    
    #output A and B
    save(A,file = paste(dir2,"/A",sep=""))
    #write.table(A, paste(dir,"/A.txt",sep=""), row.names=F, col.names=F, sep="\t")
    save(B, RSS, file = paste(dir2,"/B-RSS",sep=""))
  }
  dir.create(dir)
  repnum = 20
  registerDoParallel(repnum)
  gc()
  gc()
  gc()
  gc()
  gc()
  gc()
  gc()
  foreach (i=1:repnum) %dopar% scode(i,4,dir)
  
  load(paste(dir,"/out1/A",sep=""))
  meanA <- A
  for(i in 2:repnum){
    cat(i,"\n")
    load(paste(dir,"/out",i,"/A",sep=""))
    tmp <- A
    meanA <- meanA + tmp
    gc(gc())
  }
  meanA <- meanA / repnum
  
  # Check reproducibility of A and RSS
  cors = NULL
  for(a in 1:repnum){
    cat(a,"\n")
    load(paste(dir,"/out",a,"/A",sep=""))
    cors = c(cors,cor(matrix(A,ncol = 1),matrix(meanA,ncol = 1)))
    gc(gc(gc()))
  }
  top10 = order(cors, decreasing = T)[1:10]
  ggplot(data.frame("cor" = cors),aes(x = 1:20,y = cor))+
    geom_point()+
    theme_bw()
  gc(gc())
  
  
  # Calculate average A
  meanA = NULL
  for(i in top10){
    cat(i,"\n")
    if(i == top10[1]){
      load(paste(dir,"/out1/A",sep=""))
      meanA <- A
    } else {
      load(paste(dir,"/out",i,"/A",sep=""))
      tmp <- A
      meanA <- meanA + tmp
    }
    gc(gc())
  }
  meanA <- meanA / repnum
  res.scode = as.matrix(meanA)
  gnl = NULL
  for(gn in rownames(scode.1)){
    tmp = unique(transcript2gene.mouse$GeneName[transcript2gene.mouse$GeneID==gn])
    if(length(tmp)>0){
      names(tmp) = gn
      gnl = c(gnl,tmp)
    }
  }
  colnames(res.scode) = rownames(scode.1)
  rownames(res.scode) = rownames(scode.1)
  res.scode = res.scode[names(gnl),names(gnl)]
  colnames(res.scode) = gnl
  rownames(res.scode) = gnl
  
  save(res.scode, file = sprintf("res-scode-%s",tissue))
  
  # ROC
  load("TF2DNA")
  TF.target =  matrix(0,ncol = ncol(res.scode), nrow = nrow(res.scode))
  colnames(TF.target) = colnames(res.scode)
  rownames(TF.target) = rownames(res.scode)
  for(gn in TF.list){
    TF.target[unique(answer$target_name[answer$tf_name==gn])[is.element(unique(answer$target_name[answer$tf_name==gn]),rownames(res.scode))],gn] = 1
  }
  library(ROCR)
  
  aucs = NULL
  for(gn in TF.list){
    cat(gn,"\n")
    scode = abs(res.scode[,gn])
    answer.TF = TF.target[,gn]
    pred <- prediction(scode, answer.TF)
    perf <- performance(pred, "tpr", "fpr")
    auc.tmp <- performance(pred,"auc")
    auc <- as.numeric(auc.tmp@y.values)
    aucs = c(aucs, auc)
  }
  save(aucs,file = sprintf("AUC-%s",tissue))
}

library(ggbeeswarm)
data = NULL
for(tissue in unique(at.cardoso$tissue)){
  load(sprintf("AUC-%s",tissue))
  data = rbind(data, data.frame("AUC"=aucs,"TF" = TF.list,"tissue" = rep(tissue, length(tissue))))
}

pdf("AUC-cardoso.pdf")
g = ggplot(data, aes(y = AUC, x = 1))+
  geom_beeswarm()+
  facet_wrap(~tissue)+
  theme_bw()
print(g)

# ggplot(data, aes(y = AUC, x =TF))+
#   geom_beeswarm()+
#   geom_vline(xintercept = seq(0,450,50),color = "gray")+
#   geom_beeswarm()+
#   facet_wrap(~tissue,ncol = 1)+
#   theme_bw()
# 
# ggplot(data, aes(y = AUC, x =TF))+
#   geom_beeswarm()+
#   geom_vline(xintercept = seq(0,450,50),color = "gray")+
#   geom_hline(yintercept = seq(0.5,0.9,0.1),color = "red")+
#   geom_beeswarm()+
#   facet_wrap(~tissue,ncol = 7)+
#   theme_bw()

data = NULL
for(tissue in unique(at.cardoso$tissue)){
  load(sprintf("AUC-%s",tissue))
  data = cbind(data, aucs)
}
colnames(data) = unique(at.cardoso$tissue)
cor = cor(data,method = "spearman")
d = dist(cor)
h = hclust(d)
plot(h)
dev.off()

cal.knee = function(totals){
  totals = abs(totals)
  o <- order(totals, decreasing = TRUE)
  stuff <- rle(totals[o])
  run.rank <- cumsum(stuff$lengths) - (stuff$lengths - 1)/2
  run.totals <- stuff$values
  y <- (run.totals)
  x <- (run.rank)
  fit <- nls(y ~ (a * x + b), start = c(a=0, b=0))
  top = which((y - predict(fit,x = x))<0)[1]-1
  knee <- predict(fit)[top]
  return(knee)
}
## ligand recepter
library(RVenn)
for(tissue in unique(at.cardoso$tissue)){
  load(sprintf("res-scode-%s",tissue))
  knees = NULL
  for(gn in colnames(res.scode)){
    cat(gn, "\n")
    if(sd(res.scode[,gn])>0){
      knees = c(knees, cal.knee(res.scode[,gn]))
    } else {
      knees = c(knees, 0)
    }
  }
  names(knees) = rownames(res.scode)  
  save(knees, file = sprintf("knees-%s",tissue))
  
  regulation.mat = matrix(0,ncol = ncol(res.scode), nrow = nrow(res.scode))
  colnames(regulation.mat) = colnames(res.scode)
  rownames(regulation.mat) = rownames(res.scode)
  for(gn in colnames(res.scode)){
    pass = abs(res.scode[,gn])>knees[gn] 
    target = colnames(res.scode)[pass]
    order = order(abs(res.scode[,gn])[target],decreasing = T)
    target = target[order]
    regulation.mat[target,gn] = 1:length(target)
  }
  
  mat = regulation.mat
  ligand.target.table = NULL
  list = read.csv("TableS2.csv")
  for (i in 1:nrow(list)) {
    a = ligand.recepter.num(list[i,],mat)
    ligand.target.table = rbind(ligand.target.table,a$Table)
  }

  data = data.frame("Unique recepter targets" = ligand.target.table$Num.recepter.target-ligand.target.table$Num.overlap,"Common targets" = ligand.target.table$Num.overlap,"Unique ligand targets" = ligand.target.table$Num.ligand.target-ligand.target.table$Num.overlap, "pair" = sprintf("%s-%s",ligand.target.table$Ligand.name,ligand.target.table$Recepter.name))
  library(reshape2)
  data = melt(data)
  pdf(sprintf("ligand-recepter-%s.pdf",tissue))
  g = ggplot(data, aes(x = value, y = pair, fill = variable))+
    geom_bar(stat = "identity", color = "black")+
    scale_fill_manual(values = c("gray40","white","gray"))+
    ggtitle(tissue)+
    theme_bw()
  print(g)
  dev.off()
}

# Microwell-seq
rawcnt.Microwell = read.table("/mnt/hdd/hdd/data/species/M.musculus/mm9/ensembl/Microwell-Seq_mouce/Figure2-batch-removed.txt")
seurat.Microwell = CreateSeuratObject(rawcnt.Microwell)
seurat.Microwell = NormalizeData(seurat.Microwell)
seurat.Microwell =ScaleData(seurat.Microwell)
library(openxlsx)
at = read.xlsx("/mnt/hdd/hdd/data/species/M.musculus/mm9/ensembl/Microwell-Seq_mouce/MCA_Figure2_Cell.Info.xlsx",sheet = 1)
Idents(seurat.Microwell) = at$ClusterID
ex.ave = AverageExpression(seurat.Microwell,slot = "data")
ex.ave = ex.ave$RNA
VlnPlot(seurat.Microwell,features = c("Notch1","Dll1"), ncol = 1)+NoLegend()
load("res-scode")
knees = NULL
for(gn in colnames(res.scode)){
  cat(gn, "\n")
  knees = c(knees, cal.knee(res.scode[,gn]))
}
names(knees) = rownames(res.scode)  

regulation.mat = matrix(0,ncol = ncol(res.scode), nrow = nrow(res.scode))
colnames(regulation.mat) = colnames(res.scode)
rownames(regulation.mat) = rownames(res.scode)
for(gn in colnames(res.scode)){
  pass = abs(res.scode[,gn])>knees[gn] 
  target = colnames(res.scode)[pass]
  order = order(abs(res.scode[,gn])[target],decreasing = T)
  target = target[order]
  regulation.mat[target,gn] = 1:length(target)
}

list = read.csv("TableS2.csv")
regulation.mat = regulation.mat[is.element(rownames(regulation.mat),rownames(rawcnt.Microwell)),]
gnl = c(unlist(lapply(list$Ligand,FUN = function(x){
  x = unlist(strsplit(x, ". "))
  return(x)
})),
unlist(lapply(list$Receptor,FUN = function(x){
  x = unlist(strsplit(x, ". "))
  return(x)
})),
TF.list)
gnl = gnl[is.element(gnl, rownames(rawcnt.Microwell))]
co.ex.target = NULL
for(i in 1:length(gnl)){
  cat(as.character(Sys.time()),"")
  gn = gnl[i]
  cat(sprintf("% 5d",i),"/",length(gnl),gn,"\t")
  targets = rownames(regulation.mat)[regulation.mat[,gn]>0]
  targets = targets[targets !=gn]
  cors = apply(ex.ave[targets,],1,FUN = function(x){
    return(cor(ex.ave[gn,],x))
  })
  # colSums(cors>0.4)
  # a = colnames(rawcnt.Microwell)[rawcnt.Microwell[gn,]>0]
  # b = rowSums(rawcnt.Microwell[targets,a]>0)
  # tmp = data.frame("upstream" = gn,"downstream" = names(b),"up" = length(a), "down" = b)
  # for(target in targets){
  #   # cat(target)
  #   # cat(as.character(Sys.time()),"\n")
  #   b = sum(rawcnt.Microwell[target,a]>0)
  tmp = data.frame("upstream" = gn,"downstream" = targets, "cor" = cors)
  #   # cors = c(cors,cor(seurat.Microwell@assays$RNA@data[gn,],seurat.Microwell@assays$RNA@data[target,]))
  co.ex.target = rbind(co.ex.target,tmp)
  # }
  cat("\n")
  # tmp  = list(cors)
  # names(tmp) = gn
  # # cat(mean(cors),"\n")
  # co.ex.target = c(co.ex.target,tmp)
}

tmp = rep(F, nrow(co.ex.target))
tmp[is.element(co.ex.target$upstream,unlist(lapply(list$Ligand,FUN = function(x){
  x = unlist(strsplit(x, ". "))
  return(x)
})))] = "ligand"

tmp[is.element(co.ex.target$upstream,TF.list)] = "TF"

co.ex.target = cbind(co.ex.target, data.frame("is.ligand" = tmp))
ggplot(co.ex.target, aes(x  = is.ligand, y = ratio.co.ex))+
  geom_boxplot()+
  theme_bw()
mean(co.ex.target$cor[co.ex.target$is.ligand=="TF"])
mean(co.ex.target$cor[co.ex.target$is.ligand=="ligand"])
mean(co.ex.target$cor[co.ex.target$is.ligand==F])
t.test(co.ex.target$ratio.co.ex[co.ex.target$is.ligand=="TF"],co.ex.target$ratio.co.ex[co.ex.target$is.ligand=="ligand"])
t.test(co.ex.target$ratio.co.ex[co.ex.target$is.ligand==F],co.ex.target$ratio.co.ex[co.ex.target$is.ligand=="ligand"])
t.test(co.ex.target$ratio.co.ex[co.ex.target$is.ligand=="TF"],co.ex.target$ratio.co.ex[co.ex.target$is.ligand==F])
library(dplyr)
a = co.ex.target %>% group_by(is.ligand, upstream) %>%縲filter(cor > 0.4) %>% summarise(n())
b = co.ex.target %>% group_by(is.ligand, upstream)  %>% summarise(n())
High.cor.pair = cbind(a,data.frame("High.cor.pair"= a$`n()`/b$`n()`))
library(ggplot2)
pdf("correlated-target-rate-ligand-receptor-TF.pdf")
ggplot(High.cor.pair, aes(x = is.ligand, y = High.cor.pair))+
  geom_boxplot()+
  theme_bw()
dev.off()
t.test(High.cor.pair$High.cor.pair[High.cor.pair$is.ligand==F],High.cor.pair$High.cor.pair[High.cor.pair$is.ligand=="ligand"])
t.test(High.cor.pair$High.cor.pair[High.cor.pair$is.ligand==F],High.cor.pair$High.cor.pair[High.cor.pair$is.ligand=="TF"])
t.test(High.cor.pair$High.cor.pair[High.cor.pair$is.ligand=="TF"],High.cor.pair$High.cor.pair[High.cor.pair$is.ligand=="ligand"])
gn = "Rela"
targets = rownames(regulation.mat)[regulation.mat[,gn]>0]
targets = targets[targets !=gn]
cors = apply(ex.ave[targets,],1,FUN = function(x){
  return(cor(ex.ave[gn,],x))
})
sort(cors)[1:10]
sort(cors, decreasing = T)[1:10]
pdf("correlated-target-rate-example.pdf")
plot(ex.ave[gn,], ex.ave["Slirp",],main = sprintf("Slirp: %s", round(cor(ex.ave[gn,], ex.ave["Slirp",]),3)))
plot(ex.ave[gn,], ex.ave["Rps6ka3",],main = sprintf("Rps6ka3: %s", round(cor(ex.ave[gn,], ex.ave["Rps6ka3",]),3)))
plot(ex.ave[gn,], ex.ave[targets[i],],main = sprintf("%s: %s", targets[i],round(cor(ex.ave[gn,], ex.ave[targets[i],]),3)))
dev.off()


pdf("Microwell-Seq-ligand-recepter-co-expression.pdf")
for(i in 1:nrow(list)){
  tmp = list[i,]
  gnl.l = gsub(" ","",unlist(strsplit(tmp$Ligand,",")))
  gnl.r = gsub(" ","",unlist(strsplit(tmp$Receptor,",")))
  l.posi = colnames(rawcnt.Microwell)[colSums(rawcnt.Microwell[gnl.l,]>0)>0]
  r.posi = colnames(rawcnt.Microwell)[colSums(rawcnt.Microwell[gnl.r,]>0)>0]
  toy = list("Ligand" = l.posi, "Recepter" = r.posi)
  toy = Venn(toy)
  g = ggvenn(toy)+ggtitle(tmp$Signaling.pathway)+theme_void()
  print(g)
}
dev.off()

# Cao 2019
library(monocle3)
cds = readRDS("/mnt/hdd/hdd/data/species/M.musculus/mm9/ensembl/Cao2019/sub_trajectory_summary/Lung_epithelial_trajectory_cds.RDS")
rawcnt = as.matrix(cds@assayData$exprs)
tmp = rownames(rawcnt)
tmp = unlist(lapply(tmp, FUN = function(x){
  x=strsplit(x,"\\.")
return(x[[1]][1])
  }))
tmp2 = NULL
for(gn in tmp){
  tmp3 = unique(transcript2gene.mouse$GeneName[transcript2gene.mouse$GeneID==gn])
  if(length(tmp3)==0){
    tmp3 = gn
  }
  tmp2 = rbind(tmp2,data.frame("ID"=gn,"gn" = tmp3))
}
rownames(rawcnt) = tmp2$gn
seurat = CreateSeuratObject(counts = rawcnt)
seurat$pseudotime = cds$Pseudotime
seurat = NormalizeData(seurat)
seurat = ScaleData(seurat)
seurat = FindVariableFeatures(seurat)
seurat = RunPCA(seurat)
ElbowPlot(seurat)
seurat = RunUMAP(seurat,dims = 1:10)
Idents(seurat) = cds$Cluster
DimPlot(seurat)
FeaturePlot(seurat,features = "pseudotime")
scode.1 =  as.matrix(seurat@assays$RNA@data)
scode.2 = data.frame("id" = 1:ncol(scode.1), "pseudotime" = seurat@meta.data$pseudotime)
colnames(scode.1) = scode.2$id
  
###SCODE
# Infering gene regulatory network with SCODE with D = 5
library(MASS)
tfnum <- nrow(scode.1) #gene_number
cnum <- ncol(scode.1)# sample
maxite <- 100 # iter

maxB <- 2.0
minB <- -10.0
dir = sprintf("SCODE_Cao-Lung")
scode = function(i,pnum,dir=dir){
  RSSs = NULL
  dir2 = sprintf("%s/out%s",dir,i)
  dir.create(dir2)
  X <- as.matrix(scode.1)
  W <- matrix(rep(0,tfnum*pnum), nrow=tfnum, ncol=pnum)
  Z <- matrix(rep(0,pnum*cnum), nrow=pnum, ncol=cnum)
  WZ <- matrix(nrow=tfnum, ncol=cnum)
  
  #read pseudo-time and normalize pseudo-time
  pseudotime <- scode.2$pseudotime
  pseudotime <- pseudotime/max(pseudotime)
  
  new_B <- rep(0, pnum)
  old_B <- rep(0, pnum)
  
  #initialization
  RSS <- Inf
  for(i in 1:pnum){
    new_B[i] <- runif(1, min=minB, max=maxB)
    old_B[i] <- new_B[i]
  }
  
  #function to sample Z
  sample_Z <- function(){
    for(i in 1:pnum){
      for(j in 1:cnum){
        Z[i,j] <<- exp(new_B[i]*pseudotime[j]) + runif(1, min=-0.001, max=0.001)
      }
    }
  }
  
  #optimize W and B iteratively
  for(ite in 1:maxite){
    cat(as.character(Sys.time()), "\t",RSS, "\t", ite,"\t",pnum,"\n")
    #sampling B
    target <- floor(runif(1, min=1, max=pnum+1))
    new_B[target] <- runif(1, min=minB, max=maxB)
    
    #for last calculation
    if(ite == maxite){
      for(i in 1:pnum){
        new_B[i] <- old_B[i]
      }
    }
    
    #sample Z from new B
    sample_Z()
    
    #regression
    for(i in 1:tfnum){
      X.lm <- lm(X[i,] ~ t(Z)-1)
      for(j in 1:pnum){
        W[i,j] <- X.lm$coefficients[j]
      }
      WZ[i,] <- W[i,] %*% Z
    }
    
    #RSS
    tmp_RSS <- sum((X-WZ)**2)
    if(tmp_RSS < RSS){
      RSS <- tmp_RSS
      old_B[target] <- new_B[target]
    }
    else{
      new_B[target] <- old_B[target]
    }
    RSSs = c(RSSs, RSS)
  }
  
  #output RSS
  write.table(RSS, paste(dir2,"/RSS.txt",sep=""), row.names=F, col.names=F, sep="\t")
  
  #output W
  write.table(W, paste(dir2,"/W.txt",sep=""), row.names=F, col.names=F, sep="\t")
  
  #infer A
  B <- matrix(rep(0,pnum*pnum), nrow=pnum, ncol=pnum)
  for(i in 1:pnum){
    B[i,i] <- new_B[i]
  }
  invW <- ginv(W)
  A <- W %*% B %*% invW
  
  #output A and B
  save(A,file = paste(dir2,"/A",sep=""))
  #write.table(A, paste(dir,"/A.txt",sep=""), row.names=F, col.names=F, sep="\t")
  save(B, RSS, file = paste(dir2,"/B-RSS",sep=""))
}
dir.create(dir)
repnum = 16
library(doParallel)
registerDoParallel(repnum)
foreach (i=1:repnum) %dopar% scode(i,4,dir)

load(paste(dir,"/out1/A",sep=""))
meanA <- A
for(i in 2:repnum){
  cat(i,"\n")
  load(paste(dir,"/out",i,"/A",sep=""))
  tmp <- A
  meanA <- meanA + tmp
  gc(gc())
}
meanA <- meanA / repnum

# Check reproducibility of A and RSS
cors = NULL
for(a in 1:repnum){
  cat(a,"\n")
  load(paste(dir,"/out",a,"/A",sep=""))
  cors = c(cors,cor(matrix(A,ncol = 1),matrix(meanA,ncol = 1)))
  gc(gc(gc()))
}
top10 = order(cors, decreasing = T)[1:10]
ggplot(data.frame("cor" = cors),aes(x = 1:20,y = cor))+
  geom_point()+
  theme_bw()
gc(gc())


# Calculate average A
meanA = NULL
for(i in top10){
  cat(i,"\n")
  if(i == top10[1]){
    load(paste(dir,"/out1/A",sep=""))
    meanA <- A
  } else {
    load(paste(dir,"/out",i,"/A",sep=""))
    tmp <- A
    meanA <- meanA + tmp
  }
  gc(gc())
}
meanA <- meanA / 10
res.scode = as.matrix(meanA)
# gnl = NULL
# for(gn in rownames(scode.1)){
#   tmp = unique(transcript2gene.mouse$GeneName[transcript2gene.mouse$GeneID==gn])
#   if(length(tmp)>0){
#     names(tmp) = gn
#     gnl = c(gnl,tmp)
#   }
# }
colnames(res.scode) = rownames(scode.1)
rownames(res.scode) = rownames(scode.1)
# res.scode = res.scode[names(gnl),names(gnl)]
save(res.scode, file = sprintf("res-scode-Cao-Lung"))
# ROC
load("TF2DNA")
TF.target =  matrix(0,ncol = ncol(res.scode), nrow = nrow(res.scode))
colnames(TF.target) = colnames(res.scode)
rownames(TF.target) = rownames(res.scode)
for(gn in TF.list[is.element(TF.list, rownames(res.scode))]){
  TF.target[unique(answer$target_name[answer$tf_name==gn])[is.element(unique(answer$target_name[answer$tf_name==gn]),rownames(res.scode))],gn] = 1
}
library(ROCR)

aucs = NULL
for(gn in TF.list[is.element(TF.list, rownames(res.scode))]){
  cat(gn,"\n")
  scode = abs(res.scode[,gn])
  answer.TF = TF.target[,gn]
  pred <- prediction(scode, answer.TF)
  perf <- performance(pred, "tpr", "fpr")
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  aucs = c(aucs, auc)
}
save(aucs,file = sprintf("AUC-Cao-Lung"))
load("AUC-Cao-Lung")
g1 = ggplot(data.frame("AUC"=aucs), aes(y = AUC, x =1))+
  geom_beeswarm()+
  ggtitle("Lung lineage")+
  theme_bw()
load("AUC-Cao-Endocardium")
g2 = ggplot(data.frame("AUC"=aucs), aes(y = AUC, x =1))+
  geom_beeswarm()+
  ggtitle("Endocardium lineage")+
  theme_bw()
pdf("AUC-single-cell.pdf")
print(g1)
print(g2)
dev.off()

cal.knee = function(totals){
  totals = abs(totals)
  o <- order(totals, decreasing = TRUE)
  stuff <- rle(totals[o])
  run.rank <- cumsum(stuff$lengths) - (stuff$lengths - 1)/2
  run.totals <- stuff$values
  y <- (run.totals)
  x <- (run.rank)
  fit <- nls(y ~ (a * x + b), start = c(a=0, b=0))
  top = which((y - predict(fit,x = x))<0)[1]-1
  knee <- predict(fit)[top]
  return(knee)
}
## ligand recepter
library(RVenn)
knees = NULL
for(gn in colnames(res.scode)){
  cat(gn, "\n")
  if(sd(res.scode[,gn])>0){
    knees = c(knees, cal.knee(res.scode[,gn]))
  } else {
    knees = c(knees, 0)
  }
}
names(knees) = rownames(res.scode)  

regulation.mat = matrix(0,ncol = ncol(res.scode), nrow = nrow(res.scode))
colnames(regulation.mat) = colnames(res.scode)
rownames(regulation.mat) = rownames(res.scode)
for(gn in colnames(res.scode)){
  pass = abs(res.scode[,gn])>knees[gn] 
  target = colnames(res.scode)[pass]
  order = order(abs(res.scode[,gn])[target],decreasing = T)
  target = target[order]
  regulation.mat[target,gn] = 1:length(target)
}

mat = regulation.mat
ligand.target.table = NULL
list = read.csv("TableS2.csv")
for (i in 1:nrow(list)) {
  a = ligand.recepter.num(list[i,],mat)
  ligand.target.table = rbind(ligand.target.table,a$Table)
}

data = data.frame("Unique recepter targets" = ligand.target.table$Num.recepter.target-ligand.target.table$Num.overlap,"Common targets" = ligand.target.table$Num.overlap,"Unique ligand targets" = ligand.target.table$Num.ligand.target-ligand.target.table$Num.overlap, "pair" = sprintf("%s-%s",ligand.target.table$Ligand.name,ligand.target.table$Recepter.name))
library(reshape2)
data = melt(data)
pdf(sprintf("ligand-recepter-Cao-Lung.pdf"))
g = ggplot(data, aes(x = value, y = pair, fill = variable))+
  geom_bar(stat = "identity", color = "black")+
  scale_fill_manual(values = c("gray40","white","gray"))+
  ggtitle("Cao-Lung")+
  theme_bw()
print(g)
dev.off()

#################################################################
library(monocle3)
cds = readRDS("/mnt/hdd/hdd/data/species/M.musculus/mm9/ensembl/Cao2019/sub_trajectory_summary/Endocardium_trajectory_cds.RDS")
rawcnt = as.matrix(cds@assayData$exprs)
tmp = rownames(rawcnt)
tmp = unlist(lapply(tmp, FUN = function(x){
  x=strsplit(x,"\\.")
  return(x[[1]][1])
}))
tmp2 = NULL
for(gn in tmp){
  tmp3 = unique(transcript2gene.mouse$GeneName[transcript2gene.mouse$GeneID==gn])
  if(length(tmp3)==0){
    tmp3 = gn
  }
  tmp2 = rbind(tmp2,data.frame("ID"=gn,"gn" = tmp3))
}
rownames(rawcnt) = tmp2$gn
library(Seurat)
seurat = CreateSeuratObject(counts = rawcnt)
seurat$pseudotime = cds$Pseudotime
seurat = NormalizeData(seurat)
seurat = ScaleData(seurat)
seurat = FindVariableFeatures(seurat)
seurat = RunPCA(seurat)
ElbowPlot(seurat)
seurat = RunUMAP(seurat,dims = 1:10)
Idents(seurat) = cds$Cluster
DimPlot(seurat)
plot(seurat@meta.data$nCount_RNA)
FeaturePlot(seurat,features = "pseudotime")
scode.1 =  as.matrix(seurat@assays$RNA@data)
scode.2 = data.frame("id" = 1:ncol(scode.1), "pseudotime" = seurat@meta.data$pseudotime)
colnames(scode.1) = scode.2$id

###SCODE
# Infering gene regulatory network with SCODE with D = 5
library(MASS)
tfnum <- nrow(scode.1) #gene_number
cnum <- ncol(scode.1)# sample
maxite <- 100 # iter

dir = sprintf("SCODE_Cao-Endocardium")
dir.create(dir)
repnum = 16
library(doParallel)
registerDoParallel(repnum)
foreach (i=1:repnum) %dopar% scode(i,4,dir)

load(paste(dir,"/out1/A",sep=""))
meanA <- A
for(i in 2:repnum){
  cat(i,"\n")
  load(paste(dir,"/out",i,"/A",sep=""))
  tmp <- A
  meanA <- meanA + tmp
  gc(gc())
}
meanA <- meanA / repnum

# Check reproducibility of A and RSS
cors = NULL
for(a in 1:repnum){
  cat(a,"\n")
  load(paste(dir,"/out",a,"/A",sep=""))
  cors = c(cors,cor(matrix(A,ncol = 1),matrix(meanA,ncol = 1)))
  gc(gc(gc()))
}
top10 = order(cors, decreasing = T)[1:10]
ggplot(data.frame("cor" = cors),aes(x = 1:20,y = cor))+
  geom_point()+
  theme_bw()
gc(gc())


# Calculate average A
meanA = NULL
for(i in top10){
  cat(i,"\n")
  if(i == top10[1]){
    load(paste(dir,"/out1/A",sep=""))
    meanA <- A
  } else {
    load(paste(dir,"/out",i,"/A",sep=""))
    tmp <- A
    meanA <- meanA + tmp
  }
  gc(gc())
}
meanA <- meanA / 10
res.scode = as.matrix(meanA)
# gnl = NULL
# for(gn in rownames(scode.1)){
#   tmp = unique(transcript2gene.mouse$GeneName[transcript2gene.mouse$GeneID==gn])
#   if(length(tmp)>0){
#     names(tmp) = gn
#     gnl = c(gnl,tmp)
#   }
# }
colnames(res.scode) = rownames(scode.1)
rownames(res.scode) = rownames(scode.1)
# res.scode = res.scode[names(gnl),names(gnl)]
save(res.scode, file = sprintf("res-scode-Cao-Endocardium"))
# ROC
load("TF2DNA")
TF.target =  matrix(0,ncol = ncol(res.scode), nrow = nrow(res.scode))
colnames(TF.target) = colnames(res.scode)
rownames(TF.target) = rownames(res.scode)
for(gn in TF.list[is.element(TF.list, rownames(res.scode))]){
  TF.target[unique(answer$target_name[answer$tf_name==gn])[is.element(unique(answer$target_name[answer$tf_name==gn]),rownames(res.scode))],gn] = 1
}
library(ROCR)

aucs = NULL
for(gn in TF.list[is.element(TF.list, rownames(res.scode))]){
  cat(gn,"\n")
  scode = abs(res.scode[,gn])
  answer.TF = TF.target[,gn]
  pred <- prediction(scode, answer.TF)
  perf <- performance(pred, "tpr", "fpr")
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  aucs = c(aucs, auc)
}
save(aucs,file = sprintf("AUC-Cao-Endocardium"))
ggplot(data.frame("AUC"=aucs), aes(y = AUC, x =1))+
  geom_beeswarm()+
  theme_bw()
縲皇al.knee = function(totals){
  totals = abs(totals)
  o <- order(totals, decreasing = TRUE)
  stuff <- rle(totals[o])
  run.rank <- cumsum(stuff$lengths) - (stuff$lengths - 1)/2
  run.totals <- stuff$values
  y <- (run.totals)
  x <- (run.rank)
  fit <- nls(y ~ (a * x + b), start = c(a=0, b=0))
  top = which((y - predict(fit,x = x))<0)[1]-1
  knee <- predict(fit)[top]
  return(knee)
}
## ligand recepter
library(RVenn)
knees = NULL
for(gn in colnames(res.scode)){
  cat(gn, "\n")
  if(sd(res.scode[,gn])>0){
    knees = c(knees, cal.knee(res.scode[,gn]))
  } else {
    knees = c(knees, 0)
  }
}
names(knees) = rownames(res.scode)  

regulation.mat = matrix(0,ncol = ncol(res.scode), nrow = nrow(res.scode))
colnames(regulation.mat) = colnames(res.scode)
rownames(regulation.mat) = rownames(res.scode)
for(gn in colnames(res.scode)){
  pass = abs(res.scode[,gn])>knees[gn] 
  target = colnames(res.scode)[pass]
  order = order(abs(res.scode[,gn])[target],decreasing = T)
  target = target[order]
  regulation.mat[target,gn] = 1:length(target)
}

mat = regulation.mat
ligand.target.table = NULL
list = read.csv("TableS2.csv")
for (i in 1:nrow(list)) {
  a = ligand.recepter.num(list[i,],mat)
  ligand.target.table = rbind(ligand.target.table,a$Table)
}

data = data.frame("Unique recepter targets" = ligand.target.table$Num.recepter.target-ligand.target.table$Num.overlap,"Common targets" = ligand.target.table$Num.overlap,"Unique ligand targets" = ligand.target.table$Num.ligand.target-ligand.target.table$Num.overlap, "pair" = sprintf("%s-%s",ligand.target.table$Ligand.name,ligand.target.table$Recepter.name))
library(reshape2)
data = melt(data)
pdf(sprintf("ligand-recepter-Cao-Lung.pdf"))
g = ggplot(data, aes(x = value, y = pair, fill = variable))+
  geom_bar(stat = "identity", color = "black")+
  scale_fill_manual(values = c("gray40","white","gray"))+
  ggtitle("Cao-Lung")+
  theme_bw()
print(g)
dev.off()

