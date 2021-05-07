# Load sample attribute
library(openxlsx)
at.mouse = read.xlsx("sample_info.xlsx",sheet = 1)
at.mouse = at.mouse[!is.na(at.mouse$name),]
tmp = unlist(lapply(strsplit(at.mouse$name, " "),FUN = function(x){return(x[1])}))
at.mouse = cbind(at.mouse, data.frame("stage" = tmp))
at.mouse$stage = factor(at.mouse$stage,levels = c("E7.5","E8.5","E9.5","E10.5","E11.5","E12.5","E13.5"))

# Load count data
# RNA-Seq results are deposied in PRJNA725414 of SRA
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
seurat.mouse = FindClusters(seurat.mouse,verbose = F)
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
pdf("Fig1AB.pdf")
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

#optimize D for Fig. 2A
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

# Infering gene regulatory network with SCODE with D = 4
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
foreach (i=1:repnum) %dopar% scode(i,4)

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

# Check reproducibility of A and RSS for Fig. 2B
dir = "SCODEmouse_out_all"
cors = NULL
for(a in 1:repnum){
  cat(a,"\n")
  load(paste(dir,"/out",a,"/A",sep=""))
  cors = c(cors,cor(matrix(A,ncol = 1),matrix(meanA,ncol = 1)))
  gc(gc(gc()))
}
pdf("Fig.2B.pdf", height = 3.5)
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

# Check correlation of A and expression levels for Fig. 2D
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

pdf("Fig2D.pdf",height = 3.5)
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
pdf("Fig2EF.pdf",height = 5)
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
pdf("Fig2G.pdf", width = 9)
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

# Validation with TF2DNA database (http://fiserlab.org/pscan_files.tar.gz) for Fig. 3
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

pdf("Fig3D.pdf")
grid.arrange(g1, g4, g3, g2, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()
data = data.frame("TF" = correct.rate$gene.name,"Total" = correct.rate$Num.of.predicted.target,"Positive" = correct.rate$Num.of.predicted.target*correct.rate$correct.rate, "Negative" = correct.rate$Num.of.predicted.target-(correct.rate$Num.of.predicted.target*correct.rate$correct.rate))
data = data[base::order(data$Total, decreasing = T),]
library(reshape2)
data = melt(data, id.vars = c("TF","Total"))
data$TF = factor(data$TF, levels = as.character(unique(data$TF)))
data$variable = factor(data$variable, levels = c("Negative","Positive"))
pdf("Fig3C.pdf")
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
pdf("Fig3B.pdf",height = 3.5)
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
pdf("Fig3A.pdf",width = 3.5)
ggplot(data.frame("AUC" = aucs), aes(y = AUC, x = 1))+
  geom_beeswarm()+
  theme_bw()
dev.off()
which(aucs == max(aucs))
which(aucs == min(aucs))

# Overlapping analysis for ligand-receptor pairs for Fig. 4
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
for (i in 1:nrow(list)) {
  a = ligand.recepter.num(list[i,],mat)
  ligand.target.table = rbind(ligand.target.table,a$Table)
}

write.csv(correct.rate, file = "TableS1.csv")
write.csv(ligand.target.table, file = "ligand-target.csv")

data = data.frame("Unique recepter targets" = ligand.target.table$Num.recepter.target-ligand.target.table$Num.overlap,"Common targets" = ligand.target.table$Num.overlap,"Unique ligand targets" = ligand.target.table$Num.ligand.target-ligand.target.table$Num.overlap, "pair" = sprintf("%s-%s",ligand.target.table$Ligand.name,ligand.target.table$Recepter.name))
library(reshape2)
data = melt(data)
pdf("Fig4A.pdf")
ggplot(data, aes(x = value, y = pair, fill = variable))+
  geom_bar(stat = "identity", color = "black")+
  scale_fill_manual(values = c("gray40","white","gray"))+
  theme_bw()
dev.off()

# Top ligand and receotor
library(UpSetR)
pdf("Fig4BCD_S1.pdf", height = 3.5)
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



pdf("FigS3.pdf", height = 3.5)
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

pdf("FigS5.pdf", width = 14, height = 21)
VlnPlot(seurat.mouse, features = c("Wnt1", "Wnt2" , "Wnt2b" , "Wnt3"  , "Wnt3a", "Wnt4" , "Wnt5a" , "Wnt5b", "Wnt6" ,   "Wnt7a",   "Wnt7b" , "Wnt8a","Wnt8b", "Wnt9a" , "Wnt9b" ,  "Wnt10a","Wnt10b",      "Wnt11"  ,  "Wnt16" ),group.by = "stage", ncol = 2)
dev.off()

pdf("FigS4.pdf", width = 14, height = 10.5)
VlnPlot(seurat.mouse, features = c( "Fzd1", "Fzd2","Fzd3" ,"Fzd4","Fzd5", "Fzd6",  "Fzd7"  ,  "Fzd8","Fzd9","Fzd10"      ),group.by = "stage", ncol = 2)
dev.off()

pdf("Fig4D.pdf", height = 3.5)
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

