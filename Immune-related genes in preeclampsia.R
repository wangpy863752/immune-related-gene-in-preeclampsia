#PCA analysis
batchPCA =function(indata, batch, fig.dir, PCA.fig.title, pos="bottomright", xy=c(1,2), cols=NULL, showID=FALSE, cex=1, showLegend=T) {
# indata is a data matrix with samples in columns and genes in rows.
# batch is a vector with the order matching the order in indata.
    library(ClassDiscovery)
    
    outfile = file.path(fig.dir, paste(PCA.fig.title, ".pdf",sep=""))
    N.batch = length(unique(batch))    
    if (is.null(cols)) { 
      cols <- rainbow(N.batch) 
    }else{
      if (length(cols) != N.batch) {stop("cols length not equal to batch length")}
    }           
    
    indata=na.omit(indata)
    pca<-SamplePCA(indata, usecor=F, center=T)
    pct1 <- round (pca@variances[xy[1]]/sum(pca@variances), digits=3)*100
    pct2 <- round (pca@variances[xy[2]]/sum(pca@variances), digits=3)*100
    xlab.text = paste("Comp ", xy[1], ": ", as.character(pct1), "% variance", sep="")
    ylab.text = paste("Comp ", xy[2], ": ", as.character(pct2), "% variance", sep="")    
    
    pdf(file=outfile)
    plot(pca@scores[,xy[1]], pca@scores[,xy[2]],  cex=0.7, xlab=xlab.text, ylab=ylab.text, col=cols[factor(batch)], pch=(1:N.batch)[factor(batch)],lwd=1.5, main=PCA.fig.title)
    abline(h=0, v=0, col="brown", lty=2)
    abline(h=0, v=0, col="brown", lty=2)
    center1<-tapply(pca@scores[,xy[1]], factor(batch), mean)
    center2<-tapply(pca@scores[,xy[2]], factor(batch), mean)
    for (ii in 1:length(center1)) {
        groupi<-pca@scores[as.numeric(factor(batch))==ii, xy]
        #  print(paste("Cluster", ii))
        if (class(groupi)=="matrix") {
            for (j in (1:nrow(groupi))) {
                segments( groupi[j,1], groupi[j,2], center1[ii], center2[ii], col=cols[ii] , lwd=0.3)
            }
        }else {
            segments( groupi[1], groupi[2], center1[ii], center2[ii], col=cols[ii] , lwd=0.3)
        }
    }
    points(center1, center2, pch=7, lwd=1.5,col=cols)
    if (showID) {
      text(pca@scores[,xy[1]], pca@scores[,xy[2]], colnames(indata), lwd=1, cex=cex)
    }
    if(showLegend){
      legend(pos,legend=names(table(factor(batch))), text.col=cols, pch=(1:N.batch), col=cols, lty=1)
    }
    invisible(dev.off())
}
library(sva)
library(cluster)
library(oompaBase)

Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
source("batchPCA.R")
blue <- "#2874C5"
yellow <- "#EABF00"
green <- "#008B8A"
GSE10588<-read.table("GSE10588.txt",head=T,sep='\t',check.names = F,row.names = 1)
GSE25906<-read.table("GSE25906.txt",head=T,sep='\t',check.names = F,row.names = 1)
GSE48424<-read.table("GSE48424.txt",head=T,sep='\t',check.names = F,row.names = 1)
comgene <- intersect(intersect(rownames(GSE10588), rownames(GSE25906)), rownames(GSE48424))
combined.expr <- cbind.data.frame(GSE10588[comgene,],
                                  GSE25906[comgene,],
                                  GSE48424[comgene,])
batchPCA(indata = t(scale(t(combined.expr))),
         batch = rep(c("GSE10588","GSE25906","GSE48424"), times = c(ncol(GSE10588),ncol(GSE25906),ncol(GSE48424))),
         fig.dir = ".",
         PCA.fig.title = "Raw PCA for combined expression profile",
         cols = c(blue, yellow, green),
         showID = F,
         cex = 0.7,
         showLegend = T)

batch <- data.frame(batch = rep(c("GSE10588","GSE25906","GSE48424"), times = c(ncol(GSE10588),ncol(GSE25906),ncol(GSE48424))))
modcombat = model.matrix(~1, data = batch)
combined.expr.combat <- as.data.frame(ComBat(dat=as.matrix(combined.expr), batch=batch$batch, mod=modcombat))
write.table(combined.expr.combat, "merge.txt", quote = F,sep='\t')
batchPCA(indata = t(scale(t(combined.expr.combat))),
         batch = rep(c("GSE10588","GSE25906","GSE48424"), times = c(ncol(GSE10588),ncol(GSE25906),ncol(GSE48424))),
         fig.dir = ".",
         PCA.fig.title = "Combat PCA for combined expression profile",
         cols = c(blue, yellow, green),
         showID = F,
         cex = 0.7,
         showLegend = T)

#Diffrentiall expression analysis
logFoldChange=0                #logFC过滤阈值
adjustP=0.05                   #矫正后p值阈值
conNum=81                     #control组样品数目
treatNum=58                  #treat组样品数目

library(limma)
rt=read.table("merge.txt",sep="\t",header=T,check.names=F)    #读取输入文件

#differential
modType=c(rep("con",conNum),rep("treat",treatNum))
design <- model.matrix(~0+factor(modType))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
write.table(allDiff,file="mrnaAll.xls",sep="\t",quote=F)

#write table
diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & adj.P.Val < adjustP )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file="mrnaDiff.xls",sep="\t",quote=F,col.names=F)
write.table(diffSigOut,file="mrnaDiff.txt",sep="\t",quote=F,col.names=F)

#绘制火山图
pdf(file="mrnaVol.pdf",height=5,width=5)
xMax=max(abs(allDiff$logFC))
yMax=max(-log10(allDiff$adj.P.Val))
plot(allDiff$logFC, -log10(allDiff$adj.P.Val), xlab="logFC",ylab="-log10(adj.P.Val)",
     main="Volcano", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=0.8)
diffSub=subset(allDiff, adj.P.Val<adjustP & logFC>logFoldChange)
points(diffSub$logFC, -log10(diffSub$adj.P.Val), pch=20, col="red",cex=0.8)
diffSub=subset(allDiff, adj.P.Val<adjustP & logFC<(-logFoldChange))
points(diffSub$logFC, -log10(diffSub$adj.P.Val), pch=20, col="green",cex=0.8)
abline(v=0,lty=2,lwd=3)
dev.off()

#绘制差异基因热图
library(pheatmap)
hmExp=rt[rownames(diffSig),]
Type=c(rep("N",conNum),rep("T",treatNum))
names(Type)=colnames(rt)
Type=as.data.frame(Type)
pdf(file="mrnaHeatmap.pdf",height=12,width=15)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("purple", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,show_rownames=F,
         scale="row",
         fontsize = 12,
         fontsize_row=3,
         fontsize_col=10)
dev.off()

#Logsitic analysis
library(SimDesign)  

Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

dat <- read.table("logistic_input.txt",row.names = 1,sep = "\t",header = T,check.names = F,stringsAsFactors = F)

table(dat$Group)
colnames(dat) <- gsub("-","_",colnames(dat))

rname <- "Group" 
vname <- setdiff(colnames(dat), rname) 


or <- p <- p.lab <- c() 
p.cutoff <- 0.05 
step.dir <- "forward" 

#Univariate regression analysis
for (v in vname) {
  f <- as.formula(paste0(rname,"~",v)) 
  t <- dat[,c(rname,v)]
  

  l <- glm(f, 
           data = t, 
           family = "binomial", 
           control = list(maxit = 50), 
           na.action = na.exclude)
  
  s <- format(round(exp(cbind("OR" = coef(l), confint.default(l, level = 0.95)))[2,],3),nsmall = 3) 
  s <- paste0(s[1]," (",s[2],"-",s[3],")")
  or <- c(or,s) # odd ratio
  p <- c(p,format(round(summary(l)$coefficients[2,4],3),nsmall = 3))
  p.lab <- c(p.lab, 
             ifelse(summary(l)$coefficients[2,4] < 0.001,
                    "<0.001", format(round(summary(l)$coefficients[2,4],3),nsmall = 3)))
}

vname.sig <- vname[which(as.numeric(p) < p.cutoff)]
#Multivariate regression analysis
if(length(vname.sig) == 0) {
  cat("No significant variable found!\n") 
} else if(length(vname.sig) == 1) {
  cat("Only one significant variable found!\n") 
} else {
  cat(paste0("A total of ",length(vname.sig)," significant variables found!\n"))
  f <- as.formula(paste0(rname,"~",paste0(vname.sig,collapse = " + "))) 
  t <- dat[,c(rname,vname.sig)]
  l <- glm(f, 
           data = t, 
           family = "binomial", 
           control = list(maxit = 50), 
           na.action = na.exclude)
  l.step <- quiet(step(l,direction = step.dir, k = qchisq(p.cutoff/2,1,lower.tail=FALSE))) 
  l.step.s <- as.data.frame(format(round(exp(cbind("OR" = coef(l.step), 
                                     confint.default(l.step, level = 0.95))),3),nsmall = 3)) 
  l.step.s$p <- format(round(summary(l.step)$coefficients[,4],3),nsmall = 3)
  l.step.s$p.lab <- ifelse(summary(l.step)$coefficients[,4] < 0.001,
                       "<0.001",format(round(summary(l.step)$coefficients[,4],3),nsmall = 3))
  
  l.step.s <- l.step.s[setdiff(rownames(l.step.s),"(Intercept)"),] 
}

step.p <- ifelse(l.step.s$p < p.cutoff, l.step.s$p,"")
step.or <- ifelse(step.p == "","NA",paste0(l.step.s$OR," (",l.step.s$`2.5 %`,"-",l.step.s$`97.5 %`,")"))


step.p.lab <- rep("", length(vname)); names(step.p.lab) <- vname; step.p.lab[vname.sig] <- step.p
step.or.lab <- rep("NA", length(vname)); names(step.or.lab) <- vname; step.or.lab[vname.sig] <- step.or

outTab <- cbind.data.frame(vname, or, p.lab, step.or.lab, step.p.lab)
colnames(outTab) <- c("","Univariate analysis\nOR (95% CI)","\nP value","Multivariate analysis\nOR (95% CI)","\nP value")
rownames(outTab) <- NULL
write.table(outTab,"batch logistic results.txt",sep = "\t",row.names = F,quote = F)

#Construction of nomogram

library(rms)
library(rmda)
rt<-read.table("logistic_input.txt",head=T,sep='\t',check.names = F,row.names = 1)
ddist=datadist(rt)
options(datadist="ddist")


lrmModel=lrm(Group~ CCL24+GNAI1+LCP2+ENG+FLT3, data=rt, x=T, y=T)
nomo=nomogram(lrmModel, fun=plogis,
	fun.at=c(0.0001,0.1,0.3,0.5,0.7,0.9,0.99),
	lp=F, funlabel="Risk of Disease")

pdf("Nomo.pdf", width=8, height=6)
plot(nomo)
dev.off()


cali=calibrate(lrmModel, method="boot", B=1000)
pdf("Calibration.pdf", width=5.5, height=5.5)
plot(cali,
	xlab="Predicted probability",
	ylab="Actual probability", sub=F)
dev.off()


rt$Type=ifelse(rt$Type=="Control", 0, 1)
dc=decision_curve(Type ~ CCL24+GNAI1+LCP2+ENG+FLT3, data=rt, 
	family = binomial(link ='logit'),
	thresholds= seq(0,1,by = 0.01),
	confidence.intervals = 0.95)

pdf(file="DCA.pdf", width=5.5, height=5.5)
plot_decision_curve(dc,
	curve.names="Model",
	xlab="Threshold probability",
	cost.benefit.axis=T,
	col="red",
	confidence.intervals=FALSE,
	standardize=FALSE)
dev.off()

#Molecular subtype identification
data<-read.table("Consensus_input.txt",head=T,sep='\t',check.names = F,row.names = 1)
data<-as.matrix(data)
maxK=9
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(data,
              maxK=maxK,
              reps=50,
              pItem=0.8,
              pFeature=1,
              title=workDir,
              clusterAlg="km",
              distance="euclidean",
              seed=123456,
              plot="png")

Kvec = 2:maxK
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec))
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)}#end for i# The optimal K
optK = Kvec[which.min(PAC)]

clusterNum=optK                  
cluster=results[[clusterNum]][["consensusClass"]]
write.table(cluster,file="cluster.txt",sep="\t",quote=F,col.names=F)
