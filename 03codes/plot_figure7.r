library(Ropt)
library(data.table)
exp=fread("fpkm.allsamples.add2nd.txt", header=T,data.table=F)
exp=exp[apply(exp[,-1],1,mean)>=1,]
g=fread("anthocyanin_pathway.txt", header=T,data.table=F)[,1:2]
lig=fread("lignin_syn_reduced.txt", header=T,data.table=F)
x=fread("MPH_ratio_summary_with_p3_Yunhui.txt", header=T,data.table=F)
pos=x[x[,7]=="positive",]
neg=x[x[,7]=="negative",]
d=fread("R1_MPH_target_Gene.txt", header=T,data.table=F)
d1=fread("R1_interacted_gene_inb.txt", header=T,data.table=F)
d2=fread("R1_interacted_gene_hyb.txt", header=T,data.table=F)
d1=unique(c(d1[,1],d1[,2]))
d1=d1[d1!="Zm00001d026147"]
d2=unique(c(d2[,1],d2[,2]))
d2=d2[d2!="Zm00001d026147"]
re=fread("ZmR1.gene.FPKM_MPH_V2.txt", header=T,data.table=F)
re=re[re[,1]%in%exp[,1],]
re1=re[,c(1:6,9,8,7,12,11,10,13,16,14,17,15,18,19,21,20)]
re1[,-1]=apply(re1[,-1],2,as.numeric)
re1[,-1]=apply(re1[,-1],2,function(x) {x*100})
MEAA=apply(re1[,c(2,6)],1,mean)
MEaa=apply(re1[,c(3,7)],1,mean)

JJ=apply(re1[,c(10,14)],1,mean)
JJO1=apply(re1[,c(11,15)],1,mean)
JJO2=apply(re1[,c(12,16)],1,mean)
JJO3=apply(re1[,c(13,17)],1,mean)
re2=data.frame(re1[,1],JJ,JJO1,JJO2,JJO3,MEAA,MEaa)
colnames(re2)[1]="Gene"

re3=re2[re2[,1]%in%d[,1],]
re4=re2[re2[,1]%in%g[,2],]
re5=re2[re2[,1]%in%lig[,1],]
ov=intersect(g[,2],lig[,1])
re6=re2[re2[,1]%in%ov,]

pdf("OE_EMS_hyb_MPH_V6.pdf",width=8.6,height=2)
par(mfrow=c(1,4),mar=c(2,2,2,2))
###EMS
plot(re2[,6],re2[,7],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="MPH_EMS_AA",ylab="EMS_aa",
     xlim=c(-100,300),ylim=c(-100,300),las=1,tck=-0.03)

cor.test(re2[,6],re2[,7])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))
points(re3[,6],re3[,7],pch=16,cex=0.7,col=adjustcolor("darkred",alpha.f = 0.5))
points(re4[,6],re4[,7],pch=16,cex=0.7,col=adjustcolor("purple",alpha.f = 0.4))
points(re5[,6],re5[,7],pch=16,cex=0.7,col=adjustcolor("darkgreen",alpha.f = 0.4))
points(re6[,6],re6[,7],pch=16,cex=0.7,col=adjustcolor("orange",alpha.f = 0.8))
g2=re2[re2[,1]=="Zm00001d026147",]
points(g2[,6],g2[,7],pch=16,cex=1,col="red")
text(g2[,6],g2[,7],"ZmR1",cex=0.8,pos=4,font=3)
#legend("topleft",legend=c("bHLH targets","Anthocyanin","Lignin","Shared"),col=adjustcolor(c("darkred","purple","darkgreen","orange"),0.6),pch=16,cex=0.7,bty="n")


####R1 targets OE1
plot(re2[,2],re2[,3],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="MPH_J92/JING724",ylab="MPH_J92/OE1",
     xlim=c(-100,300),ylim=c(-100,300),las=1,tck=-0.03)
points(re3[,2],re3[,3],pch=16,cex=0.7,col=adjustcolor("darkred",alpha.f = 0.5))
cor.test(re2[,2],re2[,3])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))

points(re4[,2],re4[,3],pch=16,cex=0.7,col=adjustcolor("purple",alpha.f = 0.4))
points(re5[,2],re5[,3],pch=16,cex=0.7,col=adjustcolor("darkgreen",alpha.f = 0.4))
points(re6[,2],re6[,3],pch=16,cex=0.7,col=adjustcolor("orange",alpha.f = 0.8))
g2=re2[re2[,1]=="Zm00001d026147",]

x1=re5[,3]-re5[,2]
n=length(x1[x1>0])/length(x1)
x2=re5[,4]-re5[,2]
n2=length(x2[x2>0])/length(x2)
x3=re5[,5]-re5[,2]
n3=length(x3[x3>0])/length(x3)
x4=re5[,7]-re5[,6]
n4=length(x4[x4>0])/length(x4)
points(g2[,2],g2[,3],pch=16,cex=1,col="red")
text(g2[,2],g2[,3],"ZmR1",cex=0.8,pos=4,font=3)

####R1 targets OE2
plot(re2[,2],re2[,4],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="MPH_J92/JING724",ylab="MPH_J92/OE2",
     xlim=c(-100,300),ylim=c(-100,300),las=1,tck=-0.03)

cor.test(re2[,2],re2[,4])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))

points(re3[,2],re3[,4],pch=16,cex=0.7,col=adjustcolor("darkred",alpha.f = 0.5))
points(re4[,2],re4[,4],pch=16,cex=0.7,col=adjustcolor("purple",alpha.f = 0.4))
points(re5[,2],re5[,4],pch=16,cex=0.7,col=adjustcolor("darkgreen",alpha.f = 0.4))
points(re6[,2],re6[,4],pch=16,cex=0.7,col=adjustcolor("orange",alpha.f = 0.8))

points(g2[,2],g2[,4],pch=16,cex=1,col="red")
text(g2[,2],g2[,4],"ZmR1",cex=0.8,pos=4,font=3)

####R1 targets OE3
plot(re2[,2],re2[,5],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="MPH_J92/JING724",ylab="MPH_J92/OE3",
     xlim=c(-100,300),ylim=c(-100,300),las=1,tck=-0.03)

cor.test(re2[,2],re2[,5])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))
points(re3[,2],re3[,5],pch=16,cex=0.7,col=adjustcolor("darkred",alpha.f = 0.5))
points(re4[,2],re4[,5],pch=16,cex=0.7,col=adjustcolor("purple",alpha.f = 0.4))
points(re5[,2],re5[,5],pch=16,cex=0.7,col=adjustcolor("darkgreen",alpha.f = 0.4))
points(re6[,2],re6[,5],pch=16,cex=0.7,col=adjustcolor("orange",alpha.f = 0.8))

points(g2[,2],g2[,5],pch=16,cex=1,col="red")
text(g2[,2],g2[,5],"ZmR1",cex=0.8,pos=4,font=3)

dev.off()

####R1 co-expressed genes
pdf("OE_hyb_MPH_coexp_inb_hyb.pdf",width=8.6,height=5)
par(mfrow=c(2,3),mar=c(4,4,2,2))
##co-expressed in inbreds
####R1 targets OE1
plot(re2[,2],re2[,3],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="MPH_J92/JING724",ylab="MPH_J92/OE1",
     xlim=c(-100,300),ylim=c(-100,300),las=1)

cor.test(re2[,2],re2[,3])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))

re3=re2[re2[,1]%in%d1,]
po=re3[re3[,1]%in%pos[,1],]
ne=re3[re3[,1]%in%neg[,1],]
points(po[,2],po[,3],pch=16,cex=0.8,col=adjustcolor("darkred",alpha.f = 0.4))
points(ne[,2],ne[,3],pch=16,cex=0.8,col=adjustcolor("darkgreen",alpha.f = 0.4))
x=po[,3]-po[,2]
n=length(x[x>0])/length(x)

####R1 targets OE2
plot(re2[,2],re2[,4],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="MPH_J92/JING724",ylab="MPH_J92/OE2",
     xlim=c(-100,300),ylim=c(-100,300),las=1)

cor.test(re2[,2],re2[,4])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))

points(po[,2],po[,4],pch=16,cex=0.8,col=adjustcolor("darkred",alpha.f = 0.4))
points(ne[,2],ne[,4],pch=16,cex=0.8,col=adjustcolor("darkgreen",alpha.f = 0.4))
x=po[,4]-po[,2]
n=length(x[x>0])/length(x)
n

#legend("topright",legend=c("R1_tar_pos","R1_tar_neg","Anthocyanin"),col=adjustcolor(c("darkred","darkgreen","orange"),0.6),pch=16,cex=0.7,bty="n")

####R1 targets OE3
plot(re2[,2],re2[,5],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="MPH_J92/JING724",ylab="MPH_J92/OE3",
     xlim=c(-100,300),ylim=c(-100,300),las=1)

cor.test(re2[,2],re2[,5])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))

points(po[,2],po[,5],pch=16,cex=0.8,col=adjustcolor("darkred",alpha.f = 0.4))
points(ne[,2],ne[,5],pch=16,cex=0.8,col=adjustcolor("darkgreen",alpha.f = 0.4))
x=po[,5]-po[,2]
n=length(x[x>0])/length(x)
n

##co-expressed in hybrids
####R1 targets OE1
plot(re2[,2],re2[,3],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="MPH_J92/JING724",ylab="MPH_J92/OE1",
     xlim=c(-100,300),ylim=c(-100,300),las=1)

cor.test(re2[,2],re2[,3])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))

re3=re2[re2[,1]%in%d2,]
po=re3[re3[,1]%in%pos[,1],]
ne=re3[re3[,1]%in%neg[,1],]
points(po[,2],po[,3],pch=16,cex=0.8,col=adjustcolor("darkred",alpha.f = 0.4))
points(ne[,2],ne[,3],pch=16,cex=0.8,col=adjustcolor("darkgreen",alpha.f = 0.4))
x=po[,3]-po[,2]
n=length(x[x>0])/length(x)


####R1 targets OE2
plot(re2[,2],re2[,4],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="MPH_J92/JING724",ylab="MPH_J92/OE2",
     xlim=c(-100,300),ylim=c(-100,300),las=1)

cor.test(re2[,2],re2[,4])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))

points(po[,2],po[,4],pch=16,cex=0.8,col=adjustcolor("darkred",alpha.f = 0.4))
points(ne[,2],ne[,4],pch=16,cex=0.8,col=adjustcolor("darkgreen",alpha.f = 0.4))
x=po[,4]-po[,2]
n=length(x[x>0])/length(x)


####R1 targets OE3
plot(re2[,2],re2[,5],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="MPH_J92/JING724",ylab="MPH_J92/OE3",
     xlim=c(-100,300),ylim=c(-100,300),las=1)

cor.test(re2[,2],re2[,5])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))

points(po[,2],po[,5],pch=16,cex=0.8,col=adjustcolor("darkred",alpha.f = 0.4))
points(ne[,2],ne[,5],pch=16,cex=0.8,col=adjustcolor("darkgreen",alpha.f = 0.4))
x=po[,5]-po[,2]
n=length(x[x>0])/length(x)

legend("topright",legend=c("Positive","Negative"),col=adjustcolor(c("darkred","darkgreen"),0.6),pch=16,cex=0.7,bty="n")
dev.off()





pdf("R1_nj_MPH_V4.pdf",width=8.6,height=4)
par(mfrow=c(1,2),mar=c(4,4,2,2))
####R1 targets R1nj
plot(re2[,8],re2[,9],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="C7-2_Qi319",ylab="C7-2_Qi319nj",
     xlim=c(-100,300),ylim=c(-100,300),las=1)

cor.test(re2[,8],re2[,9])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))

re3=re2[re2[,1]%in%d[,1],]
re4=merge(re3,d,by=1)
po=re4[re4[,13]=="Positive",]
ne=re4[re4[,13]=="Negetive",]
points(po[,8],po[,9],pch=16,cex=0.7,col=adjustcolor("darkred",alpha.f = 0.4))
points(ne[,8],ne[,9],pch=16,cex=0.7,col=adjustcolor("darkgreen",alpha.f = 0.4))

points(g2[,9],g2[,10],pch=16,cex=0.7,col=g2$col)
text(g2[,9],g2[,10],g2[,2],cex=0.8,pos=4)

#####plot MB
plot(re2[,10],re2[,11],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="MO17XB73",ylab="MO17XB73nj",
     xlim=c(-100,300),ylim=c(-100,300),las=1)

cor.test(re2[,10],re2[,11])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))

re3=re2[re2[,1]%in%d[,1],]
re4=merge(re3,d,by=1)
po=re4[re4[,13]=="Positive",]
ne=re4[re4[,13]=="Negetive",]
points(po[,10],po[,11],pch=16,cex=0.7,col=adjustcolor("darkred",alpha.f = 0.4))
points(ne[,10],ne[,11],pch=16,cex=0.7,col=adjustcolor("darkgreen",alpha.f = 0.4))

points(g2[,11],g2[,12],pch=16,cex=0.7,col=g2$col)
text(g2[,11],g2[,12],g2[,2],cex=0.8,pos=4)
legend("topleft",legend=c("R1_tar_pos","R1_tar_neg","Anthocyanin"),col=adjustcolor(c("darkred","darkgreen","orange"),0.6),pch=16,cex=0.7,bty="n")
dev.off()



pdf("EMS_hyb_MPH_co_exp_inb_v2.pdf",width=5,height=5)
par(mfrow=c(1,1),mar=c(4,4,2,2))
####R1 targets EMS
plot(re2[,6],re2[,7],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="MPH_EMS_AA",ylab="EMS_aa",
     xlim=c(-100,300),ylim=c(-100,300),las=1)

cor.test(re2[,6],re2[,7])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))

re3=re2[re2[,1]%in%d1,]
po=re3[re3[,1]%in%pos[,1],]
ne=re3[re3[,1]%in%neg[,1],]
points(po[,6],po[,7],pch=16,cex=0.7,col=adjustcolor("darkred",alpha.f = 0.4))
points(ne[,6],ne[,7],pch=16,cex=0.7,col=adjustcolor("darkgreen",alpha.f = 0.4))
x=po[,7]-po[,6]
n=length(x[x>0])/length(x)

g2=merge(g,re2,by.x="Gene")
g2$col=adjustcolor("orange",alpha.f = 0.6)
g2[g2[,1]=="Zm00001d026147",13]="blue"
points(g2[,7],g2[,8],pch=16,cex=0.7,col=g2$col)
text(g2[,7],g2[,8],g2[,2],cex=0.8,pos=4)
legend("topleft",legend=c("R1_tar_pos","R1_tar_neg","Anthocyanin"),col=adjustcolor(c("darkred","darkgreen","orange"),0.6),pch=16,cex=0.7,bty="n")
dev.off()


pdf("R1_nj_MPH_co_exp_inb_v2.pdf",width=8.6,height=4)
par(mfrow=c(1,2),mar=c(4,4,2,2))
####R1 targets R1nj
plot(re2[,8],re2[,9],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="C7-2_Qi319",ylab="C7-2_Qi319nj",
     xlim=c(-100,300),ylim=c(-100,300),las=1)

cor.test(re2[,8],re2[,9])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))

re3=re2[re2[,1]%in%d1,]
po=re3[re3[,1]%in%pos[,1],]
ne=re3[re3[,1]%in%neg[,1],]
points(po[,8],po[,9],pch=16,cex=0.7,col=adjustcolor("darkred",alpha.f = 0.4))
points(ne[,8],ne[,9],pch=16,cex=0.7,col=adjustcolor("darkgreen",alpha.f = 0.4))
x=po[,9]-po[,8]
n=length(x[x>0])/length(x)


points(g2[,9],g2[,10],pch=16,cex=0.7,col=g2$col)
text(g2[,9],g2[,10],g2[,2],cex=0.8,pos=4)

#####plot MB
plot(re2[,10],re2[,11],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="MO17XB73",ylab="MO17XB73nj",
     xlim=c(-100,300),ylim=c(-100,300),las=1)

cor.test(re2[,10],re2[,11])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))

re3=re2[re2[,1]%in%d1,]
po=re3[re3[,1]%in%pos[,1],]
ne=re3[re3[,1]%in%neg[,1],]
points(po[,10],po[,11],pch=16,cex=0.7,col=adjustcolor("darkred",alpha.f = 0.4))
points(ne[,10],ne[,11],pch=16,cex=0.7,col=adjustcolor("darkgreen",alpha.f = 0.4))
x=po[,11]-po[,10]
n=length(x[x>0])/length(x)

points(g2[,11],g2[,12],pch=16,cex=0.7,col=g2$col)
text(g2[,11],g2[,12],g2[,2],cex=0.8,pos=4)
legend("topleft",legend=c("R1_tar_pos","R1_tar_neg","Anthocyanin"),col=adjustcolor(c("darkred","darkgreen","orange"),0.6),pch=16,cex=0.7,bty="n")
dev.off()

##############################################
#####plot R1 co-expressed genes in hybrid
pdf("OE_hyb_MPH_coexp_hyb_v2.pdf",width=8.6,height=3)
par(mfrow=c(1,3),mar=c(4,4,2,2))
####R1 targets OE1
plot(re2[,2],re2[,3],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="MPH_J92/JING724",ylab="MPH_J92/OE1",
     xlim=c(-100,300),ylim=c(-100,300),las=1)

cor.test(re2[,2],re2[,3])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=2, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线

re3=re2[re2[,1]%in%d2,]
po=re3[re3[,1]%in%pos[,1],]
ne=re3[re3[,1]%in%neg[,1],]
points(po[,2],po[,3],pch=16,cex=0.7,col=adjustcolor("darkred",alpha.f = 0.4))
points(ne[,2],ne[,3],pch=16,cex=0.7,col=adjustcolor("darkgreen",alpha.f = 0.4))
x=po[,3]-po[,2]
n=length(x[x>0])/length(x)

g2=merge(g,re2,by.x="Gene")
g2$col=adjustcolor("orange",alpha.f = 0.6)
g2[g2[,1]=="Zm00001d026147",13]="blue"
points(g2[,3],g2[,4],pch=16,cex=0.7,col=g2$col)
text(g2[,3],g2[,4],g2[,2],cex=0.8,pos=4)

####R1 targets OE2
plot(re2[,2],re2[,4],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="MPH_J92/JING724",ylab="MPH_J92/OE2",
     xlim=c(-100,300),ylim=c(-100,300),las=1)

cor.test(re2[,2],re2[,4])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))

points(po[,2],po[,4],pch=16,cex=0.7,col=adjustcolor("darkred",alpha.f = 0.4))
points(ne[,2],ne[,4],pch=16,cex=0.7,col=adjustcolor("darkgreen",alpha.f = 0.4))
x=po[,4]-po[,2]
n=length(x[x>0])/length(x)

points(g2[,3],g2[,5],pch=16,cex=0.7,col=g2$col)
text(g2[,3],g2[,5],g2[,2],cex=0.8,pos=4)
#legend("topright",legend=c("R1_tar_pos","R1_tar_neg","Anthocyanin"),col=adjustcolor(c("darkred","darkgreen","orange"),0.6),pch=16,cex=0.7,bty="n")

####R1 targets OE3
plot(re2[,2],re2[,5],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="MPH_J92/JING724",ylab="MPH_J92/OE3",
     xlim=c(-100,300),ylim=c(-100,300),las=1)

cor.test(re2[,2],re2[,5])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))

points(po[,2],po[,5],pch=16,cex=0.7,col=adjustcolor("darkred",alpha.f = 0.4))
points(ne[,2],ne[,5],pch=16,cex=0.7,col=adjustcolor("darkgreen",alpha.f = 0.4))
x=po[,5]-po[,2]
n=length(x[x>0])/length(x)

points(g2[,3],g2[,6],pch=16,cex=0.7,col=g2$col)
text(g2[,3],g2[,6],g2[,2],cex=0.8,pos=4)
legend("topleft",legend=c("R1_Coexp_hyb_pos","R1_Coexp_hyb_neg","Anthocyanin"),col=adjustcolor(c("darkred","darkgreen","orange"),0.6),pch=16,cex=0.7,bty="n")
dev.off()

pdf("EMS_hyb_MPH_co_exp_hyb_v2.pdf",width=5,height=4)
par(mfrow=c(1,1),mar=c(4,4,2,2))
####R1 targets EMS
plot(re2[,6],re2[,7],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="MPH_EMS_AA",ylab="EMS_aa",
     xlim=c(-100,300),ylim=c(-100,300),las=1)

cor.test(re2[,6],re2[,7])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))

re3=re2[re2[,1]%in%d2,]
po=re3[re3[,1]%in%pos[,1],]
ne=re3[re3[,1]%in%neg[,1],]
points(po[,6],po[,7],pch=16,cex=0.7,col=adjustcolor("darkred",alpha.f = 0.4))
points(ne[,6],ne[,7],pch=16,cex=0.7,col=adjustcolor("darkgreen",alpha.f = 0.4))
x=po[,7]-po[,6]
n=length(x[x>0])/length(x)

g2=merge(g,re2,by.x="Gene")
g2$col=adjustcolor("orange",alpha.f = 0.6)
g2[g2[,1]=="Zm00001d026147",13]="blue"
points(g2[,7],g2[,8],pch=16,cex=0.7,col=g2$col)
text(g2[,7],g2[,8],g2[,2],cex=0.8,pos=4)
legend("topleft",legend=c("R1_Coexp_hyb_pos","R1_Coexp_hyb_neg","Anthocyanin"),col=adjustcolor(c("darkred","darkgreen","orange"),0.6),pch=16,cex=0.7,bty="n")
dev.off()


pdf("R1_nj_MPH_co_exp_hyb_v2.pdf",width=8.6,height=4)
par(mfrow=c(1,2),mar=c(4,4,2,2))
####R1 targets R1nj
plot(re2[,8],re2[,9],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="C7-2_Qi319",ylab="C7-2_Qi319nj",
     xlim=c(-100,300),ylim=c(-100,300),las=1)

cor.test(re2[,8],re2[,9])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))

re3=re2[re2[,1]%in%d2,]
po=re3[re3[,1]%in%pos[,1],]
ne=re3[re3[,1]%in%neg[,1],]
points(po[,8],po[,9],pch=16,cex=0.7,col=adjustcolor("darkred",alpha.f = 0.4))
points(ne[,8],ne[,9],pch=16,cex=0.7,col=adjustcolor("darkgreen",alpha.f = 0.4))
x=po[,9]-po[,8]
n=length(x[x>0])/length(x)

points(g2[,9],g2[,10],pch=16,cex=0.7,col=g2$col)
text(g2[,9],g2[,10],g2[,2],cex=0.8,pos=4)

#####plot MB
plot(re2[,10],re2[,11],pch=16,cex=0.5,col=adjustcolor("gray",alpha.f = 0.4),xlab="MO17XB73",ylab="MO17XB73nj",
     xlim=c(-100,300),ylim=c(-100,300),las=1)

cor.test(re2[,10],re2[,11])
abline(a=0, b=1, col=adjustcolor("blue",alpha.f = 0.6), lwd=1, lty=2)  # 绘制从 (0,0) 到 (10,10) 的对角线
abline(v=0,h=0,lty=2,col=adjustcolor("blue",alpha.f = 0.6))

re3=re2[re2[,1]%in%d2,]
po=re3[re3[,1]%in%pos[,1],]
ne=re3[re3[,1]%in%neg[,1],]
points(po[,10],po[,11],pch=16,cex=0.7,col=adjustcolor("darkred",alpha.f = 0.4))
points(ne[,10],ne[,11],pch=16,cex=0.7,col=adjustcolor("darkgreen",alpha.f = 0.4))
x=po[,11]-po[,10]
n=length(x[x>0])/length(x)

points(g2[,11],g2[,12],pch=16,cex=0.7,col=g2$col)
text(g2[,11],g2[,12],g2[,2],cex=0.8,pos=4)
legend("topleft",legend=c("R1_Coexp_hyb_pos","R1_Coexp_hyb_neg","Anthocyanin"),col=adjustcolor(c("darkred","darkgreen","orange"),0.6),pch=16,cex=0.7,bty="n")
dev.off()
