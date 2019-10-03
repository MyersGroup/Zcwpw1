setwd("C://Users//myers.WHGC//Dropbox//ZCWPW1")

this=read.table("ssDNA_ZCWPW1_HOM_290519_rep1_type1_filtered_only_rmdup.chr1.3prime.bedgraph",as.is=T)
that=read.table("ssDNA_ZCWPW1_HOM_290519_rep1_type1_filtered_only_rmdup.chr1.5prime.bedgraph",as.is=T)

this=this[this[,1] %in% paste("chr",c(1:19,"X","Y"),sep=""),]
that=that[that[,1] %in% paste("chr",c(1:19,"X","Y"),sep=""),]

##this=this[this[,3]-this[,2]>40 & this[,3]-this[,2]<=100,]
##that=that[that[,3]-that[,2]>40 & that[,3]-that[,2]<=100,]


b6=read.table("b6.txt",as.is=T,header=T)

####overlap

startends=seq(-5000,5000,20)
totals=matrix(0,nrow=21,ncol=length(startends)-1)
curtotals=totals[1,]
hh=1
for(chrom in c(1:19,"X"))
{
print("New chrom: ")
print(chrom)
curtotals=curtotals*0
poshot=as.double(b6[b6[,1]==hh,2])
ourreads=this[this[,1]== paste("chr",chrom,sep=""),2:3]

for(j in 1:length(poshot)){
			cond=ourreads[,2]>=poshot[j]-5000 & ourreads[,1]<=poshot[j]+5000
			if(sum(cond)){
			temp=as.matrix(ourreads[cond,],ncol=2)
			print(c(j,length(temp)))

			temp=temp-poshot[j]+5000
			indicesl=floor(temp[,1]/20)+1
			indicesr=floor(temp[,2]/20)+1
			indicesl[indicesl<=0]=1
			
			indicesr[(indicesr>500)]=500
			for(k in 1:nrow(temp)) curtotals[indicesl[k]:indicesr[k]]=curtotals[indicesl[k]:indicesr[k]]+1	
			}	
		}
totals[hh,]=curtotals
hh=hh+1
}
totals3prime=totals
startends=seq(-5000,5000,20)
totals=matrix(0,nrow=21,ncol=length(startends)-1)
curtotals=totals[1,]
hh=1
for(chrom in c(1:19,"X"))
{
print("New chrom: ")
print(chrom)
curtotals=curtotals*0
poshot=as.double(b6[b6[,1]==hh,2])
ourreads=that[that[,1]== paste("chr",chrom,sep=""),2:3]

for(j in 1:length(poshot)){
			cond=ourreads[,2]>=poshot[j]-5000 & ourreads[,1]<=poshot[j]+5000
			if(sum(cond)){
			temp=as.matrix(ourreads[cond,],ncol=2)
			print(c(j,length(temp)))

			temp=temp-poshot[j]+5000
			indicesl=floor(temp[,1]/20)+1
			indicesr=floor(temp[,2]/20)+1
			indicesl[indicesl<=0]=1
			
			indicesr[(indicesr>500)]=500
			for(k in 1:nrow(temp)) curtotals[indicesl[k]:indicesr[k]]=curtotals[indicesl[k]:indicesr[k]]+1	
			}	
		}
totals[hh,]=curtotals
hh=hh+1
}

totals5prime=totals
plot(colSums(totals3prime),type="l",col=2)
lines(colSums(totals5prime),col=4)
abline(v=250,lty="dotted")
dev.copy2pdf(file="HOM_type1_strand_separated.pdf",width=6,height=4)


####HET

this=read.table("ssDNA_ZCWPW1_HET_290519_rep1_type1_filtered_only_rmdup.chr1.3prime.bedgraph",as.is=T)
that=read.table("ssDNA_ZCWPW1_HET_290519_rep1_type1_filtered_only_rmdup.chr1.5prime.bedgraph",as.is=T)

this=this[this[,1] %in% paste("chr",c(1:19,"X","Y"),sep=""),]
that=that[that[,1] %in% paste("chr",c(1:19,"X","Y"),sep=""),]
this=this[this[,3]-this[,2]<=100,]
that=that[that[,3]-that[,2]<=100,]

b6=read.table("b6.txt",as.is=T,header=T)

####overlap

startends=seq(-5000,5000,20)
totals=matrix(0,nrow=21,ncol=length(startends)-1)
curtotals=totals[1,]
hh=1
for(chrom in c(1:19,"X"))
{
print("New chrom: ")
print(chrom)
curtotals=curtotals*0
poshot=as.double(b6[b6[,1]==hh,2])
ourreads=this[this[,1]== paste("chr",chrom,sep=""),2:3]

for(j in 1:length(poshot)){
			cond=ourreads[,2]>=poshot[j]-5000 & ourreads[,1]<=poshot[j]+5000
			if(sum(cond)){
			temp=as.matrix(ourreads[cond,],ncol=2)
			print(c(j,length(temp)))

			temp=temp-poshot[j]+5000
			indicesl=floor(temp[,1]/20)+1
			indicesr=floor(temp[,2]/20)+1
			indicesl[indicesl<=0]=1
			
			indicesr[(indicesr>500)]=500
			for(k in 1:nrow(temp)) curtotals[indicesl[k]:indicesr[k]]=curtotals[indicesl[k]:indicesr[k]]+1	
			}	
		}
totals[hh,]=curtotals
hh=hh+1
}
totals3prime=totals
startends=seq(-5000,5000,20)
totals=matrix(0,nrow=21,ncol=length(startends)-1)
curtotals=totals[1,]
hh=1
for(chrom in c(1:19,"X"))
{
print("New chrom: ")
print(chrom)
curtotals=curtotals*0
poshot=as.double(b6[b6[,1]==hh,2])
ourreads=that[that[,1]== paste("chr",chrom,sep=""),2:3]

for(j in 1:length(poshot)){
			cond=ourreads[,2]>=poshot[j]-5000 & ourreads[,1]<=poshot[j]+5000
			if(sum(cond)){
			temp=as.matrix(ourreads[cond,],ncol=2)
			print(c(j,length(temp)))

			temp=temp-poshot[j]+5000
			indicesl=floor(temp[,1]/20)+1
			indicesr=floor(temp[,2]/20)+1
			indicesl[indicesl<=0]=1
			
			indicesr[(indicesr>500)]=500
			for(k in 1:nrow(temp)) curtotals[indicesl[k]:indicesr[k]]=curtotals[indicesl[k]:indicesr[k]]+1	
			}	
		}
totals[hh,]=curtotals
hh=hh+1
}

totals5prime=totals
plot(colSums(totals3prime),type="l",col=2)
lines(colSums(totals5prime),col=4)
abline(v=250,lty="dotted")
dev.copy2pdf(file="HET_type1_strand_separated.pdf",width=6,height=4)

######NO SIGNAL AT ALL

totchrom=matrix(nrow=21,ncol=3)

for(i in 1:length(totchrom)) totchrom[i,1]=sum(b6[b6[,1]==i,4])
for(i in 1:length(totchrom)) totchrom[i,2]=sum(totals3prime[i,180:320])
for(i in 1:length(totchrom)) totchrom[i,3]=sum(totals5prime[i,180:320])



#########now look at reads per hotspot

counts=1:nrow(b6)*0

curtotals=totals[1,]
hh=1
for(chrom in c(1:19,"X"))
{
print("New chrom: ")
print(chrom)
curtotals=curtotals*0
poshot=as.double(b6[b6[,1]==hh,2])
ourreads=this[this[,1]== paste("chr",chrom,sep=""),2:3]
tempcounts=1:length(poshot)*0
for(j in 1:length(poshot)){
			cond=ourreads[,2]>=poshot[j]-1250 & ourreads[,1]<=poshot[j]+1250
			tempcounts[j]=sum(cond)
					
}

counts[b6[,1]==hh]=tempcounts
hh=hh+1
}
counts3prime=counts


counts=1:nrow(b6)*0

curtotals=totals[1,]
hh=1
for(chrom in c(1:19,"X"))
{
print("New chrom: ")
print(chrom)
curtotals=curtotals*0
poshot=as.double(b6[b6[,1]==hh,2])
ourreads=that[that[,1]== paste("chr",chrom,sep=""),2:3]
tempcounts=1:length(poshot)*0
for(j in 1:length(poshot)){
			cond=ourreads[,2]>=poshot[j]-1250 & ourreads[,1]<=poshot[j]+1250
			tempcounts[j]=sum(cond)
					
}

counts[b6[,1]==hh]=tempcounts
hh=hh+1
}
counts5prime=counts

datacheck=cbind(b6[,3],counts3prime,counts5prime)
######reads vs. heat

startends=quantile(datacheck[,1],seq(0,1,0.02))
startends2=quantile(datacheck[b6[,1]==20,1],seq(0,1,.05))

mm=cbind(1:(length(startends)-1),1:(length(startends)-1))
mm=cbind(mm,mm)
for(i in 1:(length(startends)-1)){
cond=datacheck[,1]>=startends[i] & datacheck[,1]<startends[i+1] & b6[,1]!=20
if(i<=20){
 cond2=datacheck[,1]>=startends2[i] & datacheck[,1]<startends2[i+1] & b6[,1]==20
mm[i,3]=mean(datacheck[cond2,1])

mm[i,4]=mean(datacheck[cond2,2]+datacheck[cond2,3])


}

mm[i,1]=mean(datacheck[cond,1])
mm[i,2]=mean(datacheck[cond,2]+datacheck[cond,3])
}


plot(mm[,1],mm[,2],pch=19,cex=1,xlab="Heat in WT B6 (DMC1)",ylab="Heat in ZCWPW1 KO (DMC1)")
points(mm[1:20,3],mm[1:20,4],col=2,pch=19,cex=0.5)
dev.copy2pdf(file="DMC1comparisonWTvsZCWPW1KO_X_chrom_red.pdf",width=8,height=8)

load("b6full.out")
datacheck=cbind(b6full[,"SPO11"],counts3prime,counts5prime)
######reads vs. heat

startends=quantile(datacheck[,1],seq(0,1,0.04))
startends2=quantile(datacheck[b6[,1]==20,1],seq(0,1,.1))

mm=matrix(nrow=length(startends)-1,ncol=7)
mm=cbind(mm,mm)
for(i in 1:(length(startends)-1)){
cond=datacheck[,1]>=startends[i] & datacheck[,1]<startends[i+1] & b6[,1]!=20 & b6full[,"hshared"]==0
if(i<=20){
 cond2=datacheck[,1]>=startends2[i] & datacheck[,1]<startends2[i+1] & b6[,1]==20 & b6full[,"hshared"]==0
cond3=(b6full[,"hshared"]==0)


mm[i,1+7]=mean(b6full[cond2,"dmc1_orig_heat"])/mean(b6full[cond3,"dmc1_orig_heat"])
mm[i,2+7]=mean(datacheck[cond2,2]+datacheck[cond2,3])/mean(datacheck[cond3,3]+datacheck[cond3,2])
mm[i,3+7]=mean(b6full[cond2,"SPO11"])/mean(b6full[cond3,"SPO11"])
mm[i,4+7]=mean(b6full[cond2 & b6full[,"hshared"]==0,"enrichment"])/mean(b6full[cond3,"enrichment"])
mm[i,5+7]=mean(b6full[cond2,"DMC1"])/mean(b6full[cond3,"DMC1"])
mm[i,6+7]=mean(b6full[cond2,"RPA"])/mean(b6full[cond3,"RPA"])
mm[i,7+7]=mean(b6full[cond2,"RAD51"])/mean(b6full[cond3,"RAD51"])


}

mm[i,1]=mean(b6full[cond,"dmc1_orig_heat"])/mean(b6full[cond3,"dmc1_orig_heat"])
mm[i,2]=mean(datacheck[cond,2]+datacheck[cond,3])/mean(datacheck[cond3,2]+datacheck[cond3,3])
mm[i,3]=mean(b6full[cond,"SPO11"])/mean(b6full[cond3,"SPO11"])
mm[i,4]=mean(b6full[cond ,"enrichment"])/mean(b6full[cond3,"enrichment"])
mm[i,5]=mean(b6full[cond,"DMC1"])/mean(b6full[cond3,"DMC1"])
mm[i,6]=mean(b6full[cond,"RPA"])/mean(b6full[cond3,"RPA"])
mm[i,7]=mean(b6full[cond,"RAD51"])/mean(b6full[cond3,"RAD51"])


}
ttt=c("DMC1_orig","DMC1_ZCW","SPO11","H3K4","DMC1","RPA","RAD51")
ttt=c(ttt,paste(ttt,"_chrX",sep=""))
colnames(mm)=ttt
refcol=3
##plot(mm[,refcol],mm[,2],pch=19,cex=1,xlab=paste("Heat in WT B6 (",colnames(mm)[refcol],")",sep=""),ylab="Heat in ZCWPW1 KO (DMC1)")
###points(mm[1:20,refcol+7],mm[1:20,2+7],col=2,pch=19,cex=0.5)

#####bootstraps

######reads vs. heat


mmbb=array(dim=c(length(startends)-1,14,1000))
##mm=cbind(mm,mm)
for(i in 1:(length(startends)-1)){
cond=datacheck[,1]>=startends[i] & datacheck[,1]<startends[i+1] & b6[,1]!=20 & b6full[,"hshared"]==0
cond3=(b6full[,"hshared"]==0)
print(i)
for(j in 1:1000){
condb=vvv=sample(which(cond==T),replace=T)

if(i<=10){
 cond2=datacheck[,1]>=startends2[i] & datacheck[,1]<startends2[i+1] & b6[,1]==20 & b6full[,"hshared"]==0
cond2b=sample(which(cond2==T),replace=T)

mmbb[i,1+7,j]=mean(b6full[cond2b,"dmc1_orig_heat"])/mean(b6full[cond3,"dmc1_orig_heat"])
mmbb[i,2+7,j]=mean(datacheck[cond2b,2]+datacheck[cond2b,3])/mean(datacheck[cond3,3]+datacheck[cond3,2])
mmbb[i,3+7,j]=mean(b6full[cond2b,"SPO11"])/mean(b6full[cond3,"SPO11"])
mmbb[i,4+7,j]=mean(b6full[cond2b ,"enrichment"])/mean(b6full[cond3,"enrichment"])
mmbb[i,5+7,j]=mean(b6full[cond2b,"DMC1"])/mean(b6full[cond3,"DMC1"])
mmbb[i,6+7,j]=mean(b6full[cond2b,"RPA"])/mean(b6full[cond3,"RPA"])
mmbb[i,7+7,j]=mean(b6full[cond2b,"RAD51"])/mean(b6full[cond3,"RAD51"])


}
##}
##}

mmbb[i,1,j]=mean(b6full[condb,"dmc1_orig_heat"])/mean(b6full[cond3,"dmc1_orig_heat"])
mmbb[i,2,j]=mean(datacheck[condb,2]+datacheck[condb,3])/mean(datacheck[cond3,2]+datacheck[cond3,3])
mmbb[i,3,j]=mean(b6full[condb,"SPO11"])/mean(b6full[cond3,"SPO11"])
mmbb[i,4,j]=mean(b6full[condb ,"enrichment"])/mean(b6full[cond3,"enrichment"])
mmbb[i,5,j]=mean(b6full[condb,"DMC1"])/mean(b6full[cond3,"DMC1"])
mmbb[i,6,j]=mean(b6full[condb,"RPA"])/mean(b6full[cond3,"RPA"])
mmbb[i,7,j]=mean(b6full[condb,"RAD51"])/mean(b6full[cond3,"RAD51"])


}

}

mmqql=mm*0
mmqqu=mm*0
for(i in 1:nrow(mmqql)) for(j in 1:ncol(mmqql)) mmqql[i,j]=quantile(mmbb[i,j,],0.025,na.rm=T)
for(i in 1:nrow(mmqql)) for(j in 1:ncol(mmqql)) mmqqu[i,j]=quantile(mmbb[i,j,],0.975,na.rm=T)

par(mfrow=c(1,2))
refcol=3
plot(mm[,refcol],mm[,"DMC1"],pch=19,cex=1,xlab=paste("Heat in WT B6 (",colnames(mm)[refcol],")",sep=""),ylab="Heat in WT (DMC1)",type="n")
segments(x0=mm[,refcol+7],x1=mm[,refcol+7],y0=mmqql[1:20,5+7],y1=mmqqu[1:20,5+7],lty="dotted")
segments(x0=mm[,refcol],x1=mm[,refcol],y0=mmqql[,5],y1=mmqqu[,5])
points(mm[,refcol],mm[,"DMC1"],pch=19,cex=1,xlab=paste("Heat in WT B6 (",colnames(mm)[refcol],")",sep=""),ylab="Heat in WT (DMC1)")
points(mm[1:20,3+7],mm[1:20,5+7],col=2,pch=19,cex=0.5)
ww=seq(0,10,0.001)

plot(mm[,refcol],mm[,2],pch=19,cex=1,xlab=paste("Heat in WT B6 (",colnames(mm)[refcol],")",sep=""),ylab="Heat in ZCWPW1 KO (DMC1)",type="n")
segments(x0=mm[,refcol+7],x1=mm[,refcol+7],y0=mmqql[1:20,2+7],y1=mmqqu[1:20,2+7],lty="dotted")
segments(x0=mm[,refcol],x1=mm[,refcol],y0=mmqql[,2],y1=mmqqu[,2])
points(mm[,refcol],mm[,2],pch=19,cex=1,xlab=paste("Heat in WT B6 (",colnames(mm)[refcol],")",sep=""),ylab="Heat in ZCWPW1 KO (DMC1)")
points(mm[1:20,refcol+7],mm[1:20,2+7],col=2,pch=19,cex=0.5)


dev.copy2pdf(file="SPO11comparisonWTvsZCWPW1KO_X_chrom_red.pdf",width=8,height=4)

####red is chromosome X; plot is based on 1000 bootstrapped hotspots from B6. DMC1-mapped hotspots are binned by their SPO11 activity and then we compare estimated DMC1 heats.

####Note that the left plot is clearly curved, suggesting slower DSB repair for weaker hotspots. The X-chromosome is hugely elevated  suggesting persistent unrepaired DSBs there, as can be observed cytogenetically.

#####On the right a straight line apart from theg final X-chromosome point implying DMC1 removal no longer ?occurs? before arrest. And coupling between PRDM9 binding and DMC1 persistence now disappears. 

#####any other clusters?

temp=rbind(this,that)
temp=temp[order(temp[,1],temp[,2]),]
 
####clusters

clustersize=5000
clusters=1:nrow(temp)
temp2=temp[temp[,1]=="chr1",]
for(i in 1:nrow(temp)){
if(temp[i,1]!=temp2[1,1]) temp2=temp[temp[,1]==temp[i,1],]

if(!i%%100) print(i)
clusters[i]=sum(temp2[,2]<=temp[i,2]+clustersize & temp2[,2]>=temp[i,2])
}

#####define clusters/overlaps?

overlap=1:nrow(temp)*0

width=1250
hh=1
for(chrom in c(1:19,"X")){
print(chrom)
ourset=b6[b6[,1]==hh,2]

tttt=temp[temp[,1]==paste("chr",chrom,sep=""),2:3]
ww=which(temp[,1]==paste("chr",chrom,sep="") )
fill=1:length(tttt)*0
for(k in 1:length(ourset)){
	if(!k%%100) print(k)
	fill[tttt[,1]<=ourset[k]+width & tttt[,2]>=ourset[k]-width]=1
}
overlap[ww[fill==1]]=1

hh=hh+1
}


#######newdata plots

setwd("C:/Users/myers.WHGC/Dropbox/ZCWPW1")
zcwfc=read.table("ForceCallFinal.DMC1.from.ZCWPW1_KO.into.B6_composite_500_2000.txt",header=F)
zcwfc[1:5,]
b6=read.table("B6_composite.txt",header=T)

#####modestly hot hotspots, and autosomal, controlled by B6 with at least some evidence in the KO

colnames(zcwfc)=c("Chr","Pos","Heat","Heat_2","Pval")


######find bins

zcwfc[,"Heat"]=zcwfc[,"Heat"]/mean(zcwfc[zcwfc[,1]<=19,"Heat"])
b6[,"DMC1"]=b6[,"DMC1"]/mean(b6[zcwfc[,1]<=19,"DMC1"])
b6[,"SPO11"]=b6[,"SPO11"]/mean(b6[zcwfc[,1]<=19,"SPO11"])
b6[,"enrichment"]=b6[,"enrichment"]/mean(b6[zcwfc[,1]<=19,"enrichment"])

cond=zcwfc[,1]<=19&b6[,"hshared"]==0 & b6[,"allele"]=="B6" & zcwfc[,"Heat"]>0 & b6[,"SPO11"]>0 & b6[,"DMC1"]>0 & b6[,"enrichment"]>0.5


www=b6[,"DMC1"]/b6[,"SPO11"] 
plot(b6[cond,c("enrichment")],log2(www[cond]),xlab="Relative H3K4me3",ylab="DMC1 to SPO11 log-ratio (WT)",pch=19,cex=0.1,col="grey",ylim=c(-3,4),xlim=c(0,6))


www2=zcwfc[,"Heat"]/b6[,"SPO11"] 
plot(b6[cond,c("enrichment")],log2(www2[cond]),xlab="Relative H3K4me3",ylab="DMC1 to SPO11 log-ratio (WT)",pch=19,cex=0.1,col="grey",ylim=c(-3,4),xlim=c(0,6))

#####quantiles

bins=quantile(b6[cond,"enrichment"],seq(0,1,0.05))
vals=matrix(nrow=length(bins)-1,ncol=4)
vals2=vals
vals3=vals
for(i in 1:nrow(vals)) {

ourrows=which(cond & b6[,"enrichment"]>=bins[i] & b6[,"enrichment"]<bins[i+1]) 

vals[i,2]=mean(log2(www[ourrows]))
vals2[i,2]=mean(log2(www2[ourrows]))
vals3[i,2]=mean(log2(www[ourrows]/www2[ourrows]))

##vals[i,2]=log2(mean(www[ourrows]))
##vals2[i,2]=log2(mean(www2[ourrows]))
vals[i,1]=(mean(b6[ourrows,"enrichment"]))
vals2[i,1]=(mean(b6[ourrows,"enrichment"]))
vv=sd(log2(www[ourrows]))/sqrt(length(ourrows))
vv2=sd(log2(www2[ourrows]))/sqrt(length(ourrows))
vv3=sd(log2(www[ourrows]/www2[ourrows]))/sqrt(length(ourrows))

vals[i,3]=vals[i,2]-2*vv
vals[i,4]=vals[i,2]+2*vv
vals2[i,3]=vals2[i,2]-2*vv2
vals2[i,4]=vals2[i,2]+2*vv2
vals3[i,3]=vals3[i,2]-2*vv3
vals3[i,4]=vals3[i,2]+2*vv3



}
vals3[,1]=vals2[,1]

##plot(vals[,1],vals[,2])
##plot(vals2[,1],vals2[,2])
par(mfrow=c(1,2))
plot(b6[cond,c("enrichment")],log2(www[cond]),xlab="Relative H3K4me3",ylab="DMC1 to SPO11 log2-ratio (WT)",pch=19,cex=0.1,col="grey",ylim=c(-3,4),xlim=c(0,6))
points(vals[,1],vals[,2],pch=19)
segments(x0=vals[,1],x1=vals[,1],y0=vals[,3],y1=vals[,4])



plot(b6[cond,c("enrichment")],log2(www2[cond]),xlab="Relative H3K4me3",ylab="DMC1 to SPO11 log2-ratio (ZCWPW1-null)",pch=19,cex=0.1,col="grey",ylim=c(-3,4),xlim=c(0,6))
points(vals2[,1],vals2[,2],pch=19)
segments(x0=vals2[,1],x1=vals2[,1],y0=vals2[,3],y1=vals2[,4])

dev.copy2pdf(file="log2ratiodmc1vsspo11wtandnull.pdf",width=10,height=5)

par(mfrow=c(1,1))
plot(b6[cond,c("enrichment")],log2(www[cond]/www2[cond]),xlab="Relative H3K4me3",ylab="DMC1 log2-ratio (WT vs ZCWPW1-null)",pch=19,cex=0.1,col="grey",ylim=c(-3,4),xlim=c(0,6))
points(vals3[,1],vals3[,2],pch=19,type="b")
segments(x0=vals3[,1],x1=vals3[,1],y0=vals3[,3],y1=vals3[,4])

dev.copy2pdf(file="log2ratiodmc1wtvsnull.pdf",width=5,height=5)

#####4.9-fold decrease in relative heat for strongest vs. weakest binding sites, for DMC1 vs. KO. (so in KO 4.9-fold hotter)
 
####2^(1.4088991--0.8867859) 4.909871




cond=zcwfc[,1]==20 &b6[,"hshared"]==0 & b6[,"allele"]=="B6" & zcwfc[,"Heat"]>0 & b6[,"SPO11"]>0 & b6[,"DMC1"]>0 & b6[,"enrichment"]>0.5


www=b6[,"DMC1"]/b6[,"SPO11"] 
plot(b6[cond,c("enrichment")],log2(www[cond]),xlab="Relative H3K4me3",ylab="DMC1 to SPO11 log2-ratio (WT)",pch=19,cex=0.1,col="grey",ylim=c(-3,4),xlim=c(0,6))


www2=zcwfc[,"Heat"]/b6[,"SPO11"] 
plot(b6[cond,c("enrichment")],log2(www2[cond]),xlab="Relative H3K4me3",ylab="DMC1 to SPO11 log2-ratio (WT)",pch=19,cex=0.1,col="grey",ylim=c(-3,4),xlim=c(0,6))

#####quantiles

bins=quantile(b6[cond,"enrichment"],seq(0,1,0.2))
vals=matrix(nrow=length(bins)-1,ncol=4)
vals2=vals
vals3=vals
for(i in 1:nrow(vals)) {

ourrows=which(cond & b6[,"enrichment"]>=bins[i] & b6[,"enrichment"]<bins[i+1]) 

vals[i,2]=mean(log2(www[ourrows]))
vals2[i,2]=mean(log2(www2[ourrows]))
vals3[i,2]=mean(log2(www[ourrows]/www2[ourrows]))

##vals[i,2]=log2(mean(www[ourrows]))
##vals2[i,2]=log2(mean(www2[ourrows]))
vals[i,1]=(mean(b6[ourrows,"enrichment"]))
vals2[i,1]=(mean(b6[ourrows,"enrichment"]))
vv=sd(log2(www[ourrows]))/sqrt(length(ourrows))
vv2=sd(log2(www2[ourrows]))/sqrt(length(ourrows))
vv3=sd(log2(www[ourrows]/www2[ourrows]))/sqrt(length(ourrows))

vals[i,3]=vals[i,2]-2*vv
vals[i,4]=vals[i,2]+2*vv
vals2[i,3]=vals2[i,2]-2*vv2
vals2[i,4]=vals2[i,2]+2*vv2
vals3[i,3]=vals3[i,2]-2*vv3
vals3[i,4]=vals3[i,2]+2*vv3



}
vals3[,1]=vals2[,1]

##plot(vals[,1],vals[,2])
##plot(vals2[,1],vals2[,2])
par(mfrow=c(1,2))
plot(b6[cond,c("enrichment")],log2(www[cond]),xlab="Relative H3K4me3",ylab="DMC1 to SPO11 log2-ratio (WT)",pch=19,cex=0.1,col="grey",ylim=c(-3,4),xlim=c(0,6))
points(vals[,1],vals[,2],pch=19)
segments(x0=vals[,1],x1=vals[,1],y0=vals[,3],y1=vals[,4])



plot(b6[cond,c("enrichment")],log2(www2[cond]),xlab="Relative H3K4me3",ylab="DMC1 to SPO11 log2-ratio (ZCWPW1-null)",pch=19,cex=0.1,col="grey",ylim=c(-3,4),xlim=c(0,6))
points(vals2[,1],vals2[,2],pch=19)
segments(x0=vals2[,1],x1=vals2[,1],y0=vals2[,3],y1=vals2[,4])

dev.copy2pdf(file="chrXlog2ratiodmc1vsspo11wtandnull.pdf",width=10,height=5)

par(mfrow=c(1,1))
plot(b6[cond,c("enrichment")],log2(www[cond]/www2[cond]),xlab="Relative H3K4me3",ylab="DMC1 log2-ratio (WT vs ZCWPW1-null)",pch=19,cex=0.1,col="grey",ylim=c(-3,4),xlim=c(0,6))
points(vals3[,1],vals3[,2],pch=19,type="b")
segments(x0=vals3[,1],x1=vals3[,1],y0=vals3[,3],y1=vals3[,4])

dev.copy2pdf(file="chrXlog2ratiodmc1wtvsnull.pdf",width=5,height=5)


########chrX and autosomes Spo11 heats
cond=zcwfc[,1]<=19&b6[,"hshared"]==0 & b6[,"allele"]=="B6" & zcwfc[,"Heat"]>0 & b6[,"SPO11"]>0 & b6[,"DMC1"]>0 & b6[,"enrichment"]>0.5
cond2=zcwfc[,1]==20&b6[,"hshared"]==0 & b6[,"allele"]=="B6" & zcwfc[,"Heat"]>0 & b6[,"SPO11"]>0 & b6[,"DMC1"]>0 & b6[,"enrichment"]>0.5


www=b6[,"DMC1"]
plot(b6[cond,"SPO11"],www[cond],xlab="Relative SPO11",ylab="DMC1 (WT)",pch=19,cex=0.1,col="grey")


www2=zcwfc[,"Heat"]
 
plot(b6[cond,"SPO11"],www2[cond],xlab="Relative SPO11",ylab="DMC1 (ZCWPW1-null)",pch=19,cex=0.1,col="grey")

www3=b6[,"SPO11"]
 
www4=b6[,"enrichment"]
www4=www3
###www4=www

#####quantiles

bins=quantile(www4[cond],seq(0,1,0.01))
bins2=quantile(www4[cond2],seq(0,1,0.1))

vals=matrix(nrow=length(bins)-1,ncol=4)
vals2=vals
vals3=vals

valsx=vals
vals2x=vals
vals3x=vals
for(i in 1:nrow(vals)) {

ourrows=which(cond & www4>=bins[i] & www4<bins[i+1]) 

if(i<length(bins2))
ourrows2=which(cond2 & www4>=bins2[i] & www4<bins2[i+1]) 

vals[i,2]=mean((www[ourrows]))
vals2[i,2]=mean((www2[ourrows]))
vals3[i,2]=mean((www[ourrows]/www2[ourrows]))

vals[i,1]=(mean(www3[ourrows]))
vals2[i,1]=vals[i,1]
vals3[i,1]=vals[i,1]
vv=sd((www[ourrows]))/sqrt(length(ourrows))
vv2=sd((www2[ourrows]))/sqrt(length(ourrows))
vv3=sd((www[ourrows]/www2[ourrows]))/sqrt(length(ourrows))

vals[i,3]=vals[i,2]-2*vv
vals[i,4]=vals[i,2]+2*vv
vals2[i,3]=vals2[i,2]-2*vv2
vals2[i,4]=vals2[i,2]+2*vv2
vals3[i,3]=vals3[i,2]-2*vv3
vals3[i,4]=vals3[i,2]+2*vv3


valsx[i,2]=mean((www[ourrows2]))
vals2x[i,2]=mean((www2[ourrows2]))
vals3x[i,2]=mean((www[ourrows2]/www2[ourrows2]))

valsx[i,1]=(mean(www3[ourrows2]))
vals2x[i,1]=valsx[i,1]
vals3x[i,1]=valsx[i,1]
vv=sd((www[ourrows2]))/sqrt(length(ourrows2))
vv2=sd((www2[ourrows2]))/sqrt(length(ourrows2))
vv3=sd((www[ourrows2]/www2[ourrows2]))/sqrt(length(ourrows2))

valsx[i,3]=valsx[i,2]-2*vv
valsx[i,4]=valsx[i,2]+2*vv
vals2x[i,3]=vals2x[i,2]-2*vv2
vals2x[i,4]=vals2x[i,2]+2*vv2
vals3x[i,3]=vals3x[i,2]-2*vv3
vals3x[i,4]=vals3x[i,2]+2*vv3

}
vals3[,1]=vals2[,1]

vals3x[,1]=vals2x[,1]
valsx=valsx[1:(length(bins2)-1),]
vals2x=vals2x[1:(length(bins2)-1),]
vals3x=vals3x[1:(length(bins2)-1),]


colormap=heat.colors(100)

color=1+floor((b6[,"enrichment"]-0.5)*75/1.5)
color[color>75]=75

par(mfrow=c(1,2))
plot(www3[cond2],(www[cond2]),xlab="SPO11",ylab="DMC1 (WT)",pch=19,cex=0.1,col=5,ylim=c(0,15))
points(www3[cond],(www[cond]),xlab="SPO11",ylab="DMC1 (WT)",pch=19,cex=0.1,col=colormap[color[cond]])

points(vals[,1],vals[,2],pch=19)
segments(x0=vals[,1],x1=vals[,1],y0=vals[,3],y1=vals[,4])
points(valsx[,1],valsx[,2],pch=19,col=4)
segments(x0=valsx[,1],x1=valsx[,1],y0=valsx[,3],y1=valsx[,4],col=4)




plot(www3[cond2],(www2[cond2]),xlab="SPO11",ylab="DMC1 (WT)",pch=19,cex=0.1,col=5)
points(www3[cond],(www2[cond]),xlab="SPO11",ylab="DMC1 (WT)",pch=19,cex=0.1,col=colormap[color[cond]])
points(www3[cond],(www2[cond]),xlab="SPO11",ylab="DMC1 (ZCWPW1-null)",pch=19,cex=0.1,col=colormap[color[cond]])
points(vals2[,1],vals2[,2],pch=19)
segments(x0=vals2[,1],x1=vals2[,1],y0=vals2[,3],y1=vals2[,4])
points(valsx[,1],vals2x[,2],pch=19,col=4)
segments(x0=vals2x[,1],x1=vals2x[,1],y0=vals2x[,3],y1=vals2x[,4],col=4)

dev.copy2pdf(file="xydmc1vsspo11wtandnull.pdf",width=10,height=5)


plot(www3[cond],(www2[cond]),xlab="SPO11",ylab="DMC1 (ZCWPW1-null)",pch=19,cex=0.1,col=colormap[color[cond]])
points(vals2[,1],vals2[,2],pch=19)
segments(x0=vals2[,1],x1=vals2[,1],y0=vals2[,3],y1=vals2[,4])


###hotspot sharing plot (code lost!)
###looking at overlaps.
>99% of denovo hotspots are B6 ones (within 1kb of centre of a B6 hotspot)

###then stratify B6 hotspots from weakest to strongest (x-axis is heat relative to strongest) and divide hotter ones into whether force-calling in ZCWPW1 shows evidence, or not in terms of p<0.05, p<0.001,or some reads (heat>0).

####autosomal only but X-chrom. plot appears similar.

####bars on top are individual hotspot heats so can see many hotspots are being counted (and all of top several hundred are re-seen at p<0.001).





