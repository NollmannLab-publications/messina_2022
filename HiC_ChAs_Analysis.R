## 3D chromatin interactions involving Drosophila insulators are infrequent but preferential and arise before TADs and transcription


##############
##  Fig 1B  ##
##############


features_dir="/Users/FR/Documents/Messina_et_al/Features_datasets/"
network_dir="/Users/FR/Documents/Messina_et_al/Network_datasets/"


### Load libraries
library(devtools)
#devtools::install_bitbucket("eraineri/chaser",  build_vignettes=TRUE) # install ChAseR package
library(chaser)


### Load ChIP-Seq peaks
beaf=read.table(paste0(features_dir,"BEAF32_summits_GSM1535963.bed"))
cbp=read.table(paste0(features_dir,"CBP_Kc167_WT_Li_dm6_peaks.narrowPeak"))
chro=read.table(paste0(features_dir,"Chromator_Kc167_WT_Li_dm6_peaks.narrowPeak"))
cp190=read.table(paste0(features_dir,"CP190_dm6.bed"))
ctcf=read.table(paste0(features_dir,"CTCF_N_dm6.bed"))
dref=read.table(paste0(features_dir,"DREF_Kc167_WT_Li_dm6_peaks.narrowPeak"))
fs1h=read.table(paste0(features_dir,"Fs1h_Kc167_WT_Li_dm6_peaks.narrowPeak"))
gaf=read.table(paste0(features_dir,"GAF_dm6.bed"))
l3mbt=read.table(paste0(features_dir,"L3mbt_Kc167_WT_Li_dm6_peaks.narrowPeak"))
mdg4=read.table(paste0(features_dir,"MDG4_dm6.bed"))
pita=read.table(paste0(features_dir,"Pita_S2_WT_Maksimenko_dm6_peaks.narrowPeak"))
suhw=read.table(paste0(features_dir,"suHw_dm6.bed"))
z4=read.table(paste0(features_dir,"Z4_Kc167_WT_Li_dm6_peaks.narrowPeak"))
zipic=read.table(paste0(features_dir,"ZIPIC_S2_WT_Maksimenko_dm6_peaks.narrowPeak"))
zw5=read.table(paste0(features_dir,"ZW5_summits_dm6_GSM2042227.bed"))


### Load network
net=read.table(paste0(network_dir,"nc14_detect_res7.tsv"),header = T)
net$chrom1=paste0("chr",net$chrom1);net$chrom2=paste0("chr",net$chrom2)
net=make_chromnet(net[,1:6],missing=NA)


### Load features
net <- load_features(net, beaf, type="bed3",missingv=0,auxfun="mean",featnames = "BEAF-32")
net <- load_features(net, cbp, type="bed3",missingv=0,auxfun="mean",featnames = "CBP")
net <- load_features(net, chro, type="bed3",missingv=0,auxfun="mean",featnames = "CHRO")
net <- load_features(net, cp190, type="bed3",missingv=0,auxfun="mean",featnames = "CP190")
net <- load_features(net, ctcf, type="bed3",missingv=0,auxfun="mean",featnames = "dCTCF")
net <- load_features(net, dref, type="bed3",missingv=0,auxfun="mean",featnames = "DREF")
net <- load_features(net, fs1h, type="bed3",missingv=0,auxfun="mean",featnames = "Fs(1)h")
net <- load_features(net, gaf, type="bed3",missingv=0,auxfun="mean",featnames = "GAF")
net <- load_features(net, l3mbt, type="bed3",missingv=0,auxfun="mean",featnames = "L(3)MBT")
net <- load_features(net, mdg4, type="bed3",missingv=0,auxfun="mean",featnames = "Mod(mdg4)")
net <- load_features(net, pita, type="bed3",missingv=0,auxfun="mean",featnames = "Pita")
net <- load_features(net, suhw, type="bed3",missingv=0,auxfun="mean",featnames = "SU(HW)")
net <- load_features(net, z4, type="bed3",missingv=0,auxfun="mean",featnames = "Z4")
net <- load_features(net, zipic, type="bed3",missingv=0,auxfun="mean",featnames = "ZIPIC")
net <- load_features(net, zw5, type="bed3",missingv=0,auxfun="mean",featnames = "Zw5")


### Compute ChAs
chas_net <- chas(net,"chas")


### Compute Randomisations ChAs
random_net=randomize(net,1000,dist.match=T)
rm(chas_random_net)
chas_random_net=vector()
for (i in 1:length(random_net)){
  chas_random_net=c(chas_random_net,chas(random_net[[i]],"chas"))
}

rm(df_random_net,df_chas_net)
df_random_net=data.frame(Value=chas_random_net,Sample=names(chas_random_net))
df_chas_net=data.frame(Value=chas_net,Sample=names(chas_net))


### Compute Z-Scores
nb_features=length(df_chas_net$Sample)
mean_random <- aggregate(. ~ Sample, data = df_random_net, mean)
sd_random <- aggregate(. ~ Sample, data = df_random_net, sd)

rm(i,zscores)
zscores=vector()
for(i in mean_random$Sample){
  zscores=c(zscores,(df_chas_net[which(df_chas_net$Sample==i),1]-mean_random[which(mean_random$Sample==i),2])/sd_random[which(sd_random$Sample==i),2])
}
names(zscores)=mean_random$Sample


### Make the ChAs Z-Scores Figure
published_zscores=c(9.2509149,0.2834827,7.6957085,1.5531658,1.7042579,4.0366287,1.3135821,0.1754332,4.2655178,0.6520337,2.6650313,0.6622692,8.3708116,2.3071721,3.2684808)
names(published_zscores)=c("BEAF-32","CBP","CHRO","CP190","dCTCF","DREF","Fs(1)h","GAF","L(3)MBT","Mod(mdg4)","Pita","SU(HW)","Z4","ZIPIC","Zw5")
all_colors=colorRampPalette(c("#016600","#00cc00"))(17)
all_colors=ifelse(published_zscores<2,"darkgrey",all_colors)

b=barplot(published_zscores,col=all_colors, ylab="",main="",ylim=c(0,20),xaxt="n",axes=0)
text(b, published_zscores, labels=as.character(round(published_zscores,1)), cex= 1.1,pos=3)
mtext("ChAs Z-Score", side = 2, line = 2,cex=1.4)
axis(2, at=c(0,20),las=2,cex.axis=1.4,lwd=1.2,font=1)
abline(h=2,lty=3,col="darkred",lwd=1.5)
text(cex=1.2, x=b+.4, y=-1.5, names(published_zscores), xpd=T, srt=45,pos=2)






##############
##  Fig 2B  ##
##############


features_dir="/Users/FR/Documents/Messina_et_al/Features_datasets/"
network_dir="/Users/FR/Documents/Messina_et_al/Network_datasets/"


### Load libraries
library(devtools)
#devtools::install_bitbucket("eraineri/chaser",  build_vignettes=TRUE) # install ChAseR package
library(chaser)


### Load TAD boundaries
tad=read.table(paste0(features_dir,"consensus_boundaries_5kb.tsv"),header=T)[,c(1:3,6)]
tad$chrom=paste0("chr",tad$chrom)

### Load Network
net=read.table(paste0(network_dir,"nc14_detect_res7.tsv"),header = T)
net$chrom1=paste0("chr",net$chrom1);net$chrom2=paste0("chr",net$chrom2)
net=make_chromnet(net[,1:6],missing=NA)


### Load Features
net <- load_features(net, tad, type="bed3",missingv=0,auxfun="mean",featnames = "TAD borders")


### Compute ChAs
chas_net <- chas(net,"chas")


### Compute Randomisations ChAs
random_net=randomize(net,1000,dist.match=T)
rm(chas_random_net)
chas_random_net=vector()
for (i in 1:length(random_net)){
  chas_random_net=c(chas_random_net,chas(random_net[[i]],"chas"))
}

rm(df_random_net,df_chas_net)
df_random_net=data.frame(Value=chas_random_net,Sample=names(chas_random_net))
df_chas_net=data.frame(Value=chas_net,Sample=names(chas_net))


### Compute Z-Scores
nb_features=length(df_chas_net$Sample)
mean_random <- aggregate(. ~ Sample, data = df_random_net, mean)
sd_random <- aggregate(. ~ Sample, data = df_random_net, sd)

rm(i,zscores)
zscores=vector()
for(i in mean_random$Sample){
  zscores=c(zscores,(df_chas_net[which(df_chas_net$Sample==i),1]-mean_random[which(mean_random$Sample==i),2])/sd_random[which(sd_random$Sample==i),2])
}
names(zscores)=mean_random$Sample


### Make the TAD borders ChAs Figure
plot(density(df_random_net$Value),ylim=c(1.3,33.8),xlim=c(-0.043,0.143),lwd=4,ylab="",main="",xlab="",xaxt="n",yaxt="n")
mtext("Density", side = 2, line = 2,cex=1.4)
axis(2, at=c(0,35),las=2,cex.axis=1.4,lwd=1.2,font=1)
mtext("TAD borders ChAs", side = 1, line = 3,cex=1.4)
axis(1, at=c(-0.05,0.05,0.15),las=1,cex.axis=1.4,lwd=1.2,font=1)
abline(v=df_chas_net$Value,lty=2,col="blue",lwd=2)
legend("topright", legend=c("Randomized", "TAD borders"),col=c("black", "blue"), lty=1:2, cex=0.8,text.font=2, bg="lightgrey",lwd=2)






###################
##  Supp Fig 1A  ##
###################


features_dir="/Users/FR/Documents/Messina_et_al/Features_datasets/"
network_dir="/Users/FR/Documents/Messina_et_al/Network_datasets/"


### Load libraries
library(devtools)
#devtools::install_bitbucket("eraineri/chaser",  build_vignettes=TRUE) # install ChAseR package
library(chaser)

### Load ChIP-Seq peaks
rad21=read.table(paste0(features_dir,"RAD21_Kc167_WT_Li_dm6_peaks.narrowPeak"))
pc=read.table(paste0(features_dir,"PC_GSE60428_4_12h_peaks.bed"))
ph=read.table(paste0(features_dir,"PH_GSE60428_4_12h_peaks.bed"))
s5p=read.table(paste0(features_dir,"PolII_pSer5_nc14L_WT_Blythe_dm6_peaks.narrowPeak"))
zelda=read.table(paste0(features_dir,"ZLD_3h_WT_Harrison_dm6_noinput_peaks.narrowPeak"))


### Load Network
net=read.table(paste0(network_dir,"nc14_detect_res7.tsv"),header = T)
net$chrom1=paste0("chr",net$chrom1);net$chrom2=paste0("chr",net$chrom2)
net=make_chromnet(net[,1:6],missing=NA)


#Load Features
net <- load_features(net, rad21, type="bed3",missingv=0,auxfun="mean",featnames = "Rad21")
net <- load_features(net, pc, type="bed3",missingv=0,auxfun="mean",featnames = "Pc")
net <- load_features(net, ph, type="bed3",missingv=0,auxfun="mean",featnames = "Ph")
net <- load_features(net, s5p, type="bed3",missingv=0,auxfun="mean",featnames = "S5P")
net <- load_features(net, zelda, type="bed3",missingv=0,auxfun="mean",featnames = "Zelda")


### Compute ChAs
chas_net <- chas(net,"chas")


### Compute Randomisations ChAs
random_net=randomize(net,1000,dist.match=T)
rm(chas_random_net)
chas_random_net=vector()
for (i in 1:length(random_net)){
  chas_random_net=c(chas_random_net,chas(random_net[[i]],"chas"))
}

rm(df_random_net,df_chas_net)
df_random_net=data.frame(Value=chas_random_net,Sample=names(chas_random_net))
df_chas_net=data.frame(Value=chas_net,Sample=names(chas_net))


### Compute Z-Scores
nb_features=length(df_chas_net$Sample)
mean_random <- aggregate(. ~ Sample, data = df_random_net, mean)
sd_random <- aggregate(. ~ Sample, data = df_random_net, sd)

rm(i,zscores)
zscores=vector()
for(i in df_chas_net$Sample){
  zscores=c(zscores,(df_chas_net[which(df_chas_net$Sample==i),1]-mean_random[which(mean_random$Sample==i),2])/sd_random[which(sd_random$Sample==i),2])
}
names(zscores)=df_chas_net$Sample


### Make the ChAs Z-Scores Figure
published_zscores=c(0.6853735,5.7036930,2.2631745,15.5249395,16.3786160) 
names(published_zscores)=c("Rad21","Pc","Ph","PS5","Zelda") 
all_colors=c("#ff00ff","#0099ff","#1300ff","#ff4c4c","#fc0000")

b=barplot(published_zscores,col=all_colors, ylab="",main="",ylim=c(0,20),xaxt="n",axes=0)
text(b, published_zscores, labels=as.character(round(published_zscores,1)), cex= 1.1,pos=3)
mtext("ChAs Z-Score", side = 2, line = 2,cex=1.4)
axis(2, at=c(0,20),las=2,cex.axis=1.4,lwd=1.2,font=1)
abline(h=2,lty=3,col="darkred",lwd=1.5)
text(cex=1.2, x=b+.4, y=-1.5, names(published_zscores), xpd=T, srt=45,pos=2)





