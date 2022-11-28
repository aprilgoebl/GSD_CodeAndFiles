## Great Sand Dune (GSD), Helianthus petiolaris, reciprocal transplant data analysis script
## Script that uses GBS (genotyping-by-sequencing) data from plants in 2012 GSD reciprocal transplant
## Uses PCA to look at genetic structure of samples.  
## Looks at relationships between PC loadings, allele frequency (AF) differences between ecotypes, 
## and AF change during the experiment.
## Assesses these relationships in inversion loci and non-inversion loci. 
## Plots the amount of AF change at each locus as a function of genomic position. 




rm(list=ls())
dev.off()
def_par = par()
windowsFonts(A = windowsFont("Times New Roman"))





## LOAD PACKAGES AND FUNCTIONS --------------------------------------------------------------------
library(pcadapt)
library(EnvStats)
library(stringr)
library(dplyr)
library(tidyr)
library(plotrix)
library(svglite)
## ------------------------------------------------------------------------------------------------




## READ IN AND FORMAT DATA -----------------------------------------------------------------------------

## Full VCF table for PCA analysis
path_to_vcf <- "KOP3_4_5.Ha412HOv2.filt.noleafSmpls.vcf" 


## Trait data (seed number data used; obtained from 10.5061/dryad.223p4)
datLate <- read.table(file ='2012_RT_Plants_full.csv', sep=',', header = TRUE)
datLate <- datLate[datLate$suspected_volunteer != 'Y',]       #Remove suspected volunteers 
## When zero heads collected, but plant had inflorescences, change zero seed counts to NA  
datLate$n.seeds.counted[datLate$n.seeds.counted==0 & datLate$heads.counted==0] <- NA


## Genotypes (012 format) for calculating allele frequencies
## Obtained from full vcf table above using vcftools
smpl_gts <- read.csv("KOP3_4_5.Ha412HOv2.filt.gts.csv", header=TRUE, stringsAsFactors=FALSE)    


## List of chromosome and genome positions for each SNP
pos <- read.table("KOP3_4_5.Ha412HOv2.filt.pos")       
colnames(pos) <- c("CHROM", "POS")
pos$CHROM <- as.character(pos$CHROM)


## List of chromosomal inversion start and end positions from Huang et al, 2020
invs <- read.table("InversionRanges_Ha412HOv2_Huang2020.txt", header=TRUE, fileEncoding="UCS-2LE")
invs$CHROM <- as.character(invs$CHROM)  


## VCF tables based only on inversion loci (obtained by filtering full vcf table above using vcftools)
## For analyzing each inversion as a single locus 
path_to_file.5 <- "KOP3_4_5.Ha412HOv2.filt.noLeafSmpls.invPet05.01.vcf"
path_to_file.7 <- "KOP3_4_5.Ha412HOv2.filt.noLeafSmpls.invPet07.01.vcf"
path_to_file.9.1 <- "KOP3_4_5.Ha412HOv2.filt.noLeafSmpls.invPet09.01.vcf"
path_to_file.9.2 <- "KOP3_4_5.Ha412HOv2.filt.noLeafSmpls.invPet09.02.vcf"
path_to_file.11 <- "KOP3_4_5.Ha412HOv2.filt.noLeafSmpls.invPet11.01.vcf"
path_to_file.14 <- "KOP3_4_5.Ha412HOv2.filt.noLeafSmpls.invPet14.01.vcf"
path_to_file.17 <- "KOP3_4_5.Ha412HOv2.filt.noLeafSmpls.invPet17.01.vcf"
## ----------------------------------------------------------------------------------------------------






## PCA ANALYSIS TO LOOK AT STURCTURE OF PRE- AND POST-SELECTION SAMPLES ------------------------------ 
vcf_tbl <- read.pcadapt(path_to_vcf, type="vcf")                            #Read in data from vcf table
x <- pcadapt(vcf_tbl, K=5)                                                  #Run PCAdapt 


## Obtain list of samples  
smpls <- smpl_gts %>% dplyr::select(FILE, PLT_CODE, POP_ID, SOURCE_ID)

## Remove wild, non-experimenal samples from list
smpls <- smpls[smpls$SOURCE_ID!="post_leafdune_d" & smpls$SOURCE_ID!="post_leafnon_n",]


## Visualize scree plot and PCA
plot(x, option="screeplot")                              #Scree plot
plot(x,option="scores", pop=smpls$SOURCE_ID)             #PC1 and PC2

PC1.perctVar <- round(((x$singular.values[1])^2) * 100, digit=2)
PC2.perctVar <- round(((x$singular.values[2])^2) * 100, digit=2)


## Extract PCA scores  
xScores <- x$scores                                  

dfScores <- as.data.frame(xScores[,1:2])             #Isolate scores for PC1 and PC2 only
colnames(dfScores) <- c("PC1score","PC2score")
dfScores$SOURCE_ID <- smpls$SOURCE_ID                #Add a column with source-selection ID
dfScores$POP_ID <- smpls$POP_ID                      #Add a column with type/pop ID
dfScores$PLT_CODE <- smpls$PLT_CODE                  #Add a column with unique plant code
dfScores <- dfScores[!is.na(dfScores$SOURCE_ID),]    #Remove rows with data for control plts



## Calculate mean PC scores for each source-selection group
mean_PCs <- aggregate(PC1score ~ SOURCE_ID, dfScores, mean)    #Calculate PC1 mean of each category 
mean_PC2 <- aggregate(PC2score ~ SOURCE_ID, dfScores, mean)    #Calculate PC2 mean of each category
mean_PCs$PC2score <- mean_PC2$PC2score                         #Combine mean values for PC1 and PC2



## Replace mean pre-hybrid PC scores after performing random sampling of hybrid types
## This will allow equal represenatation of genotypes from each hybrid type 
## including BCN which is missing from pre selection samples

## Combine different non-dune populations into single group
dfScores$POP_ID[grepl("pre_nondune1", dfScores$POP_ID)] <- "pre_nondune" 
dfScores$POP_ID[grepl("pre_nondune2", dfScores$POP_ID)] <- "pre_nondune" 
dfScores$POP_ID[grepl("pre_nondune3", dfScores$POP_ID)] <- "pre_nondune" 

## Loop over plant types of interest and select random smpl in designated proportions
type.list <- c("pre_BCD","pre_F1D","pre_F1N","pre_nondune")  #Plt types of interest
type.props <- c(20, 25, 25, 10)                              #Designated proportions
num.PC <- 2                                                  #Num of PC axes of interest
reps <- 1000
store.means <- as.data.frame(matrix(NA, reps, (num.PC)))
colnames(store.means) <- c("PC1","PC2")

for (rr in 1:reps) {
    store.scores <- NULL
    
    for (tt in 1:length(type.list)) {
        type <- subset(dfScores, POP_ID==type.list[tt])
        type$RAND <- runif(nrow(type), 0, 1)
        type.ord <- type[order(type$RAND),]
        type.rand <- type.ord[1:type.props[tt],]
        store.scores <- as.data.frame(rbind(store.scores, type.rand))
    }
    
    store.means[rr,1] <- mean(store.scores$PC1score)
    store.means[rr,2] <- mean(store.scores$PC2score)
}

## Replace pre-hyb PC scores
mean_PCs$PC1score[grepl("pre_hyb", mean_PCs$SOURCE_ID)] <- mean(store.means$PC1)  
mean_PCs$PC2score[grepl("pre_hyb", mean_PCs$SOURCE_ID)] <- mean(store.means$PC2)


# Calculate standard error of mean 
calcSE <- function(x){sd(x)/sqrt(length(x))}
PC1_SE <- aggregate(PC1score ~ SOURCE_ID, dfScores, calcSE)
PC2_SE <- aggregate(PC2score ~ SOURCE_ID, dfScores, calcSE)


## Calculate CIs for pre-selection samples
preDune <- subset(dfScores, SOURCE_ID=="pre_dune")
errPC1.preDune <- t.test(preDune$PC1score, conf.level=0.95)$conf.int
errPC2.preDune <- t.test(preDune$PC2score, conf.level=0.95)$conf.int

preNon <- subset(dfScores, SOURCE_ID=="pre_non")
errPC1.preNon <- t.test(preNon$PC1score, conf.level=0.95)$conf.int
errPC2.preNon <- t.test(preNon$PC2score, conf.level=0.95)$conf.int

errPC1.preHyb <- t.test(store.means$PC1, conf.level = 0.95)$conf.int
errPC2.preHyb <- t.test(store.means$PC2, conf.level = 0.95)$conf.int




## Generate figures ---------------------------------------------------------------- 

## Colours
colD <- rgb(t(col2rgb("orange")),alpha=230,maxColorValue = 255)
colH <- rgb(t(col2rgb("#33a02c")),alpha=220,maxColorValue = 255)
colN <- rgb(t(col2rgb("royalblue")),alpha=220,maxColorValue = 255)
cols <- c(colD,colH,colN)

colD_bg <- rgb(t(col2rgb("orange")),alpha=70,maxColorValue = 255)
colH_bg <- rgb(t(col2rgb("#33a02c")),alpha=70,maxColorValue = 255)
colN_bg <- rgb(t(col2rgb("royalblue")),alpha=70,maxColorValue = 255)
cols_bg <- c(colD_bg,colH_bg,colN_bg)


symbs_pre <- 17
symbs <- c(19,19,19)
par(mar=c(5,6,3,4))  


## Figure S6
## Plot pre-selection samples
par(mfrow=c(1,1)) 

plot(NA, NA, ylim=c(-0.13, 0.12), xlim=c(-0.15, 0.1), xlab=paste("PC1", " ", "(",PC1.perctVar,"%)", sep=""), 
     ylab=paste("PC2"," ", "(",PC2.perctVar,"%)", sep=""),cex.axis=1.25, cex.lab=2, family="A")
points(mean_PCs$PC1score[7], mean_PCs$PC2score[7], col=cols[1], pch=symbs_pre, cex=3)
points(mean_PCs$PC1score[8], mean_PCs$PC2score[8], col=cols[2], pch=symbs_pre, cex=3)
points(mean_PCs$PC1score[9], mean_PCs$PC2score[9], col=cols[3], pch=symbs_pre, cex=3)

## Add error bars (95% CIs)
arrows(errPC1.preDune[1], mean_PCs$PC2score[7], errPC1.preDune[2], mean_PCs$PC2score[7],
       col=cols[1],length=0.05,angle=90,code=3)
arrows(mean_PCs$PC1score[7],errPC2.preDune[1], mean_PCs$PC1score[7], errPC2.preDune[2], 
       col=cols[1],length=0.05,angle=90,code=3)

arrows(errPC1.preHyb[1], mean_PCs$PC2score[8], errPC1.preHyb[2], mean_PCs$PC2score[8],
       col=cols[2],length=0.05,angle=90,code=3)
arrows(mean_PCs$PC1score[8],errPC2.preHyb[1], mean_PCs$PC1score[8], errPC2.preHyb[2], 
       col=cols[2],length=0.05,angle=90,code=3)

arrows(errPC1.preNon[1], mean_PCs$PC2score[9], errPC1.preNon[2], mean_PCs$PC2score[9],
       col=cols[3],length=0.05,angle=90,code=3)
arrows(mean_PCs$PC1score[9],errPC2.preNon[1], mean_PCs$PC1score[9], errPC2.preNon[2], 
       col=cols[3],length=0.05,angle=90,code=3)

source_names <- c("pre_dune","pre_hyb","pre_non")

for (i in 1:length(source_names)) {
    points(dfScores$PC1score[dfScores$SOURCE_ID==source_names[i]], 
           dfScores$PC2score[dfScores$SOURCE_ID==source_names[i]], 
           col=cols_bg[i], pch=symbs_pre,cex=1.25)
}

legend("topleft", c("Dune source", "Hybrid source", "Non-dune source"), col=cols, pch=symbs_pre,cex=1.75,
       inset=c(0,0), xpd=TRUE, horiz=FALSE, bty="y")
## ----------------------------------------------------------------------



## Figure 3
## Plot post-selection in dunes
pdf('GSD_fig3_221111.pdf', width=10, height=8)
par(mfrow=c(2,1))  #Plot the following two PCAs in 2 rows and 1  column
par(mar=c(6,5,2.6,2)) 


plot(NA, NA, ylim=c(-0.1, 0.1), xlim=c(-0.09, 0.09), xlab=NA, ylab="PC2 (2.0%)",cex.axis=1.25,cex.lab=1.5)
title("Dune habitat", line = 0.3, cex.main=1.75, font.main=6)
text(-0.115,0.14, "A", cex=2,xpd=TRUE,font=1.75)
legend("bottom", c("Dune source", "Hybrid source", "Non-dune source"), col=cols, pch=symbs,cex=1.25,
       inset=c(0,-0.58), xpd=TRUE, horiz=TRUE, bty="n")

points(mean_PCs$PC1score[1], mean_PCs$PC2score[1], col=cols[1], pch=symbs[1], cex=2.5) #post_dune_d
points(mean_PCs$PC1score[3], mean_PCs$PC2score[3], col=cols[2], pch=symbs[2], cex=2.5) #post_hyb_d
points(mean_PCs$PC1score[5], mean_PCs$PC2score[5], col=cols[3], pch=symbs[3], cex=2.5) #post_non_d

## Add standard error bars 
plotCI(x=mean_PCs$PC1score[1], y=mean_PCs$PC2score[1], uiw=PC1_SE[1,2],err="x", col=cols, add=T)
plotCI(x=mean_PCs$PC1score[1], y=mean_PCs$PC2score[1], uiw=PC2_SE[1,2],err="y", col=cols, add=T)
plotCI(x=mean_PCs$PC1score[3], y=mean_PCs$PC2score[3], uiw=PC1_SE[4,2],err="x", col=cols[2], add=T)
plotCI(x=mean_PCs$PC1score[3], y=mean_PCs$PC2score[3], uiw=PC2_SE[4,2],err="y", col=cols[2], add=T)
plotCI(x=mean_PCs$PC1score[5], y=mean_PCs$PC2score[5], uiw=PC1_SE[7,2],err="x", col=cols[3], add=T)
plotCI(x=mean_PCs$PC1score[5], y=mean_PCs$PC2score[5], uiw=PC2_SE[7,2],err="y", col=cols[3], add=T)

postd_names <- c("post_dune_d","post_hyb_d","post_non_d")

for (i in 1:length(postd_names)) {
    points(dfScores$PC1score[dfScores$SOURCE_ID==postd_names[i]], 
           dfScores$PC2score[dfScores$SOURCE_ID==postd_names[i]], 
           col=cols_bg[i], pch=symbs[1],cex=1)
}

## Add arrows showing direction and magnitide of PC1 shift
arrows(mean_PCs$PC1score[7], 0.06, mean_PCs$PC1score[1], 0.06, angle=30,code=2, col=colD,lwd=2.25,length=0.15)
arrows(mean_PCs$PC1score[8], 0.06, mean_PCs$PC1score[3], 0.06, angle=30,code=2, col=colH,lwd=2.25,length=0.15)
arrows(mean_PCs$PC1score[9], 0.06, mean_PCs$PC1score[5], 0.06, angle=30,code=2, col=colN,lwd=2.25,length=0.15)


## ----------- Add pre-selection mean points 
points(mean_PCs$PC1score[7], mean_PCs$PC2score[7], col=cols[1], pch=17, cex=2.5)
points(mean_PCs$PC1score[8], mean_PCs$PC2score[8], col=cols[2], pch=17, cex=2.5)
points(mean_PCs$PC1score[9], mean_PCs$PC2score[9], col=cols[3], pch=17, cex=2.5)

## Add error bars (95% CIs)
arrows(errPC1.preDune[1], mean_PCs$PC2score[7], errPC1.preDune[2], mean_PCs$PC2score[7],
       col=cols[1],length=0.05,angle=90,code=3)
arrows(mean_PCs$PC1score[7],errPC2.preDune[1], mean_PCs$PC1score[7], errPC2.preDune[2], 
       col=cols[1],length=0.05,angle=90,code=3)

arrows(errPC1.preHyb[1], mean_PCs$PC2score[8], errPC1.preHyb[2], mean_PCs$PC2score[8],
       col=cols[2],length=0.05,angle=90,code=3)
arrows(mean_PCs$PC1score[8],errPC2.preHyb[1], mean_PCs$PC1score[8], errPC2.preHyb[2], 
       col=cols[2],length=0.05,angle=90,code=3)

arrows(errPC1.preNon[1], mean_PCs$PC2score[9], errPC1.preNon[2], mean_PCs$PC2score[9],
       col=cols[3],length=0.05,angle=90,code=3)
arrows(mean_PCs$PC1score[9],errPC2.preNon[1], mean_PCs$PC1score[9], errPC2.preNon[2], 
       col=cols[3],length=0.05,angle=90,code=3)
## ----------------



## Plot post-selection in nondunes 
plot(NA, NA, ylim=c(-0.1, 0.1), xlim=c(-0.09, 0.09), xlab="PC1 (4.2%)", ylab="PC2 (2.0%)",cex.axis=1.25,cex.lab=1.5)
title("Non-dune habitat", line = 0.3, cex.main=1.75, font.main=6)

points(mean_PCs$PC1score[2], mean_PCs$PC2score[2], col=cols[1], pch=symbs[1], cex=2.5) #post_dune_n
points(mean_PCs$PC1score[4], mean_PCs$PC2score[4], col=cols[2], pch=symbs[2], cex=2.5) #post_hyb_n
points(mean_PCs$PC1score[6], mean_PCs$PC2score[6], col=cols[3], pch=symbs[3], cex=2.5) #post_non_n

## Add standard error bars
plotCI(x=mean_PCs$PC1score[2], y=mean_PCs$PC2score[2], uiw=PC1_SE[2,2],err="x", col=cols, add=T)
plotCI(x=mean_PCs$PC1score[2], y=mean_PCs$PC2score[2], uiw=PC2_SE[2,2],err="y", col=cols, add=T)
plotCI(x=mean_PCs$PC1score[4], y=mean_PCs$PC2score[4], uiw=PC1_SE[5,2],err="x", col=cols[2], add=T)
plotCI(x=mean_PCs$PC1score[4], y=mean_PCs$PC2score[4], uiw=PC2_SE[5,2],err="y", col=cols[2], add=T)
plotCI(x=mean_PCs$PC1score[6], y=mean_PCs$PC2score[6], uiw=PC1_SE[8,2],err="x", col=cols[3], add=T)
plotCI(x=mean_PCs$PC1score[6], y=mean_PCs$PC2score[6], uiw=PC2_SE[8,2],err="y", col=cols[3], add=T)

postn_names <- c("post_dune_n","post_hyb_n","post_non_n")

for (i in 1:length(postn_names)) {
    points(dfScores$PC1score[dfScores$SOURCE_ID==postn_names[i]], 
           dfScores$PC2score[dfScores$SOURCE_ID==postn_names[i]], 
           col=cols_bg[i], pch=symbs[1],cex=1)
}

## Add arrows showing direction and magnitide of PC1 shift
arrows(mean_PCs$PC1score[7], 0.05, mean_PCs$PC1score[2], 0.05, angle=30,code=2, col=colD,lwd=2.25,length=0.15)
arrows(mean_PCs$PC1score[8], 0.05, mean_PCs$PC1score[4], 0.05, angle=30,code=2, col=colH,lwd=2.25,length=0.15)
arrows(mean_PCs$PC1score[9], 0.05, mean_PCs$PC1score[6], 0.05, angle=30,code=2, col=colN,lwd=2.25,length=0.15)



## ----------- Add pre-selection mean points 
points(mean_PCs$PC1score[7], mean_PCs$PC2score[7], col=cols[1], pch=17, cex=2.5)
points(mean_PCs$PC1score[8], mean_PCs$PC2score[8], col=cols[2], pch=17, cex=2.5)
points(mean_PCs$PC1score[9], mean_PCs$PC2score[9], col=cols[3], pch=17, cex=2.5)

## Add error bars (95% CIs)
arrows(errPC1.preDune[1], mean_PCs$PC2score[7], errPC1.preDune[2], mean_PCs$PC2score[7],
       col=cols[1],length=0.05,angle=90,code=3)
arrows(mean_PCs$PC1score[7],errPC2.preDune[1], mean_PCs$PC1score[7], errPC2.preDune[2], 
       col=cols[1],length=0.05,angle=90,code=3)

arrows(errPC1.preHyb[1], mean_PCs$PC2score[8], errPC1.preHyb[2], mean_PCs$PC2score[8],
       col=cols[2],length=0.05,angle=90,code=3)
arrows(mean_PCs$PC1score[8],errPC2.preHyb[1], mean_PCs$PC1score[8], errPC2.preHyb[2], 
       col=cols[2],length=0.05,angle=90,code=3)

arrows(errPC1.preNon[1], mean_PCs$PC2score[9], errPC1.preNon[2], mean_PCs$PC2score[9],
       col=cols[3],length=0.05,angle=90,code=3)
arrows(mean_PCs$PC1score[9],errPC2.preNon[1], mean_PCs$PC1score[9], errPC2.preNon[2], 
       col=cols[3],length=0.05,angle=90,code=3)

text(-0.115,0.14, "B", cex=2, xpd=TRUE,font=1.75)
text(-0.075, -0.17, "\"Nondune-like\"", cex=1.5,xpd=TRUE)
text(0.075, -0.17, "\"Dune-like\"", cex=1.5,xpd=TRUE)
dev.off()
## ---------------------------------------------------------------------------------------------------







## Figure 5 ------------------------------------------------------------------------------------------
## Fecundity adjusted PC1 scores

## Add column with combined groups (ie. dune_sel, nondune_sel, pre_sel)
dfScores$SELECTION[grepl("post_dune_d", dfScores$SOURCE_ID)] <- "dune_sel"
dfScores$SELECTION[grepl("post_hyb_d", dfScores$SOURCE_ID)] <- "dune_sel"
dfScores$SELECTION[grepl("post_non_d", dfScores$SOURCE_ID)] <- "dune_sel"

dfScores$SELECTION[grepl("post_dune_n", dfScores$SOURCE_ID)] <- "nondune_sel"
dfScores$SELECTION[grepl("post_hyb_n", dfScores$SOURCE_ID)] <- "nondune_sel"
dfScores$SELECTION[grepl("post_non_n", dfScores$SOURCE_ID)] <- "nondune_sel"

dfScores$SELECTION[grepl("pre_*", dfScores$SOURCE_ID)] <- "pre_sel"



## Calculate mean pre-selection PC1 score using designated proportions of hybrid, nondune, & dune
## Designated proportions equal the proportions of each source planted in the field
## This will also incorporate equal represenatation of genotypes from each hybrid type (as above) 

## Combine different dune populations into single group
dfScores$POP_ID[grepl("pre_dune1", dfScores$POP_ID)] <- "pre_dune" 
dfScores$POP_ID[grepl("pre_dune2", dfScores$POP_ID)] <- "pre_dune" 
dfScores$POP_ID[grepl("pre_dune3", dfScores$POP_ID)] <- "pre_dune" 

## Loop over pre-selection sources and select random sample in designated proportions
pre.list <- c("pre_BCD","pre_F1D","pre_F1N","pre_nondune", "pre_dune")  
type.props <- c(20, 25, 25, 70, 60)                          #Designated proportions
type.props <- round(type.props * (56/60))                    #Adjust proportions based on num dune (limiting smpl group)
num.PC <- 1                                                  #Num of PC axes of interest
reps <- 1000
store.means <- as.data.frame(matrix(NA, reps, num.PC))

for (rr in 1:reps) {
    store.scores <- NULL
    
    for (tt in 1:length(pre.list)) {
        type <- subset(dfScores, POP_ID==pre.list[tt])
        type$RAND <- runif(nrow(type), 0, 1)
        type.ord <- type[order(type$RAND),]
        type.rand <- type.ord[1:type.props[tt],]
        store.scores <- as.data.frame(rbind(store.scores, type.rand))
    }
    store.means[rr,1] <- mean(store.scores$PC1score)
}

## Calculate the mean PC1 for pre-selection samples
mean.preSel <- mean(store.means$V1)



## Calculate mean post-selection PC1 scores 
mean_PC1 <- aggregate(PC1score ~ SELECTION, dfScores, mean)

## Replace mean PC score for pre-selection group with value from random sampling
mean_PC1$PC1score[grepl("pre_sel", mean_PC1$SELECTION)] <- mean.preSel




## Combine trait data with PC scores 
scores.post <- dfScores[dfScores$SELECTION!="pre_sel",]   #Subset post-selection samples from scores data
    
## Modify unique plant code data type to be suitable for following conditional comparison 
scores.post$PLT_CODE <- as.character(scores.post$PLT_CODE)                                                            
datLate$plt_code <- as.character(datLate$plt_code)

## Loop over samples in scores dataframe and look for corresponding plant in trait dataframe   
scoresAndTraits <- as.data.frame(matrix(NA,nrow=nrow(scores.post),ncol=(ncol(scores.post)+ncol(datLate))))
colnames(scoresAndTraits) <- c(colnames(scores.post),colnames(datLate))

for (ss in 1:nrow(scores.post)) {
    for (pp in 1:nrow(datLate)) {
        if (scores.post$PLT_CODE[ss] == datLate$plt_code[pp]) {
            scoresAndTraits[ss,] <- cbind(scores.post[ss,], datLate[pp,])
        } 
    }
}




## Calculate new PC scores for post-selection smpls by weighting with fecundity  
scoresAndTraits$PC1scoreNEW <- scoresAndTraits$PC1score * scoresAndTraits$n.seeds.counted   #New PC score is product of old score and # of seeds produced 
sum_inds_NEW <- aggregate(n.seeds.counted ~ SELECTION, scoresAndTraits, sum)                #Total number of seeds produced by selection type
sum_PC1_NEW <- aggregate(PC1scoreNEW ~ SELECTION, scoresAndTraits, sum)                     #Sum new PC scores & divide by total seeds (below) to get average 
mean_PC1_NEW <- sum_PC1_NEW[,2] / sum_inds_NEW[,2]                                          #Calc new PC means by dividing score sum by seed sum for each sel type



## Generate figure
## Plot selection stage vs PC1 (Figure 5) 
pdf('GSD_fig5_221111.pdf', width=10, height=7.5)
par(mfrow=c(1,1))  
par(mar=c(6,6,2,3), xpd=TRUE)  

## Look at order of data
mean_PC1

y_axis <- c(3,2,1)
PC1_d <- c(mean_PC1$PC1score[3], mean_PC1$PC1score[1], mean_PC1_NEW[1])
PC1_nd <- c(mean_PC1$PC1score[3], mean_PC1$PC1score[2], mean_PC1_NEW[2])

plot(PC1_d[1],y_axis[1],cex=3,xlim=c(-0.05,0.055),ylim=c(0.5,3.7),pch=1,yaxt='n',xlab="PC1",
     cex.lab=2,col="black",ylab=NA,cex.axis=1.5)

arrows(PC1_d[1], y_axis[1], PC1_d[2], y_axis[2], angle=30,code=2,col="grey80",lwd=3.5,length=.15)
arrows(PC1_d[2], y_axis[2], PC1_d[3], y_axis[3], angle=30,code=2,col="grey80",lwd=3.5,length=.15,lty="twodash")
arrows(PC1_nd[1], y_axis[1], PC1_nd[2], y_axis[2], angle=30,code=2,col="black",lwd=3.5,length=.15)
arrows(PC1_nd[2], y_axis[2], PC1_nd[3], y_axis[3], angle=30,code=2,col="black",lwd=3.5,length=.15,lty="twodash")

text(-0.061,3, "Pre\nselection", cex=1.4,font=6)
text(-0.061,2, "Early life\nselection", cex=1.4,font=6)
text(-0.061,1, "Late life\nselection", cex=1.4,font=6)

legend("topright", legend=c("Dune habitat","Non-dune habitat"),col=c("grey80","black"), 
       cex=1.75, lwd=4,text.font=6)

text(-0.04, -0.05, "\"Nondune-like\"", cex=1.75,xpd=TRUE)
text(0.049, -0.05, "\"Dune-like\"", cex=1.75,xpd=TRUE)
dev.off()
## ------------------------------------------------------------------------------------------







## EXTRACT PC1 LOADINGS ---------------------------------------------------------------------
PC1_load <- x$loadings[,1]
PC1_load.abs <- abs(PC1_load)
## ------------------------------------------------------------------------------------------







## MAKE CONTINUOUS LIST OF GENOME POSITIONS -------------------------------------------------------
chroms <- as.character(unique(pos$CHROM))              #Make unique list of chromosome names
pos_on_genome <- NULL                                  #Make a storage variable for new positions
addVal  <- 0                                           #Create a variable to cumulatively add last values of each chromosome
spacer <- 10                                           #Add gap between chromosomes
adder <- max(pos$POS)*(spacer) / length(chroms)

for( i in 1:length(chroms) ){
    newPos <- adder + addVal + subset(pos, CHROM == chroms[i])$POS
    pos_on_genome <- c(pos_on_genome, newPos)
    addVal <- tail(pos_on_genome, n=1)
}
## ------------------------------------------------------------------------------------------------







## IDENTIFY SNPS THAT ARE LOCATED IN KNOWN CHROMOSOMAL INVERSIONS ---------------------------------
inv.sites <- as.data.frame(matrix(NA, nrow=nrow(pos), ncol=2))

for (ii in 1:nrow(pos)) {
    for (jj in 1:nrow(invs)) {
        if (pos$CHROM[ii] == invs$CHROM[jj] && pos$POS[ii] >= invs$POS_START[jj]
            && pos$POS[ii] <= invs$POS_END[jj]) {
            inv.sites[ii,1] <- pos$CHROM[ii]
            inv.sites[ii,2] <- pos$POS[ii]
        }
    }
}
## ----------------------------------------------------------------------------------------------







## CALCULATE ALLELE FREQUENCIES (AFs) ------------------------------------------------------------
## For both ecotypes, pre-selection dune, and post-selection dune and hybrid sources in dune habitat
group.list <- c("post_leafdune_d","post_leafnon_n","pre_dune","post_dune_d","post_hyb_d")  #List groups of interest

AFs <- as.data.frame(matrix(NA, nrow(pos), length(group.list)))            #Empty matrix to hold AFs
numChrs <- as.data.frame(matrix(NA, nrow(pos), length(group.list)))        #Empty matrix to hold num chrms with data for each SNP
colnames(AFs) <- group.list
colnames(numChrs) <- paste(group.list, rep("chrs", length(group.list)), sep="_")

for (gg in 1:length(group.list)) {                                          #Loop over each group
    group <- subset(smpl_gts, smpl_gts$SOURCE_ID==group.list[gg])
    group.gts <- group %>% dplyr::select(-c(FILE,PLT_CODE,POP_ID,SOURCE_ID))
    
    for (ss in 1:(ncol(group.gts))) {                                       #Loop over each SNP
        SNP <- group.gts[,ss]
        AFs[ss,gg] <- sum(SNP, na.rm=TRUE) / (length(which(!is.na(SNP)))*2) #Calculate AF at each SNP
        numChrs[ss,gg] <- length(which(!is.na(SNP)))*2                      #Calculate num chroms w data at each SNP 
    }
}



## For pre-selection hybrids 
## do this separate from groups above so can use random sampling to obtain equal proportions of each hyb type
smpl_gts$POP_ID <- str_replace(smpl_gts$POP_ID, "pre_nondune1", "pre_nondune")
smpl_gts$POP_ID <- str_replace(smpl_gts$POP_ID, "pre_nondune2", "pre_nondune")
smpl_gts$POP_ID <- str_replace(smpl_gts$POP_ID, "pre_nondune3", "pre_nondune")

## Loop over plant types of interest and select random smpl in designated proportions
type.list <- c("pre_BCD","pre_F1D","pre_F1N","pre_nondune")  #Plant types of interest
type.props <- c(20, 25, 25, 10)                              #Designated proportions
num.loci <- nrow(pos)                                        #Number of SNPs
reps <- 1000
AFs.prehyb <- as.data.frame(matrix(NA, reps, num.loci))         
numChrs.prehyb <- as.data.frame(matrix(NA, reps, num.loci))

for (rr in 1:reps) {
    store.gts <- NULL
    
    for (tt in 1:length(type.list)) {
        type <- subset(smpl_gts, POP_ID==type.list[tt])
        type$RAND <- runif(nrow(type), 0, 1)
        type.ord <- type[order(type$RAND),]
        type.rand <- type.ord[1:type.props[tt],]
        store.gts <- as.data.frame(rbind(store.gts, type.rand))
        gts.only <- store.gts %>% dplyr::select(-c(FILE,PLT_CODE,POP_ID,SOURCE_ID,RAND))
    }
    
    for (ss in 1:(ncol(gts.only))) {                                                   #Loop over each SNP
        SNP <- gts.only[,ss]
        AFs.prehyb[rr,ss] <- sum(SNP, na.rm=TRUE) / (length(which(!is.na(SNP)))*2)     #Calculate AF at each SNP
        numChrs.prehyb[rr,ss] <- length(which(!is.na(SNP)))*2                          #Calculate num chrs with data at each SNP
    }
}

#saveRDS(AFs.prehyb, file="AFsPrehyb.rds")
#saveRDS(AFs, file="AFs.rds")


## IF NEEDED, LOAD FILES FROM ABOVE ---------------------------------------------------------------------------------
#AFs.prehyb <- readRDS("AFsPrehyb.rds")
#AFs <- readRDS("AFs.rds")
## ------------------------------------------------------------------------------------------------------------------




## Calculate mean pre-selction hybrid AFs from random sampling 
for (ll in 1:ncol(AFs.prehyb)) {                       #Loop over each SNP
    AFs$pre_hyb_randSamp[ll] <- mean(AFs.prehyb[,ll])
}

## Calculate mean pre-selection hybrid chromosome numbers from random sampling 
for (ll in 1:ncol(numChrs.prehyb)) {                   #Loop over each SNP
    numChrs$pre_hyb_randSamp_chrs[ll] <- mean(numChrs.prehyb[,ll])
}


## Calculate AFs of the reference allele
AFsRef.prehyb <- 1-AFs.prehyb                                      
AFs.ref <- 1-AFs                                      
## ------------------------------------------------------------------------------------------------







## CALCULATE ALLELE FREQ DIFFERENCES AT EACH SNP ----------------------------------------------------
## Ecotype differences
AFdiffs_ecoT.abs <- abs(AFs.ref$post_leafdune_d - AFs.ref$post_leafnon_n)
AFdiffs_ecoT <- AFs.ref$post_leafdune_d - AFs.ref$post_leafnon_n

## Pre- vs post-selection differences
AFdiffs_D_d.abs <- abs(AFs.ref$post_dune_d - AFs.ref$pre_dune)         #Dune samples in dune habitat
AFdiffs_D_d <- AFs.ref$post_dune_d - AFs.ref$pre_dune

AFdiffs_H_d.abs <- abs(AFs.ref$post_hyb_d - AFs.ref$pre_hyb_randSamp)  #Hybrid samples in dune habitat
AFdiffs_H_d <- AFs.ref$post_hyb_d - AFs.ref$pre_hyb_randSamp
## ------------------------------------------------------------------------------------------------







## MODEL POST-SELECTION AF CHANGE AT EACH SNP USING LOGISTIC REGRESSION ----------------------------

## Dune source in dune habitat
num_comb <- 2                                                    #Number of selection types (ie. pre vs post)
num_col <- 3                                                     #Number of columns in input data
pVals.d <- rep(NA, nrow(AFs.ref))                                #Storage vector to hold model outputs

for(i in 1:nrow(AFs.ref)) {
    dat <- as.data.frame(matrix(NA,nrow=num_comb,ncol=num_col))              #Create dataframe for input data
    colnames(dat) <- c("Num.ref","Num.alt","Selection")
    
    dat$Selection <- as.factor(rep(c("pre","post")))
    
    dat[1,1] <- round((AFs.ref$pre_dune[i] * numChrs$pre_dune_chrs[i]))      #Assign allele counts 
    dat[1,2] <- round((AFs$pre_dune[i] * numChrs$pre_dune_chrs)[i])
    
    dat[2,1] <- round((AFs.ref$post_dune_d[i] * numChrs$post_dune_d_chrs)[i])
    dat[2,2] <- round((AFs$post_dune_d[i] * numChrs$post_dune_d_chrs)[i])
    
    glm.d <- glm(cbind(Num.ref, Num.alt) ~ Selection, data=dat, family=binomial(link="logit"))
    
    summary.d <- summary(glm.d)
    
    pVals.d[i] <- as.numeric(summary.d$coefficients[2,4])
}




## Hybrid source in dune habitat 
pVals.h <- rep(NA, nrow(AFs.ref))                                       #Storage vector to hold model outputs                           

for(i in 1:nrow(AFs.ref)) {
    dat <- as.data.frame(matrix(NA,nrow=num_comb,ncol=num_col))         #Create dataframe for input data
    colnames(dat) <- c("Num.ref","Num.alt","Selection")
   
    dat$Selection <- as.factor(rep(c("pre","post")))
    
    dat[1,1] <- round(mean((AFsRef.prehyb[,i] * numChrs.prehyb[,i])))   #Assign allele counts 
    dat[1,2] <- round(mean((AFs.prehyb[,i] * numChrs.prehyb[,i]))) 
    
    dat[2,1] <- round((AFs.ref$post_hyb_d[i] * numChrs$post_hyb_d_chrs)[i])
    dat[2,2] <- round((AFs$post_hyb_d[i] * numChrs$post_hyb_d_chrs)[i])
    
    glm.h <- glm(cbind(Num.ref, Num.alt) ~ Selection, data=dat, family=binomial(link="logit"))
    
    summary.h <- summary(glm.h)
    
    pVals.h[i] <- as.numeric(summary.h$coefficients[2,4])
}
## --------------------------------------------------------------------------------







## MAKE COMBINED DATA FRAME ----------------------------------------------------------------------
comb.dat <- as.data.frame(cbind(pos, pos_on_genome, inv.sites,
                                AFdiffs_D_d, AFdiffs_H_d, AFdiffs_D_d.abs, AFdiffs_H_d.abs,
                                AFdiffs_ecoT, AFdiffs_ecoT.abs,
                                PC1_load, PC1_load.abs, pVals.d, pVals.h))

colnames(comb.dat) <- c("CHROM","POS","POS_ON_GENOME", "INV_CHROM","INV_POS",
                        "DIFFS_Dd","DIFFS_Hd","DIFFS_Dd_ABS","DIFFS_Hd_ABS",
                        "DIFFS_ECOT", "DIFFS_ECOT_ABS", 
                        "PC1_LOAD", "PC1_LOAD_ABS", "PVALS_D", "PVALS_H")


## Subset data to keep only non-inversion SNPs
comb.nonInv <- comb.dat[is.na(comb.dat$INV_CHROM),]


## Subset data to keep only SNPs in inversions
comb.inv <- comb.dat[complete.cases(comb.dat$INV_CHROM), ]  #Keep only non-NA lines
## ------------------------------------------------------------------------------------------------







## LOOK AT RELATIONSHIPS BETWEEN PC1 LOADINGS, ECOTYPE AF DIFFERENCES, & AF CHANGE DURING EXPERIMENT 

## Figure S7
## PC1 loadings vs ecotype AF differences
par(pty="s")
par(mfrow=c(1,1)) 

plot(comb.dat$PC1_LOAD_ABS, comb.dat$DIFFS_ECOT_ABS, pch=16, cex=0.75, family="A", cex.lab=1.75, 
     cex.axis=1.5, xlab="PC1 loadings", ylab="AF difference between ecotypes (dune - nondune)")
lmod.1 <- lm(comb.dat$DIFFS_ECOT_ABS ~ comb.dat$PC1_LOAD_ABS)
summary(lmod.1)
abline(lmod.1, xpd=FALSE)

points(comb.inv$PC1_LOAD_ABS, comb.inv$DIFFS_ECOT_ABS, pch=16, cex=0.75, col="red")
lmod.1.inv <- lm(comb.inv$DIFFS_ECOT_ABS ~ comb.inv$PC1_LOAD_ABS)
summary(lmod.1.inv)
abline(lmod.1.inv, col="red", xpd=FALSE)

lmod.1.nonInv <- lm(comb.nonInv$DIFFS_ECOT_ABS ~ comb.nonInv$PC1_LOAD_ABS)
summary(lmod.1.nonInv)
## ---------------------------------------------------------------------------------





## Figure 4 and S8
## Ecotype allele freq differences vs allele freq change during experiment 

## Hybrid source (figure 4)
pdf('GSD_fig4_221125.pdf', width=8, height=8)
par(mar=c(4,5,3,2))
plot(comb.dat$DIFFS_ECOT, comb.dat$DIFFS_Hd, pch=16, cex=0.5, cex.lab=1.5, cex.main=1.9,
     xlab="AF difference between ecotypes", col="black", font.main=1, asp=1, cex.axis=1.2,
     ylab="AF change between post- & pre-selection", main="Hybrid source in dune habitat")

lmod.h <- lm(comb.dat$DIFFS_Hd ~ comb.dat$DIFFS_ECOT)
summary(lmod.h)
abline(lmod.h, xpd=FALSE, col="black", lwd=1.9)

points(comb.inv$DIFFS_ECOT, comb.inv$DIFFS_Hd, pch=16, cex=0.5, col="grey70")
lmod.h.inv <- lm(comb.inv$DIFFS_Hd ~ comb.inv$DIFFS_ECOT)
summary(lmod.h.inv)
abline(lmod.h.inv, xpd=FALSE, col="grey70", lwd=1.9)

lmod.h.nonInv <- lm(comb.nonInv$DIFFS_Hd ~ comb.nonInv$DIFFS_ECOT)
summary(lmod.h.nonInv)
legend("bottomright", c("Non-inversion loci", "Inversion loci"), col=c("black","grey70"), pch=16,cex=1.4,
       horiz=FALSE, bty="y")
dev.off()


## Calculate correlation coefficients for this relationship
cor.h.all <- cor(comb.dat$DIFFS_ECOT, comb.dat$DIFFS_Hd, use="complete.obs")
cor.h.inv <- cor(comb.inv$DIFFS_ECOT, comb.inv$DIFFS_Hd, use="complete.obs")
cor.h.non <- cor(comb.nonInv$DIFFS_ECOT, comb.nonInv$DIFFS_Hd, use="complete.obs")




## Dune source (figure S8)
plot(comb.dat$DIFFS_ECOT, comb.dat$DIFFS_Dd, pch=16, cex=0.25, cex.lab=1.5, family="A", cex.main=1.5,
     xlab="AF difference between ecotypes (dune - nondune)", col="black", font.main=1, asp=1,
     ylab="AF change between post-\nand pre-selection (post - pre)", main="Dune source in dune habitat")

lmod.d <- lm(comb.dat$DIFFS_Dd ~ comb.dat$DIFFS_ECOT)
summary(lmod.d)
abline(lmod.d, xpd=FALSE, col="black", lwd=1.75)

points(comb.inv$DIFFS_ECOT, comb.inv$DIFFS_Dd, pch=16, cex=0.25, col="red")
lmod.d.inv <- lm(comb.inv$DIFFS_Dd ~ comb.inv$DIFFS_ECOT)
summary(lmod.d.inv)
abline(lmod.d.inv, xpd=FALSE, col="red", lwd=1.75)

lmod.d.nonInv <- lm(comb.nonInv$DIFFS_Dd ~ comb.nonInv$DIFFS_ECOT)
summary(lmod.d.nonInv)


## Calculate correlation coefficients for this relationship
cor.d.all <- cor(comb.dat$DIFFS_ECOT, comb.dat$DIFFS_Dd, use="complete.obs")
cor.d.inv <- cor(comb.inv$DIFFS_ECOT, comb.inv$DIFFS_Dd, use="complete.obs")
cor.d.non <- cor(comb.nonInv$DIFFS_ECOT, comb.nonInv$DIFFS_Dd, use="complete.obs")




## RANDOMIZATION TEST TO ASSESS SIGNIFICANCE OF ECOTYPE DIFF VS AF CHANGE RELATIONSHIP
## Figure S9

## Hybrid source
## Assign indivs to one of 4 groups (wild dune (x27), wild non (x53), pre hyb (x81), post hyb-dune (x21), 
assignmentHyb <- c(rep("post_leafdune_d",27), rep("post_leafnon_n",53), rep("pre_hyb",81),
                   rep("post_hyb_d",21))

## Subset genotype data to have only groups of interest
unique(smpl_gts$SOURCE_ID)
smpl_gtsForRand.hyb <- smpl_gts[smpl_gts$SOURCE_ID!="pre_dune" & !is.na(smpl_gts$SOURCE_ID) &
                                smpl_gts$SOURCE_ID!="post_non_n" & smpl_gts$SOURCE_ID!="pre_non" &
                                smpl_gts$SOURCE_ID!="post_non_d" & smpl_gts$SOURCE_ID!="post_hyb_n" &
                                smpl_gts$SOURCE_ID!="post_dune_n" & smpl_gts$SOURCE_ID!="post_dune_d",]

group.list <- unique(assignmentHyb)               #List groups of interest

rep <- 1000
store.r.hyb <- as.data.frame(matrix(NA, rep, 3))
colnames(store.r.hyb) <- c("allLoci", "invSites", "nonInvSites")

for (ll in 1:rep) {
    smpl_gtsForRand.hyb$SAMPLE <- sample(assignmentHyb)
    
    AFs.rand <- as.data.frame(matrix(NA, nrow(pos), length(group.list)))             #Empty matrix to hold AFs
    colnames(AFs.rand) <- group.list
    
    for (gg in 1:length(group.list)) {                                               #Loop over each group
        group <- subset(smpl_gtsForRand.hyb, smpl_gtsForRand.hyb$SAMPLE==group.list[gg])
        group.gts <- group %>% dplyr::select(-c(FILE,PLT_CODE,POP_ID,SOURCE_ID,SAMPLE))
        
        for (ss in 1:(ncol(group.gts))) {                                            #Loop over each SNP
            SNP <- group.gts[,ss]
            AFs.rand[ss,gg] <- sum(SNP, na.rm=TRUE) / (length(which(!is.na(SNP)))*2) #Calculate AF at each locus
        }
    }
    
    AFsRand.ref <- 1-AFs.rand                                                        #Calculate AFs of the ref allele
    
    ## Calculate AF differences
    AFdiffs_ecoT <- AFsRand.ref$post_leafdune_d - AFsRand.ref$post_leafnon_n         #Ecotype differences
    AFdiffs_H_d <- AFsRand.ref$post_hyb_d - AFsRand.ref$pre_hyb                      #Pre vs post differences

    
    ## Make combined data frame
    comb.datRand <- as.data.frame(cbind(AFdiffs_H_d, AFdiffs_ecoT, inv.sites))
    colnames(comb.datRand) <- c("DIFFS_Hd", "DIFFS_ECOT", "INV_CHROM", "INV_POS")
    
    ## Subset to keep only non-inversion sites
    comb.nonInvRand <- comb.datRand[is.na(comb.datRand$INV_CHROM),]
    
    ## Subset to keep only sites that are in inversions
    comb.invRand <- comb.datRand[complete.cases(comb.datRand$INV_POS), ]  
    
    
    ## Allele freq diffs relationship, store correlation coefficients
    store.r.hyb[ll,1] <- cor(comb.datRand$DIFFS_ECOT, comb.datRand$DIFFS_Hd, use="complete.obs")
    store.r.hyb[ll,2] <- cor(comb.invRand$DIFFS_ECOT, comb.invRand$DIFFS_Hd, use="complete.obs")
    store.r.hyb[ll,3] <- cor(comb.nonInvRand$DIFFS_ECOT, comb.nonInvRand$DIFFS_Hd, use="complete.obs")
}




## Dune source
## Assign indivs to one of 4 groups (wild dune (x27), wild non (x53), pre dune (x56), post dune-dune (x70))
assignmentDune <- c(rep("post_leafdune_d",27), rep("post_leafnon_n",53), rep("pre_dune",56), rep("post_dune_d",70))

## Subset genotype data to have only groups of interest
smpl_gtsForRand.dune <- smpl_gts[smpl_gts$SOURCE_ID!="post_hyb_d" & !is.na(smpl_gts$SOURCE_ID) 
                            & smpl_gts$SOURCE!="post_non_n" & smpl_gts$SOURCE_ID!="pre_hyb"
                            & smpl_gts$SOURCE_ID!="pre_non" & smpl_gts$SOURCE_ID!="post_non_d" 
                            & smpl_gts$SOURCE_ID!="post_hyb_n" & smpl_gts$SOURCE_ID!="post_dune_n",]
                                
group.list <- unique(assignmentDune)            #List groups of interest

rep <- 1000
store.r.dune <- as.data.frame(matrix(NA, rep, 3))
colnames(store.r.dune) <- c("allLoci", "invSites", "nonInvSites")


for (ll in 1:rep) {
    smpl_gtsForRand.dune$SAMPLE <- sample(assignmentDune)                             #Random assignment
    
    ## Calculate AFs
    AFs.rand <- as.data.frame(matrix(NA, ncol(pos), length(group.list)))              #Empty matrix to hold AFs
    colnames(AFs.rand) <- group.list
    
    for (gg in 1:length(group.list)) {                                                #Loop over each group
        group <- subset(smpl_gtsForRand.dune, smpl_gtsForRand.dune$SAMPLE==group.list[gg])
        group.gts <- group %>% dplyr::select(-c(FILE,PLT_CODE,POP_ID,SOURCE_ID,SAMPLE))

        for (ss in 1:(ncol(group.gts))) {                                             #Loop over each SNP
            SNP <- group.gts[,ss]
            AFs.rand[ss,gg] <- sum(SNP, na.rm=TRUE) / (length(which(!is.na(SNP)))*2)  #Calculate AF at each locus
        }
    }
    
    AFsRand.ref <- 1-AFs.rand                                                         #Calculate AFs of the ref allele
    
    ## Calculate AF differences
    AFdiffs_ecoT <- AFsRand.ref$post_leafdune_d - AFsRand.ref$post_leafnon_n          #Ecotype differences
    AFdiffs_D_d <- AFsRand.ref$post_dune_d - AFsRand.ref$pre_dune                     #Pre vs post differences

    
    ## Make combined data frame
    comb.datRand <- as.data.frame(cbind(AFdiffs_D_d, AFdiffs_ecoT, inv.sites))
    colnames(comb.datRand) <- c("DIFFS_Dd", "DIFFS_ECOT", "INV_CHROM", "INV_POS")
    
    ## Subset to keep only non-inversion sites
    comb.nonInvRand <- comb.datRand[is.na(comb.datRand$INV_CHROM),]
    
    ## Subset to keep only sites that are in inversions
    comb.invRand <- comb.datRand[complete.cases(comb.datRand$INV_POS), ] 
    
    
    ## Allele freq diffs relationship,  store correlation coefficients
    store.r.dune[ll,1] <- cor(comb.datRand$DIFFS_ECOT, comb.datRand$DIFFS_Dd, use="complete.obs")
    store.r.dune[ll,2] <- cor(comb.invRand$DIFFS_ECOT, comb.invRand$DIFFS_Dd, use="complete.obs")
    store.r.dune[ll,3] <- cor(comb.nonInvRand$DIFFS_ECOT, comb.nonInvRand$DIFFS_Dd, use="complete.obs")
}




## Plot distributions of correlation coefficeints (r) from randomization tests
par(def_par) 
par(mfrow=c(3,2))  
col_bar <- rgb(t(col2rgb("grey")), alpha=180, maxColorValue=255)
par(mar=c(4.5,5,3.5,2))

hist(store.r.hyb$allLoci, breaks=80, main="All loci- hybrid source", xlab="r", font.main=1,
     cex.lab=1.75, cex.main=2, cex.axis=1.25, border=FALSE, col=col_bar, family="A")
abline(v=quantile(store.r.hyb$allLoci, probs=c(0.025, 0.975)))
abline(v=cor.h.all, col="red", lty=2)

text(-0.15,69, "A", cex=2.5, family="A",xpd=TRUE,font=2)

hist(store.r.dune$allLoci, breaks=90, main="All loci- dune source", xlab="r", family="A",
     cex.axis=1.25, border=FALSE, col=col_bar, font.main=1,cex.lab=1.75,cex.main=2)
abline(v=quantile(store.r.dune$allLoci, probs=c(0.025, 0.975)))
abline(v=cor.d.all, col="red",lty=2)

text(-0.2,59, "B", cex=2.5, family="A",xpd=TRUE,font=2)

hist(store.r.hyb$invSites, breaks=80, main="Inversion loci only- hybrid source", xlab="r",
     cex.lab=1.75, cex.main=2, cex.axis=1.25, border=FALSE, col=col_bar, family="A", font.main=1)
abline(v=quantile(store.r.hyb$invSites, probs=c(0.025, 0.975)))
abline(v=cor.h.inv, col="red", lty=2)

hist(store.r.dune$invSites, breaks=70, main="Inversion loci only- dune source", xlab="r",
     cex.axis=1.25, border=FALSE, col=col_bar, family="A",cex.lab=1.75,cex.main=2,font.main=1,)
abline(v=quantile(store.r.dune$invSites, probs=c(0.025, 0.975)))
abline(v=cor.d.inv, col="red",lty=2)

hist(store.r.hyb$nonInvSites, breaks=90, main="Non-inversion loci only- hybrid source", xlab="r",
     cex.lab=1.75, cex.main=2, cex.axis=1.25, border=FALSE, col=col_bar, family="A", font.main=1)
abline(v=quantile(store.r.hyb$nonInvSites, probs=c(0.025, 0.975)))
abline(v=cor.h.non, col="red", lty=2)

hist(store.r.dune$nonInvSites, breaks=90, main="Non-inversion loci only- dune source", xlab="r",
     cex.axis=1.25, border=FALSE, col=col_bar, family="A",font.main=1,cex.lab=1.75,cex.main=2)
abline(v=quantile(store.r.dune$nonInvSites, probs=c(0.025, 0.975)))
abline(v=cor.d.non, col="red",lty=2)





## Calculate p values (one-sided test) using randomization test distributions 

## Hybrid source
p.allLoci.h <- ((length(store.r.hyb$allLoci[store.r.hyb$allLoci > cor.h.all]) + 1) / 
                     (length(store.r.hyb$allLoci) + 1))

p.inv.h <- ((length(store.r.hyb$invSites[store.r.hyb$invSites > cor.h.inv]) + 1) / 
                 (length(store.r.hyb$invSites) + 1))

p.nonInv.h <- ((length(store.r.hyb$nonInvSites[store.r.hyb$nonInvSites > cor.h.non]) + 1) / 
                   (length(store.r.hyb$nonInvSites) + 1))



## Dune source 
p.allLoci.d <- ((length(store.r.dune$allLoci[store.r.dune$allLoci < cor.d.all]) + 1) / 
                    (length(store.r.dune$allLoci) + 1))

p.inv.d <- ((length(store.r.dune$invSites[store.r.dune$invSites < cor.d.inv]) + 1) / 
                (length(store.r.dune$invSites) + 1))

p.nonInv.d <- ((length(store.r.dune$nonInvSites[store.r.dune$nonInvSites < cor.d.non]) + 1) / 
                   (length(store.r.dune$nonInvSites) + 1))
## ------------------------------------------------------------------------------------------------







## AF change at each SNP across genome

## Figure S10 (post minus pre-selection AF differences)
par(def_par) 
par(mfrow=c(2,1))  
par(mar=c(5,6,5,1))  #Adjust plot window so axis labels don't get cut off


## Dune source
plot(comb.dat$POS_ON_GENOME, comb.dat$DIFFS_Dd_ABS, ylim=c(0,0.65), xlab=NA, pch=16, cex=0.6,
     family="A", xaxt='n', main="Dune source in dune habitat", cex.main=1.5,
     ylab="AF change between post-\nand pre-selection (post - pre)", cex.lab=1.5, font.main=1)

## Highlight SNPs that are in inversions
points(comb.inv$POS_ON_GENOME, comb.inv$DIFFS_Dd_ABS, col="red", pch=16, cex=0.6)

text(-1e7,0.8, "A", cex=1.75, family="A",xpd=TRUE,font=2)

## Add chromosome numbers to x-axis
num.pos <- rep(NA, length(chroms))      #Define position of label

for (zz in 1:length(chroms)) {
    sub.chr <- subset(comb.dat, CHROM==chroms[zz])
    num.pos[zz] <- mean(sub.chr$POS_ON_GENOME) 
}

num.lab <- c(1:17)                      #Define label
axis(side=1, at=num.pos, labels=num.lab, tck=0, family="A", cex=1.5)




## Hybrid source
plot(comb.dat$POS_ON_GENOME, comb.dat$DIFFS_Hd_ABS, ylim=c(0,0.65), 
     xlab="Position on genome (chromosome number)", 
     pch=16, cex=0.6, ylab="AF change between post-\nand pre-selection (post - pre)", 
     cex.lab=1.5, family="A", font.main=1,
     main="Hybrid source in dune habitat", cex.main=1.5, xaxt='n')

## Highlight SNPs that are in inversions
points(comb.inv$POS_ON_GENOME, comb.inv$DIFFS_Hd_ABS, col="red", pch=16, cex=0.6)

text(-1e8,0.85, "B", cex=1.75, family="A",xpd=TRUE,font=2)

## Add chromosome numbers to x-axis
axis(side=1, at=num.pos, labels=num.lab, tck=0, family="A", cex=1.5)






## Figure S11 (p-values from model of allele counts as a function of selection)

## Dune source
plot(comb.dat$POS_ON_GENOME, -log(comb.dat$PVALS_D), ylim=c(0,30), xaxt='n',
     ylab="-log(p)",pch=16,cex.lab=1.5, family="A",cex=0.6, col="black", font.main=1,
     main="Dune source in dune habitat", cex.main=1.5, xlab=NA)

## Highlight SNPs that are in inversions
points(comb.inv$POS_ON_GENOME, -log(comb.inv$PVALS_H), col="red", pch=16, cex=0.6)

## Add chromosome numbers to x-axis
axis(side=1, at=num.pos, labels=num.lab, tck=0, family="A", cex=1.5)

text(-5e8,34, "A", cex=1.75, family="A",xpd=TRUE,font=2)




## Hybrid source
plot(comb.dat$POS_ON_GENOME, -log(comb.dat$PVALS_H), ylim=c(0,30), xaxt='n',
     xlab="Position on genome (chromosome number)", font.main=1,
     ylab="-log(p)",pch=16,cex.lab=1.5, family="A",cex=0.6, col="black",
     main="Hybrid source in dune habitat", cex.main=1.5)#, type="n") 

## Highlight SNPs that are in inversions
points(comb.inv$POS_ON_GENOME, -log(comb.inv$PVALS_H), col="red", pch=16, cex=0.6)

## Add chromosome numbers to x-axis
axis(side=1, at=num.pos, labels=num.lab, tck=0, family="A", cex=1.5)

text(-5e8,34, "B", cex=1.75, family="A",xpd=TRUE,font=2)




## Fig S12
## Relationship between p-values & pre-post selection AF differences
par(pty="s")

## Dune source
plot(comb.dat$DIFFS_Dd_ABS, -log(comb.dat$PVALS_D), cex=0.6, pch=16, ylab="-log(p)",
     xlab="AF change between post- and pre-selection (post - pre)", family='A',
     main="Dune source in dune habitat", font.main=1, cex.lab=1.25, cex.main=1.5)

text(-0.21,19.5, "A", cex=2, family="A",xpd=TRUE,font=2)

summary(lm((-log(comb.dat$PVALS_D) ~ comb.dat$DIFFS_Dd_ABS)))


plot(comb.dat$DIFFS_Hd_ABS, -log(comb.dat$PVALS_H), cex=0.6, pch=16, ylab="-log(p)",
     xlab="AF change between post- and pre-selection (post - pre)", family='A',
     main="Hybrid source in dune habitat", font.main=1, cex.lab=1.25, cex.main=1.5)

text(-0.32,30, "B", cex=2, family="A",xpd=TRUE,font=2)

summary(lm((-log(comb.dat$PVALS_H)) ~ comb.dat$DIFFS_Hd_ABS))
## ------------------------------------------------------------------------------------------------







##  INVERSIONS INDIVIDUALLY -----------------------------------------------------------------------

## Add column with inversion name to inversion data frame
inv.list <- c("pet05.01", "pet07.01", "pet9.comb", "pet11.01", "pet14.01", "pet17.01")

chrom.list <- unique(as.character(comb.inv$CHROM))

for (dd in 1:length(chrom.list)) {
    comb.inv$INV_CODE[grepl(chrom.list[dd], comb.inv$CHROM)] <- inv.list[dd]
}

## Assign labels to the two inversion in chromosome 9 (based on positions in Huang et al 2020)
comb.inv$INV_CODE[comb.inv$INV_CHROM=="Ha412HOChr09" & comb.inv$INV_POS<=140632318] <- "pet09.01"
comb.inv$INV_CODE[comb.inv$INV_CHROM=="Ha412HOChr09" & comb.inv$INV_POS>=171481816] <- "pet09.02"





## Plot distributions of AF change  
par(def_par)
par(mar=c(6,5,3,3)) 

## Dune source
layout(matrix(c(1,9,2,6,3,7,4,8,5,10), 2, 5, byrow=FALSE))

## All non-inversion loci
hist(comb.nonInv$DIFFS_Dd_ABS, breaks=70, col="black",border=FALSE,
     family="A",cex.lab=1.75,main="All non-inversion loci- Dune samples", cex.axis=1.25,
     xlab=NA, font.main=1, cex.main=1.75)


## Each inversion separately
break.list <- c(10,20,30,20,30,30,20)
inv.codes <- unique(comb.inv$INV_CODE)

for (hh in 1:length(inv.codes)) {
    subsetDat <- subset(comb.inv, comb.inv$INV_CODE==inv.codes[hh])
    hist(subsetDat$DIFFS_Dd_ABS, breaks=break.list[hh], col="grey40", border=FALSE, 
         family="A", cex.lab=1.75, main=inv.codes[hh],
         cex.axis=1.25,xlab=NA, xlim=c(0,0.4), font.main=1, cex.main=1.75)
}


## All inversion loci
hist(comb.inv$DIFFS_Dd_ABS, breaks=40, col="grey40", border=FALSE, family="A",
     cex.lab=1.75, main="All inversion loci", cex.main=1.75, font.main=1,
     cex.axis=1.25,xlab=NA, xlim=c(0,0.4))

mtext("AF change (absolute value)", side=1, outer=TRUE, cex=1.5, family="A", at=0.5,-2)




## Rank test of AF change distribution differences
rankTest.d <- as.data.frame(matrix(NA, nrow=length(inv.codes), ncol=3))
colnames(rankTest.d) <- c("Inversion", "z_statistic", "p_value")
rankTest.d$Inversion <- inv.codes

for (rr in 1:length(inv.codes)) {
    subsetDat <- subset(comb.inv, comb.inv$INV_CODE==inv.codes[rr])
    rankTest.d$z_statistic[rr] <- twoSampleLinearRankTest(comb.nonInv$DIFFS_Dd_ABS, subsetDat$DIFFS_Dd_ABS, 
                                shift.type = "location")$statistic
    rankTest.d$p_value[rr] <- twoSampleLinearRankTest(comb.nonInv$DIFFS_Dd_ABS, subsetDat$DIFFS_Dd_ABS, 
                            shift.type = "location")$p.value
}




## Hybrid source
layout(matrix(c(1,9,2,6,3,7,4,8,5,10), 2, 5, byrow=FALSE))

## All non-inversion loci
hist(comb.nonInv$DIFFS_Hd_ABS, breaks=70, col="black",border=FALSE,
     family="A",cex.lab=1.75,main="All non-inversion loci- Hybrid samples", cex.axis=1.25,
     xlab=NA, font.main=1, cex.main=1.75)


## Each inversion separately
break.list <- c(30,30,40,20,30,20,20)
inv.codes <- unique(comb.inv$INV_CODE)

for (hh in 1:length(inv.codes)) {
    subsetDat <- subset(comb.inv, comb.inv$INV_CODE==inv.codes[hh])
    hist(subsetDat$DIFFS_Hd_ABS, breaks=break.list[hh], col="grey40", border=FALSE, 
         family="A", cex.lab=1.75, main=inv.codes[hh],
         cex.axis=1.25, xlab=NA, font.main=1, cex.main=1.75)
}


## All inversion loci
hist(comb.inv$DIFFS_Hd_ABS, breaks=40, col="grey40", border=FALSE, family="A",
     cex.lab=1.75, main="All inversion loci", font.main=1, cex.main=1.75,
     cex.axis=1.25,xlab=NA, xlim=c(0,0.6))

mtext("AF change (absolute value)", side=1, outer=TRUE, cex=1.5, family="A", at=0.5,-2)




## Rank test of AF change distribution differences
rankTest.h <- as.data.frame(matrix(NA, nrow=length(inv.codes), ncol=3))
colnames(rankTest.h) <- c("Inversion", "z_statistic", "p_value")
rankTest.h$Inversion <- inv.codes

for (rr in 1:length(inv.codes)) {
    subsetDat <- subset(comb.inv, comb.inv$INV_CODE==inv.codes[rr])
    rankTest.h$z_statistic[rr] <- twoSampleLinearRankTest(comb.nonInv$DIFFS_Hd_ABS, subsetDat$DIFFS_Hd_ABS, 
                                  shift.type = "location")$statistic
    rankTest.h$p_value[rr] <- twoSampleLinearRankTest(comb.nonInv$DIFFS_Hd_ABS, subsetDat$DIFFS_Hd_ABS, 
                              shift.type = "location")$p.value
}
## ------------------------------------------------------------------------------------------------







## CALCULATE AF CHANGE IN EACH INVERSION WHEN INVERSIONS TREATED AS SINGLE LOCI --------------------
## Load VCF tables for each inversion and run PCAdapt
filename.5 <- read.pcadapt(path_to_file.5, type="vcf")
filename.7 <- read.pcadapt(path_to_file.7, type="vcf")
filename.9.1 <- read.pcadapt(path_to_file.9.1, type="vcf")
filename.9.2 <- read.pcadapt(path_to_file.9.2, type="vcf")
filename.11 <- read.pcadapt(path_to_file.11, type="vcf")
filename.14 <- read.pcadapt(path_to_file.14, type="vcf")
filename.17 <- read.pcadapt(path_to_file.17, type="vcf")

x.5 <- pcadapt(filename.5,K=5)
x.7 <- pcadapt(filename.7,K=5)
x.9.1 <- pcadapt(filename.9.1,K=5)
x.9.2 <- pcadapt(filename.9.2,K=5)
x.11 <- pcadapt(filename.11,K=5)
x.14 <- pcadapt(filename.14,K=5)
x.17 <- pcadapt(filename.17,K=5)




## Visualize clustering of samples
smpls$ECOT[grepl("*dune*", smpls$SOURCE_ID)] <- "Dune"    #Add ecotype label
smpls$ECOT[grepl("*hyb*", smpls$SOURCE_ID)] <- "Hybrid"
smpls$ECOT[grepl("*non*", smpls$SOURCE_ID)] <- "Nondune"

plot(x.5,option="scores", pop=smpls$ECOT)             
plot(x.7,option="scores", pop=smpls$ECOT)             
plot(x.9.1,option="scores", pop=smpls$ECOT)                      
plot(x.9.2,option="scores", pop=smpls$ECOT)                   
plot(x.11,option="scores", pop=smpls$ECOT)           
plot(x.14,option="scores", pop=smpls$ECOT)                   
plot(x.17,option="scores", pop=smpls$ECOT)            



## Use kmeans function to cluster individuals into 3 clusters suggesting genotype 
## (homozygote ref, heterozygote, homozygote alt) based on PC1 score.
## kmeans will calculate cluster means and assignment

## Run kmeans multiple times to assess consistency of clusters
store.kmeans <- as.data.frame(matrix(NA, nrow=1000, ncol=21))

colnames(store.kmeans) <- paste(rep(inv.codes,1,each=3), rep(c(1,2,3), length(inv.codes)), sep="_")
                                
for (kk in 1:1000) {
    kmean.5 <- kmeans(x.5$scores[,1], 3)
    kmean.7 <- kmeans(x.7$scores[,1], 3)
    kmean.9.1 <- kmeans(x.9.1$scores[,1], 3)
    kmean.9.2 <- kmeans(x.9.2$scores[,1], 3)
    kmean.11 <- kmeans(x.11$scores[,1], 3)
    kmean.14 <- kmeans(x.14$scores[,1], 3)
    kmean.17 <- kmeans(x.17$scores[,1], 3)
    
    store.kmeans[kk,1:3] <- kmean.5$size
    store.kmeans[kk,4:6] <- kmean.7$size
    store.kmeans[kk,7:9] <- kmean.9.1$size
    store.kmeans[kk,10:12] <- kmean.9.2$size
    store.kmeans[kk,13:15] <- kmean.11$size
    store.kmeans[kk,16:18] <- kmean.14$size
    store.kmeans[kk,19:21] <- kmean.17$size
}

rbind(sort(unique(store.kmeans$pet05.01_1)), sort(unique(store.kmeans$pet05.01_2)), sort(unique(store.kmeans$pet05.01_3)))
rbind(sort(unique(store.kmeans$pet07.01_1)), sort(unique(store.kmeans$pet07.01_2)), sort(unique(store.kmeans$pet07.01_3)))
rbind(sort(unique(store.kmeans$pet09.01_1)), sort(unique(store.kmeans$pet09.01_2)), sort(unique(store.kmeans$pet09.01_3)))
rbind(sort(unique(store.kmeans$pet09.02_1)), sort(unique(store.kmeans$pet09.02_2)), sort(unique(store.kmeans$pet09.02_3)))
rbind(sort(unique(store.kmeans$pet11.01_1)), sort(unique(store.kmeans$pet11.01_2)), sort(unique(store.kmeans$pet11.01_3)))
rbind(sort(unique(store.kmeans$pet14.01_1)), sort(unique(store.kmeans$pet14.01_2)), sort(unique(store.kmeans$pet14.01_3)))
rbind(sort(unique(store.kmeans$pet17.01_1)), sort(unique(store.kmeans$pet17.01_2)), sort(unique(store.kmeans$pet17.01_3)))

## All inversions except pet09.02 and pet14.01 have consistent cluster assignment.
## Based on visualization of pet14.01, it looks like only 2 invidiuals are likely homozygous for the dune allele
i <- 0
while (i == 0 ) {
    kmean.14 <- kmeans(x.14$scores[,1], 3)
    if (min(kmean.14$size)==2) {
        kmean.14 <- kmean.14
        i <- 1
    } else { i <- 0 }
}

## Based on visualization of pet09.02, separation of samples along PC1 is not distinct and therefore
## clustering based on PC1 to assign genotypes may not be appropriate. 





## Combine columns with cluster assignment and PC1 scores for each inversion
smpls <- smpls %>% mutate(CLUSTER.5=kmean.5$cluster, CLUSTER.7=kmean.7$cluster, CLUSTER.9.1=kmean.9.1$cluster,
                          CLUSTER.9.2=kmean.9.2$cluster, CLUSTER.11=kmean.11$cluster, CLUSTER.14=kmean.14$cluster,
                          CLUSTER.17=kmean.17$cluster,
                          SCORES.5=x.5$scores[,1], SCORES.7=x.7$scores[,1], SCORES.9.1=x.9.1$scores[,1], 
                          SCORES.9.2=x.9.2$scores[,1], SCORES.11=x.11$scores[,1], SCORES.14=x.14$scores[,1],
                          SCORES.17=x.17$scores[,1])


## Determine which cluster is most dune- vs nondune-like, and assign 'nondune' cluster as genotype 0 (homozygote reference)
## Inv pet5
if (mean(smpls$SCORES.5[smpls$ECOT=="Nondune"], na.rm=TRUE) > mean(smpls$SCORES.5[smpls$ECOT=="Dune"], na.rm=TRUE)) {
    smpls$GT.5[grepl(which.max(kmean.5$centers), smpls$CLUSTER.5)] <- 0
    smpls$GT.5[grepl(which.min(kmean.5$centers), smpls$CLUSTER.5)] <- 2
    smpls$GT.5[is.na(smpls$GT.5)] <- 1
} else {
    smpls$GT.5[grepl(which.min(kmean.5$centers), smpls$CLUSTER.5)] <- 0
    smpls$GT.5[grepl(which.max(kmean.5$centers), smpls$CLUSTER.5)] <- 2
    smpls$GT.5[is.na(smpls$GT.5)] <- 1
}

## Inv pet7
if (mean(smpls$SCORES.7[smpls$ECOT=="Nondune"], na.rm=TRUE) > mean(smpls$SCORES.7[smpls$ECOT=="Dune"], na.rm=TRUE)) {
    smpls$GT.7[grepl(which.max(kmean.7$centers), smpls$CLUSTER.7)] <- 0
    smpls$GT.7[grepl(which.min(kmean.7$centers), smpls$CLUSTER.7)] <- 2
    smpls$GT.7[is.na(smpls$GT.7)] <- 1
} else {
    smpls$GT.7[grepl(which.min(kmean.7$centers), smpls$CLUSTER.7)] <- 0
    smpls$GT.7[grepl(which.max(kmean.7$centers), smpls$CLUSTER.7)] <- 2
    smpls$GT.7[is.na(smpls$GT.7)] <- 1
}

## Inv pet9.1
if (mean(smpls$SCORES.9.1[smpls$ECOT=="Nondune"], na.rm=TRUE) > mean(smpls$SCORES.9.1[smpls$ECOT=="Dune"], na.rm=TRUE)) {
    smpls$GT.9.1[grepl(which.max(kmean.9.1$centers), smpls$CLUSTER.9.1)] <- 0
    smpls$GT.9.1[grepl(which.min(kmean.9.1$centers), smpls$CLUSTER.9.1)] <- 2
    smpls$GT.9.1[is.na(smpls$GT.9.1)] <- 1
} else {
    smpls$GT.9.1[grepl(which.min(kmean.9.1$centers), smpls$CLUSTER.9.1)] <- 0
    smpls$GT.9.1[grepl(which.max(kmean.9.1$centers), smpls$CLUSTER.9.1)] <- 2
    smpls$GT.9.1[is.na(smpls$GT.9.1)] <- 1
}

## Inv pet9.2
if (mean(smpls$SCORES.9.2[smpls$ECOT=="Nondune"], na.rm=TRUE) > mean(smpls$SCORES.9.2[smpls$ECOT=="Dune"], na.rm=TRUE)) {
    smpls$GT.9.2[grepl(which.max(kmean.9.2$centers), smpls$CLUSTER.9.2)] <- 0
    smpls$GT.9.2[grepl(which.min(kmean.9.2$centers), smpls$CLUSTER.9.2)] <- 2
    smpls$GT.9.2[is.na(smpls$GT.9.2)] <- 1
} else {
    smpls$GT.9.2[grepl(which.min(kmean.9.2$centers), smpls$CLUSTER.9.2)] <- 0
    smpls$GT.9.2[grepl(which.max(kmean.9.2$centers), smpls$CLUSTER.9.2)] <- 2
    smpls$GT.9.2[is.na(smpls$GT.9.2)] <- 1
}

## Inv pet11
if (mean(smpls$SCORES.11[smpls$ECOT=="Nondune"], na.rm=TRUE) > mean(smpls$SCORES.11[smpls$ECOT=="Dune"], na.rm=TRUE)) {
    smpls$GT.11[grepl(which.max(kmean.11$centers), smpls$CLUSTER.11)] <- 0
    smpls$GT.11[grepl(which.min(kmean.11$centers), smpls$CLUSTER.11)] <- 2
    smpls$GT.11[is.na(smpls$GT.11)] <- 1
} else {
    smpls$GT.11[grepl(which.min(kmean.11$centers), smpls$CLUSTER.11)] <- 0
    smpls$GT.11[grepl(which.max(kmean.11$centers), smpls$CLUSTER.11)] <- 2
    smpls$GT.11[is.na(smpls$GT.11)] <- 1
}

## Inv pet14
if (mean(smpls$SCORES.14[smpls$ECOT=="Nondune"], na.rm=TRUE) > mean(smpls$SCORES.14[smpls$ECOT=="Dune"], na.rm=TRUE)) {
    smpls$GT.14[grepl(which.max(kmean.14$centers), smpls$CLUSTER.14)] <- 0
    smpls$GT.14[grepl(which.min(kmean.14$centers), smpls$CLUSTER.14)] <- 2
    smpls$GT.14[is.na(smpls$GT.14)] <- 1
} else {
    smpls$GT.14[grepl(which.min(kmean.14$centers), smpls$CLUSTER.14)] <- 0
    smpls$GT.14[grepl(which.max(kmean.14$centers), smpls$CLUSTER.14)] <- 2
    smpls$GT.14[is.na(smpls$GT.14)] <- 1
}

## Inv pet17
if (mean(smpls$SCORES.17[smpls$ECOT=="Nondune"], na.rm=TRUE) > mean(smpls$SCORES.17[smpls$ECOT=="Dune"], na.rm=TRUE)) {
    smpls$GT.17[grepl(which.max(kmean.17$centers), smpls$CLUSTER.17)] <- 0
    smpls$GT.17[grepl(which.min(kmean.17$centers), smpls$CLUSTER.17)] <- 2
    smpls$GT.17[is.na(smpls$GT.17)] <- 1
} else {
    smpls$GT.17[grepl(which.min(kmean.17$centers), smpls$CLUSTER.17)] <- 0
    smpls$GT.17[grepl(which.max(kmean.17$centers), smpls$CLUSTER.17)] <- 2
    smpls$GT.17[is.na(smpls$GT.17)] <- 1
}





## Calculate allele frequences for each inversion locus
## Subset genotypes by groups of interest
preN <- subset(smpls, SOURCE_ID=="pre_non") %>% dplyr::select(starts_with("GT."))
preD <- subset(smpls, SOURCE_ID=="pre_dune") %>% dplyr::select(starts_with("GT."))
postD.d <- subset(smpls, SOURCE_ID=="post_dune_d") %>% dplyr::select(starts_with("GT."))
postH.d <- subset(smpls, SOURCE_ID=="post_hyb_d") %>% dplyr::select(starts_with("GT."))

num.invs <- length(inv.codes)
cols <- c("AFpreN", "AFpreD", "AFpreH", "AFpostD", "AFpostH","AFchngD", "AFchngH")
AFs.inv <- as.data.frame(matrix(NA, num.invs, length(cols)))
colnames(AFs.inv) <- cols

for (gg in (1:num.invs)) {
    AFs.inv$AFpreN[gg] <- sum(preN[,gg])/(nrow(preN)*2)                #Calculate allele freq for each inversion locus
    AFs.inv$AFpreD[gg] <- sum(preD[,gg])/(nrow(preD)*2)
    AFs.inv$AFpostD[gg] <- sum(postD.d[,gg])/(nrow(postD.d)*2)
    AFs.inv$AFpostH[gg] <- sum(postH.d[,gg])/(nrow(postH.d)*2)
}

AFs.inv <- AFs.inv %>% mutate(AFpreD_chrs=nrow(preD)*2, AFpostD_chrs=nrow(postD.d)*2,
                              AFpostH_chrs=nrow(postH.d)*2)



## For pre-selection hybrids, do random sampling
## Loop over plt types of interest and select random smpl in designated proportions

## Combine different non-dune populations into single group
smpls$POP_ID[grepl("pre_nondune1", smpls$POP_ID)] <- "pre_nondune" 
smpls$POP_ID[grepl("pre_nondune2", smpls$POP_ID)] <- "pre_nondune" 
smpls$POP_ID[grepl("pre_nondune3", smpls$POP_ID)] <- "pre_nondune" 

reps <- 1000
AFsInv.prehyb <- as.data.frame(matrix(NA, reps, num.invs))
colnames(AFsInv.prehyb) <- inv.codes

for (rr in 1:reps) {
    store.gts <- NULL
    
    for (tt in 1:length(type.list)) {
        type <- subset(smpls, POP_ID==type.list[tt])
        type$RAND <- runif(nrow(type), 0, 1)
        type.ord <- type[order(type$RAND),]
        type.rand <- type.ord[1:type.props[tt],]
        store.gts <- as.data.frame(rbind(store.gts, type.rand))
    }
    
    AFsInv.prehyb[rr,1] <- sum(store.gts$GT.5) / (nrow(store.gts)*2)
    AFsInv.prehyb[rr,2] <- sum(store.gts$GT.7) / (nrow(store.gts)*2)
    AFsInv.prehyb[rr,3] <- sum(store.gts$GT.9.1) / (nrow(store.gts)*2)
    AFsInv.prehyb[rr,4] <- sum(store.gts$GT.9.2) / (nrow(store.gts)*2)
    AFsInv.prehyb[rr,5] <- sum(store.gts$GT.11) / (nrow(store.gts)*2)
    AFsInv.prehyb[rr,6] <- sum(store.gts$GT.14) / (nrow(store.gts)*2)
    AFsInv.prehyb[rr,7] <- sum(store.gts$GT.17) / (nrow(store.gts)*2)
}


## Calculate mean AFs from random sampling 
AFs.inv$AFpreH <- c(mean(AFsInv.prehyb$pet05.01), mean(AFsInv.prehyb$pet07.01), 
                     mean(AFsInv.prehyb$pet09.01), mean(AFsInv.prehyb$pet09.02), 
                     mean(AFsInv.prehyb$pet11.01), mean(AFsInv.prehyb$pet14.01), 
                     mean(AFsInv.prehyb$pet17.01))

## Add in chrom numbers for pre-hyb
AFs.inv$AFpreH_chrs <- nrow(store.gts)*2


## Calculate AF change from pre to post
AFs.inv$AFchngD <- (AFs.inv$AFpostD - AFs.inv$AFpreD)
AFs.inv$AFchngH <- (AFs.inv$AFpostH - AFs.inv$AFpreH)






## FISHER'S EXACT TEST TO TEST SIGNIFICANCE OF AF CHANGE IN INVERSIONS 

## Dune source
allele1.preDune <- round(AFs.inv$AFpreD * AFs.inv$AFpreD_chrs)
allele2.preDune <- round((1-AFs.inv$AFpreD) * AFs.inv$AFpreD_chrs)
allele1.postDune <- round(AFs.inv$AFpostD * AFs.inv$AFpostD_chrs)
allele2.postDune <- round((1-AFs.inv$AFpostD) * AFs.inv$AFpostD_chrs)

for (i in 1:num.invs) {
    mat <- matrix(c(allele1.preDune[i],allele2.preDune[i],allele1.postDune[i],allele2.postDune[i]), nrow=2)
    print(mat)
    print(fisher.test(mat)[1])
}


## Hybrid source
allele1.preHyb <- round(AFs.inv$AFpreH * AFs.inv$AFpreH_chrs)
allele2.preHyb <- round((1-AFs.inv$AFpreH) * AFs.inv$AFpreH_chrs)
allele1.postHyb <- round(AFs.inv$AFpostH * AFs.inv$AFpostH_chrs)
allele2.postHyb <- round((1-AFs.inv$AFpostH) * AFs.inv$AFpostH_chrs)

for (i in 1:num.invs) {
    mat <- matrix(c(allele1.preHyb[i],allele2.preHyb[i],allele1.postHyb[i],allele2.postHyb[i]), nrow=2)
    print(mat)
    print(fisher.test(mat)[1])
}





## COMPARE INVERSION AF CHANGE TO DISTRIBUTION OF AF CHANGE IN NON-INVERSION LOCI
percnts.d <- rep(NA, nrow(AFs.inv))

for (pp in 1:nrow(AFs.inv)) {
    percnts.d[pp] <- round((length(comb.nonInv$DIFFS_Dd_ABS[comb.nonInv$DIFFS_Dd_ABS < abs(AFs.inv$AFchngD[pp])]) / 
                                length(comb.nonInv$DIFFS_Dd_ABS)) * 100)
}


percnts.h <- rep(NA, nrow(AFs.inv))

for (pp in 1:nrow(AFs.inv)) {
    percnts.h[pp] <- round((length(comb.nonInv$DIFFS_Dd_ABS[comb.nonInv$DIFFS_Dd_ABS < abs(AFs.inv$AFchngH[pp])]) / 
                                length(comb.nonInv$DIFFS_Dd_ABS)) * 100)
}
## ------------------------------------------------------------------------------------------------



