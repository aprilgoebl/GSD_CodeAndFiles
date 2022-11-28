## Great Sand Dune (GSD), Helianthus petiolaris, reciprocal transplant data analysis script
## Script that uses 2012 (& 2010, 2011) GSD reciprocal transplant data in demographic models
## Estimates fitness components, lambda and performs life table response experiment (LTRE): 
## 1. Uses logistic regression to model emergence and survival rates 
## 2. Uses negative binomial to model fecundity (number of seeds or inflorescences produced)
## 3. Estimates lambdas (population growth rate)
## 4. Plots modeled fitness components and lambda estimates 
## 5. Performs and plots results for LTRE




rm(list=ls())
dev.off()
windowsFonts(A = windowsFont("Times New Roman"))




## LOAD PACKAGES ----------------------------------------------------------------------------
library(plotrix)
library(lme4)
library(car)
library(MuMIn)
library(popbio)
library(MASS)
library(dplyr)
library(ggeffects)
## -----------------------------------------------------------------------------------------------------



## Note: All data obtained from 10.5061/dryad.223p4 
## READ IN AND FORMAT 2012 DATA ------------------------------------------------------------------------

## Emergence and survival
datEarly <- read.table(file ='2012_RT_Plots.csv', sep=',', header = TRUE)  
datEarly$num.planted <- rep(10, nrow(datEarly))                #Add column with total seeds planted per plot
datEarly <- subset(datEarly, datEarly$ecotype != "control")    #Remove controls
datEarly$pop[grepl("1300", datEarly$type)] <- "D1"             #Add columns with population names, group hybrid types 
datEarly$pop[grepl("1701", datEarly$type)] <- "D2"
datEarly$pop[grepl("1547", datEarly$type)] <- "D3"
datEarly$pop[grepl("2001", datEarly$type)] <- "N1"
datEarly$pop[grepl("1791", datEarly$type)] <- "N2"
datEarly$pop[grepl("1363", datEarly$type)] <- "N3"
datEarly$pop[grepl("F1D", datEarly$type)] <- "Hybrid"
datEarly$pop[grepl("F1N", datEarly$type)] <- "Hybrid"
datEarly$pop[grepl("BCD", datEarly$type)] <- "Hybrid"
datEarly$pop[grepl("BCN", datEarly$type)] <- "Hybrid"
datEarly$plot <- as.factor(datEarly$plot)                      #Change data types for plot and pop
datEarly$pop <- as.factor(datEarly$pop)                        

datEarly <- datEarly[datEarly$suspected_volunteer != 'Y',]     #Remove suspected volunteers (natural recruits)
## Note: Plants labelled as survivors in their foreign habitat that produced seeds outside of their ecotypic 
## weight range and inside the range of the local ecotype (95% CIs) were designated as volunteeers (natural recruits)



## Fecundity (number of seeds produced)
datLate <- read.table(file ='2012_RT_Plants_full.csv', sep=',', header = TRUE)
datLate <- subset(datLate, datLate$ecotype != "control")      #Remove controls
datLate$pop[grepl("1300", datLate$type)] <- "D1"              #Add columns w population names, group hybrid types 
datLate$pop[grepl("1701", datLate$type)] <- "D2"
datLate$pop[grepl("1547", datLate$type)] <- "D3"
datLate$pop[grepl("2001", datLate$type)] <- "N1"
datLate$pop[grepl("1791", datLate$type)] <- "N2"
datLate$pop[grepl("1363", datLate$type)] <- "N3"
datLate$pop[grepl("F1D", datLate$type)] <- "Hybrid"
datLate$pop[grepl("F1N", datLate$type)] <- "Hybrid"
datLate$pop[grepl("BCD", datLate$type)] <- "Hybrid"
datLate$pop[grepl("BCN", datLate$type)] <- "Hybrid"
datLate$plot <- as.factor(datLate$plot)                       #Change data types for plot and pop
datLate$pop <- as.factor(datLate$pop)                         

datLate <- datLate[datLate$suspected_volunteer != 'Y',]       #Remove suspected volunteers  

## When zero heads collected, but plant had inflorescences, change zero seed counts to NA  
datLate$n.seeds.counted[datLate$n.seeds.counted==0 & datLate$heads.counted==0] <- NA

#Adjust number of seeds based on number of inflorescences (sometimes not all heads were collected)
datLate$numSdsAdj <- rep(NA, nrow(datLate))
for (ss in 1:nrow(datLate)) {                                 
  if (datLate$heads.counted[ss] < datLate$n.flowers[ss] && datLate$heads.counted[ss]>0) {        
    datLate$numSdsAdj[ss] <- (datLate$n.seeds.counted[ss] / datLate$heads.counted[ss]) * datLate$n.flowers[ss]
  } else {
    datLate$numSdsAdj[ss] <- datLate$n.seeds.counted[ss]
  }  
}

## Adjust for proportion of seeds eaten by larva (pre-dispersal predation)
datLate$prop.seeds.eaten[is.na(datLate$prop.seeds.eaten)] <- 0
datLate$numSdsAdj <- datLate$numSdsAdj - (datLate$numSdsAdj*datLate$prop.seeds.eaten) 
## ------------------------------------------------------------------------------------------------------



## READ IN AND FORMAT 2011 DATA ------------------------------------------------------------------------
dat2011 <- read.table(file ='2011_RT.csv', sep=',', header = TRUE)  
str(dat2011)
dat2011$ecotype <- as.character(dat2011$ecotype)
dat2011$ecotype[dat2011$ecotype=="non-dune"] <- "non.dune"

dat2011$plot <- rep(1:50, each=6) #Rename plot numbers
dat2011$plot <- as.factor(dat2011$plot) 
dat2011$ecotype <- as.factor(dat2011$ecotype)
dat2011$num.planted <- rep(15, nrow(dat2011))              #Add column with seed planted per subplot
dat2011 <- dat2011[dat2011$ecotype!="control",]            #Remove control rows 

dat2011.early <- dat2011 %>% group_by(plot, pop, ecotype) %>% summarise(SUM=sum(num.planted), 
                                                              COUNT=(n()), EMRG=sum(emg)) 
## -----------------------------------------------------------------------------------------------------



## READ IN AND FORMAT 2010 DATA -----------------------------------------------------------------------
dat2010 <- read.table(file ='2010_RT_long.csv', sep=',', header = TRUE)  
str(dat2010)
dat2010$ecotype <- as.character(dat2010$ecotype)
dat2010$habitat <- as.character(dat2010$habitat) 
dat2010$ecotype[dat2010$ecotype=="intermediate"] <- "hybrid"
dat2010$ecotype[dat2010$ecotype=="non-dune"] <- "non.dune"
dat2010$habitat[dat2010$habitat=="non"] <- "non.dune"

dat2010$plot <- as.factor(dat2010$plot) 
dat2010$pop <- as.factor(dat2010$pop) 
dat2010$ecotype <- as.factor(dat2010$ecotype)
dat2010$habitat <- as.factor(dat2010$habitat) 
dat2010$num.planted <- rep(1, nrow(dat2010))                 #Add column with seed planted per indiv
dat2010 <- dat2010[dat2010$treatment=="4degrees",]           #Remove scarified treatment 

dat2010.early <- dat2010 %>% group_by(habitat, plot, pop, ecotype) %>% summarise(SUM=sum(num.planted), 
                                                           COUNT=(n()), EMRG=sum(emg), SURV=sum(surv)) 

## Subset late (reproduction) data
dat2010.late <- dat2010[dat2010$surv==1,]    
## ------------------------------------------------------------------------------------------------------





## MODEL EMERGENCE, SURVIVAL, AND FECUNDITY DATA --------------------------------------------------------
## Emergence 
## Logistic regression 
## 2012
fitEmrg.mm <- glmer(cbind(emrg, num.planted-emrg) ~ ecotype * habitat + (1|pop), 
                    data=datEarly, family=binomial(link="logit"))

## 2011
fitEmrg.mm2011 <- glmer(cbind(EMRG, SUM-EMRG) ~ ecotype + (1|pop), data=dat2011.early,
                       family=binomial(link="logit"))

## 2010
fitEmrg.2010 <- glm(cbind(EMRG, SUM-EMRG) ~ ecotype * habitat, 
                        data=dat2010.early, family=binomial(link="logit"))




## Seedling-to-adult survival 
## Logistic regression 
## 2012
fitSurv.mm <- glmer(cbind(surv, emrg-surv) ~ ecotype * habitat + (1|pop), 
                    data=datEarly, family=binomial(link="logit"))
#isSingular(fitSurv.mm, tol=1e-05)

## No survival data from 2011

## 2010
fitSurv.2010 <- glm(cbind(SURV, EMRG-SURV) ~ ecotype * habitat, 
                        data=dat2010.early, family=binomial(link="logit"))




## Number of seeds (or inflorescences) produced (fecundity)
## Negative binomial regression
fitSds.mm <- glmer.nb(round(numSdsAdj) ~ ecotype * habitat + (1|pop), data=datLate)

## No fecundity data from 2011

## 2010
fitSds.2010 <- glm.nb(n.flower ~ ecotype * habitat, data=dat2010.late)
## --------------------------------------------------------------------------------------------------------





## Obtain model statistics and summaries -----------------------------------------------------------------
smry.E.mm <- summary(fitEmrg.mm)
smry.S.mm <- summary(fitSurv.mm)
smry.F.mm <- summary(fitSds.mm)

smry.E.mm2011 <- summary(fitEmrg.mm2011)

Anova(fitEmrg.mm)
Anova(fitSurv.mm)
Anova(fitSds.mm)

Anova(fitEmrg.2010)
Anova(fitSurv.2010)
Anova(fitSds.2010)

Anova(fitEmrg.mm2011)
## --------------------------------------------------------------------------------------------------------



## ESTIMATE FITNESS COMPONENTS, CIs, AND LAMBDA FOR EACH SOURCE IN EACH HABITAT

## Using ggeffects package to get estimates and CIs
## 2012
prEmrg <- ggpredict(fitEmrg.mm, c("ecotype", "habitat"))
prSurv <- ggpredict(fitSurv.mm, c("ecotype", "habitat"))
prSds <- ggpredict(fitSds.mm, c("ecotype", "habitat"))

## Calculate CIs as 1.96xSE for seed estimates (ggpredict can't compute CIs for this model)
## assume normal approximation of negative binomial given large mean estimated values
## Calculate variance using the estimated mean and each raw data point
prSds
calcVar <- function(x,m,n) {(sum((x-m)^2, na.rm=TRUE))/(n-1)}
dd <- subset(datLate, habitat=="dune" & ecotype=="dune")
dh <- subset(datLate, habitat=="dune" & ecotype=="hybrid")
dn <- subset(datLate, habitat=="dune" & ecotype=="non.dune")
nd <- subset(datLate, habitat=="non.dune" & ecotype=="dune")
nh <- subset(datLate, habitat=="non.dune" & ecotype=="hybrid")
nn <- subset(datLate, habitat=="non.dune" & ecotype=="non.dune")
ddCI <- ((sqrt(calcVar(dd$numSdsAdj, prSds$predicted[1], nrow(dd))))/(nrow(dd)))*1.96
dhCI <- ((sqrt(calcVar(dh$numSdsAdj, prSds$predicted[3], nrow(dh))))/(nrow(dh)))*1.96
dnCI <- ((sqrt(calcVar(dn$numSdsAdj, prSds$predicted[5], nrow(dn))))/(nrow(dn)))*1.96
ndCI <- ((sqrt(calcVar(nd$numSdsAdj, prSds$predicted[2], nrow(nd))))/(nrow(nd)))*1.96
nhCI <- ((sqrt(calcVar(nh$numSdsAdj, prSds$predicted[4], nrow(nh))))/(nrow(nh)))*1.96
nnCI <- ((sqrt(calcVar(nn$numSdsAdj, prSds$predicted[6], nrow(nn))))/(nrow(nn)))*1.96


## Estimate lambda
lam <- prEmrg$predicted * prSurv$predicted * prSds$predicted
## Adjust lambda to include estimates of seed survival
sdSurv <- 0.3       #Set values for seed survival rate (eg. predation, fungal attack, etc.)
lam.mod <- lam * sdSurv


## 2010
prEmrg.2010 <- ggpredict(fitEmrg.2010, c("ecotype", "habitat"))
prSurv.2010 <- ggpredict(fitSurv.2010, c("ecotype", "habitat"))
prSds.2010 <- ggpredict(fitSds.2010, c("ecotype", "habitat"))


##2011
prEmrg.mm2011 <- ggpredict(fitEmrg.mm2011, ("ecotype"))






## PLOT RESULTS AS REACTION NORMS --------------------------------------------------------------------------------
## Assign colours
colD <- rgb(t(col2rgb("orange")),alpha=215,maxColorValue = 255)
colN <- rgb(t(col2rgb("royalblue")),alpha=215,maxColorValue = 255)
colH <- rgb(t(col2rgb("#33a02c")),alpha=215,maxColorValue = 255)
cols <- c(colD,colH,colN,colD,colH,colN)

colD_ln <- rgb(t(col2rgb("orange")),alpha=100,maxColorValue = 255)
colH_ln <- rgb(t(col2rgb("#33a02c")),alpha=100,maxColorValue = 255)
colN_ln <- rgb(t(col2rgb("royalblue")),alpha=100,maxColorValue = 255)


svg('GSD_fig1.svg', width=17, height=12)
layout(rbind(c(1,1,2,2,3,3), c(5,4,4,4,4,6), c(7,4,4,4,4,8)))
par(mar=c(5,6,3,1))




## Plot fitness component estimates with confidence intervals

##2012 mixed model

## Emergence 
## Reaction norm 
## ggeffects results
plot(NA,NA,ylim=c(0,0.2),xlim=c(1,2),frame.plot=TRUE,ylab="Emergence rate",xaxt='n',
     xlab=NA,cex.lab=2.5,cex.axis=2.25, family="A")
axis(side=1, at=c(1.05), labels = "Dune", cex.axis=2.5, tick=FALSE, family="A")
axis(side=1, at=c(1.92), labels = "Non-dune", cex.axis=2.5, tick=FALSE, family="A")
points(c(prEmrg$predicted[1],prEmrg$predicted[2]),col=colD,lwd=8)
points(c(prEmrg$predicted[3],prEmrg$predicted[4]),col=colH,lwd=8)
points(c(prEmrg$predicted[5],prEmrg$predicted[6]), col=colN,lwd=8) 

lines(c(prEmrg$predicted[1],prEmrg$predicted[2]),col=colD_ln,lwd=2)
lines(c(prEmrg$predicted[3],prEmrg$predicted[4]),col=colH_ln,lwd=2)
lines(c(prEmrg$predicted[5],prEmrg$predicted[6]), col=colN_ln,lwd=2) 

## Add CIs
arrows(1,prEmrg$predicted[1],1, (prEmrg$conf.low[1]), lwd = 1.5, angle = 90, code = 3, length=0, col=colD)
arrows(1,prEmrg$predicted[1],1, (prEmrg$conf.high[1]), lwd = 1.5, angle = 90, code = 3, length=0, col=colD)
arrows(2,prEmrg$predicted[2],2, (prEmrg$conf.low[2]), lwd = 1.5, angle = 90, code = 3, length=0,col=colD)
arrows(2,prEmrg$predicted[2],2, (prEmrg$conf.high[2]), lwd = 1.5, angle = 90, code = 3, length=0,col=colD)
arrows(1,prEmrg$predicted[3],1, (prEmrg$conf.low[3]), lwd = 1.5, angle = 90, code = 3, length=0,col=colH)
arrows(1,prEmrg$predicted[3],1, (prEmrg$conf.high[3]), lwd = 1.5, angle = 90, code = 3, length=0,col=colH)
arrows(2,prEmrg$predicted[4],2, (prEmrg$conf.low[4]), lwd = 1.5, angle = 90, code = 3, length=0,col=colH)
arrows(2,prEmrg$predicted[4],2, (prEmrg$conf.high[4]), lwd = 1.5, angle = 90,code = 3, length=0,col=colH)
arrows(1,prEmrg$predicted[5],1, (prEmrg$conf.low[5]), lwd = 1.5, angle = 90,code = 3, length=0,col=colN)
arrows(1,prEmrg$predicted[5],1, (prEmrg$conf.high[5]), lwd = 1.5, angle = 90,code = 3, length=0,col=colN)
arrows(2,prEmrg$predicted[6],2, (prEmrg$conf.low[6]), lwd = 1.5, angle = 90,code = 3, length=0,col=colN)
arrows(2,prEmrg$predicted[6],2, (prEmrg$conf.high[6]), lwd = 1.5, angle = 90,code = 3, length=0,col=colN)

text(0.825,0.21, "A", cex=3, family="A",xpd=TRUE,font=2)




## Survival 
## Reaction norm
## ggeffects results
plot(NA,NA,ylim=c(0,0.6),xlim=c(1,2),frame.plot=TRUE,ylab="Survival rate",xaxt='n',
     xlab=NA,cex.lab=2.5,cex.axis=2.25, family="A")
axis(side=1, at=c(1.05), labels = "Dune", cex.axis=2.5, tick=FALSE, family="A")
axis(side=1, at=c(1.92), labels = "Non-dune", cex.axis=2.5, tick=FALSE, family="A")
points(c(prSurv$predicted[1],prSurv$predicted[2]),col=colD,lwd=8)
points(c(prSurv$predicted[3],prSurv$predicted[4]),col=colH,lwd=8)
points(c(prSurv$predicted[5],prSurv$predicted[6]), col=colN,lwd=8) 

lines(c(prSurv$predicted[1],prSurv$predicted[2]),col=colD_ln,lwd=2)
lines(c(prSurv$predicted[3],prSurv$predicted[4]),col=colH_ln,lwd=2)
lines(c(prSurv$predicted[5],prSurv$predicted[6]), col=colN_ln,lwd=2) 

## Add CIs (jittered for visual clarity)
arrows(1.001,prSurv$predicted[1],1.001, (prSurv$conf.low[1]), lwd = 1.5, angle = 90, code = 3, length=0, col=colD)
arrows(1.001,prSurv$predicted[1],1.001, (prSurv$conf.high[1]), lwd = 1.5, angle = 90, code = 3, length=0, col=colD)
arrows(2.001,prSurv$predicted[2],2.001, (prSurv$conf.low[2]), lwd = 1.5, angle = 90, code = 3, length=0,col=colD)
arrows(2.001,prSurv$predicted[2],2.001, (prSurv$conf.high[2]), lwd = 1.5, angle = 90, code = 3, length=0,col=colD)
arrows(0.999,prSurv$predicted[3],0.999, (prSurv$conf.low[3]), lwd = 1.5, angle = 90, code = 3, length=0,col=colH)
arrows(0.999,prSurv$predicted[3],0.999, (prSurv$conf.high[3]), lwd = 1.5, angle = 90, code = 3, length=0,col=colH)
arrows(2,prSurv$predicted[4],2, (prSurv$conf.low[4]), lwd = 1.5, angle = 90, code = 3, length=0,col=colH)
arrows(2,prSurv$predicted[4],2, (prSurv$conf.high[4]), lwd = 1.5, angle = 90,code = 3, length=0,col=colH)
arrows(1,prSurv$predicted[5],1, (prSurv$conf.low[5]), lwd = 1.5, angle = 90,code = 3, length=0,col=colN)
arrows(1,prSurv$predicted[5],1, (prSurv$conf.high[5]), lwd = 1.5, angle = 90,code = 3, length=0,col=colN)
arrows(1.999,prSurv$predicted[6],1.999, (prSurv$conf.low[6]), lwd = 1.5, angle = 90,code = 3, length=0,col=colN)
arrows(1.999,prSurv$predicted[6],1.999, (prSurv$conf.high[6]), lwd = 1.5, angle = 90,code = 3, length=0,col=colN)

text(0.825,0.63, "B", cex=3, family="A",xpd=TRUE,font=2)



## Seed output 
## Reaction norm 
## ggeffects results
plot(NA,NA,ylim=c(1,5500),xlim=c(1,2),frame.plot=TRUE,ylab="Total seed number",xaxt='n',
     xlab=NA,cex.lab=2.5,cex.axis=2.25, family="A", log='y')
axis(side=1, at=c(1.05), labels = "Dune", cex.axis=2.5, tick=FALSE, family="A")
axis(side=1, at=c(1.92), labels = "Non-dune", cex.axis=2.5, tick=FALSE, family="A")
points(c(prSds$predicted[1],prSds$predicted[2]),col=colD,lwd=8)
points(c(prSds$predicted[3],prSds$predicted[4]),col=colH,lwd=8)
points(c(prSds$predicted[5],prSds$predicted[6]), col=colN,lwd=8) 

lines(c(prSds$predicted[1],prSds$predicted[2]),col=colD_ln,lwd=2)
lines(c(prSds$predicted[3],prSds$predicted[4]),col=colH_ln,lwd=2)
lines(c(prSds$predicted[5],prSds$predicted[6]), col=colN_ln,lwd=2) 

## Add CIs (jittered for visual clarity)
arrows(1.001,prSds$predicted[1],1.001, prSds$predicted[1]-ddCI, lwd = 1.5, angle = 90, code = 3, length=0, col=colD)
arrows(1.001,prSds$predicted[1],1.001, prSds$predicted[1]+ddCI, lwd = 1.5, angle = 90, code = 3, length=0, col=colD)
arrows(2.001,prSds$predicted[2],2.001, prSds$predicted[2]-ndCI, lwd = 1.5, angle = 90, code = 3, length=0,col=colD)
arrows(2.001,prSds$predicted[2],2.001, prSds$predicted[2]+ndCI, lwd = 1.5, angle = 90, code = 3, length=0,col=colD)
arrows(0.999,prSds$predicted[3],0.999, prSds$predicted[3]-dhCI, lwd = 1.5, angle = 90, code = 3, length=0,col=colH)
arrows(0.999,prSds$predicted[3],0.999, prSds$predicted[3]+dhCI, lwd = 1.5, angle = 90, code = 3, length=0,col=colH)
arrows(2,prSds$predicted[4],2, prSds$predicted[4]-nhCI, lwd = 1.5, angle = 90, code = 3, length=0,col=colH)
arrows(2,prSds$predicted[4],2, prSds$predicted[4]+nhCI, lwd = 1.5, angle = 90,code = 3, length=0,col=colH)
arrows(1,prSds$predicted[5],1, 1, lwd = 1.5, angle = 90,code = 3, length=0,col=colN) #Plot as 0 (log(1)) since can't plot negative value on log scale
arrows(1,prSds$predicted[5],1, prSds$predicted[5]+dnCI, lwd = 1.5, angle = 90,code = 3, length=0,col=colN)
arrows(1.999,prSds$predicted[6],1.999, prSds$predicted[6]-nnCI, lwd = 1.5, angle = 90,code = 3, length=0,col=colN)
arrows(1.999,prSds$predicted[6],1.999, prSds$predicted[6]+nnCI, lwd = 1.5, angle = 90,code = 3, length=0,col=colN)

text(0.825,8999, "C", cex=3, family="A",xpd=TRUE,font=2)



## Lambda 
## Reaction norm 
## ggeffects results
plot(NA,NA,ylim=c(0,4),xlim=c(1,2),frame.plot=TRUE,ylab="Lambda",xaxt='n',
     xlab="Reciprocal Transplant Habitat",cex.lab=2.7,cex.axis=2.5, family="A")#, log='y')
axis(side=1, at=c(1.05),labels = "Dune",cex.axis=2.7,tick=FALSE, family="A")
axis(side=1, at=c(1.92),labels = "Non-dune",cex.axis=2.75,tick=FALSE, family="A")
points(c(lam.mod[1],lam.mod[4]),col=colD,lwd=10)
points(c(lam.mod[2],lam.mod[5]),col=colH,lwd=10)
points(c(lam.mod[3],lam.mod[6]), col=colN,lwd=10) 
abline(h=1, lty=2, col="grey20")

lines(c(lam.mod[1],lam.mod[4]),col=colD_ln,lwd=3)
lines(c(lam.mod[2],lam.mod[5]),col=colH_ln,lwd=3)
lines(c(lam.mod[3],lam.mod[6]), col=colN_ln,lwd=3)

text(1, 3.5, "*", cex=3, family="A")
text(2, 2.15, "N.S.", cex=2.2, family="A")

text(0.9,4.2, "D", cex=3, family="A",xpd=TRUE,font=2)

legend("topright", c("Dune source", "Hybrid source", "Non-dune source"), col=cols,cex=2.5,
       inset=c(-0.28,0.0009), xpd=NA, horiz=FALSE,bty="y",pch=16, text.font=6)
dev.off()
## --------------------------------------------------------------------------------------------------





## 2012 simple model - barplots -------------------------------------------------------

## Re-order for plotting
index<-c(1,4,2,5,3,6)
prEmrg.ord <- prEmrg$predicted[order(index)]
prEmrgCIhi.ord <- prEmrg$conf.high[order(index)]
prEmrgCIlo.ord <- prEmrg$conf.low[order(index)]
prSurv.ord <- prSurv$predicted[order(index)]
prSurvCIhi.ord <- prSurv$conf.high[order(index)]
prSurvCIlo.ord <- prSurv$conf.low[order(index)]
prSds.ord <- prSds$predicted[order(index)]
prSdsCI.ord <- c(ddCI,dhCI,dnCI,ndCI,nhCI,nnCI)


layout(rbind(c(1,1,2,2,3,3)))
par(mar=c(4,7,3,1)+.1)  #Adjust plot window so axis labels don't get cut off
#Margins:(bottom,left,top,right)

## Emergence barplot
plotCI(barplot(prEmrg.ord, beside=T, ylab="Emergence rate",xlab=NA,border=FALSE, family="A",
               cex.axis=1.5, cex.lab=2,col=cols,ylim=c(0,0.2)),prEmrg.ord,uiw=(prEmrgCIhi.ord-prEmrg.ord),
       liw=(prEmrg.ord-prEmrgCIlo.ord), add=TRUE,pch=FALSE,sfrac=0) 
axis(side=1, at=c(1.8),labels = "Dune",cex.axis=1.5,tick=FALSE, family="A")
axis(side=1, at=c(5.5),labels = "Non-dune",cex.axis=1.5,tick=FALSE, family="A")
abline(v=3.7,lty=2,lwd=1, col="grey70")    

## Survival barplot
plotCI(barplot(prSurv.ord, beside=T, ylab="Survival rate",xlab=NA,border=FALSE, family="A",
               cex.axis=1.5, cex.lab=2,col=cols,ylim=c(0,0.8)),prSurv.ord,uiw=(prSurvCIhi.ord-prSurv.ord),
       liw=(prSurv.ord-prSurvCIlo.ord), add=TRUE,pch=FALSE,sfrac=0) 
axis(side=1, at=c(1.8),labels = "Dune",cex.axis=1.5,tick=FALSE, family="A")
axis(side=1, at=c(5.5),labels = "Non-dune",cex.axis=1.5,tick=FALSE, family="A")
abline(v=3.7,lty=2,lwd=1, col="grey70")    

## Fecundity barplot
plotCI(barplot(prSds.ord, beside=T, ylab="Total seed number",xlab=NA,border=FALSE, family="A",
               cex.axis=1.5, cex.lab=2,col=cols,ylim=c(0,7000)),prSds.ord,uiw=(prSdsCI.ord),
       liw=(prSdsCI.ord), add=TRUE,pch=FALSE,sfrac=0) 
axis(side=1, at=c(1.8),labels = "Dune",cex.axis=1.5,tick=FALSE, family="A")
axis(side=1, at=c(5.5),labels = "Non-dune",cex.axis=1.5,tick=FALSE, family="A")
abline(v=3.7,lty=2,lwd=1, col="grey70")    

legend("topleft", legend=c("Dune source", "Hybrid source", "Non-dune source"), col=c(colD,colH,colN), 
       pch=16,cex=1.5,text.font=6)
## -------------------------------------------------------------------------------------------------




## 2011 simple model - barplots --------------------------------------------------------------------

## Emergence barplot
plotCI(barplot(prEmrg.mm2011$predicted, beside=T, ylab="Emergence rate",xlab=NA,border=FALSE, family="A",
               cex.axis=1.5, cex.lab=2,col=c(cols[1],cols[3]),ylim=c(0,0.08)),prEmrg.mm2011$predicted,uiw=(prEmrg.mm2011$conf.high-prEmrg.mm2011$predicted),
       liw=(prEmrg.mm2011$predicted-prEmrg.mm2011$conf.low), add=TRUE,pch=FALSE,sfrac=0)
axis(side=1, at=c(1.3),labels = "Dune-Habitat",cex.axis=1.5,tick=FALSE, family="A")

legend("topright", legend=c("Dune source", "Non-dune source"), col=c(colD,colN), 
       pch=16,cex=1.75,text.font=6)
## ------------------------------------------------------------------------------------------------




## 2010 simple model - barplots -------------------------------------------------------------------

## Re-order for plotting
index<-c(1,4,2,5,3,6)
prEmrg.2010ord <- prEmrg.2010$predicted[order(index)]
prEmrgCIhi.2010ord <- prEmrg.2010$conf.high[order(index)]
prEmrgCIlo.2010ord <- prEmrg.2010$conf.low[order(index)]
prSurv.2010ord <- prSurv.2010$predicted[order(index)]
prSurvCIhi.2010ord <- prSurv.2010$conf.high[order(index)]
prSurvCIlo.2010ord <- prSurv.2010$conf.low[order(index)]
prSds.2010ord <- prSds.2010$predicted[order(index)]
prSdsCIhi.2010ord <- prSds.2010$conf.high[order(index)]
prSdsCIlo.2010ord <- prSds.2010$conf.low[order(index)]

## Change estimates and CIs for categories with no data to zero
prEmrgCIhi.2010ord[3] <- 0
prSurvCIhi.2010ord[2] <- 0
prSurvCIlo.2010ord[2] <- 0
prSurv.2010ord[3] <-0
prSurvCIhi.2010ord[3] <- 0
prSurvCIlo.2010ord[3] <- 0
prSds.2010ord[1:3] <-0
prSdsCIhi.2010ord[1:3] <- 0
prSdsCIlo.2010ord[1:3] <- 0

par(mar=c(4,7,3,1)+.1)  #Adjust plot window so axis labels don't get cut off

## Emergence barplot
plotCI(barplot(prEmrg.2010ord, beside=T, ylab="Emergence rate",xlab=NA,border=FALSE, family="A",
     cex.axis=1.5, cex.lab=2,col=cols,ylim=c(0,0.2)),prEmrg.2010ord,uiw=(prEmrgCIhi.2010ord-prEmrg.2010ord),
       liw=(prEmrg.2010ord-prEmrgCIlo.2010ord), add=TRUE,pch=FALSE,sfrac=0) 
axis(side=1, at=c(1.8),labels = "Dune",cex.axis=1.5,tick=FALSE, family="A")
axis(side=1, at=c(5.5),labels = "Non-dune",cex.axis=1.5,tick=FALSE, family="A")
abline(v=3.7,lty=2,lwd=1, col="grey70")    

plotCI(barplot(prSurv.2010ord, beside=T, ylab="Survival rate",xlab=NA,border=FALSE, family="A",
       cex.axis=1.5, cex.lab=2,col=cols,ylim=c(0,1)),prSurv.2010ord,uiw=(prSurvCIhi.2010ord-prSurv.2010ord),
       liw=(prSurv.2010ord-prSurvCIlo.2010ord), add=TRUE,pch=FALSE,sfrac=0) 
axis(side=1, at=c(1.8),labels = "Dune",cex.axis=1.5,tick=FALSE, family="A")
axis(side=1, at=c(5.5),labels = "Non-dune",cex.axis=1.5,tick=FALSE, family="A")
abline(v=3.7,lty=2,lwd=1, col="grey70")    

plotCI(barplot(prSds.2010ord, beside=T, ylab="Total inflorescence number",xlab=NA,border=FALSE, family="A",
       cex.axis=1.5, cex.lab=2,col=cols,ylim=c(0,14)),prSds.2010ord,uiw=(prSdsCIhi.2010ord-prSds.2010ord),
       liw=(prSds.2010ord-prSdsCIlo.2010ord), add=TRUE,pch=FALSE,sfrac=0) 
axis(side=1, at=c(1.8),labels = "Dune",cex.axis=1.5,tick=FALSE, family="A")
axis(side=1, at=c(5.5),labels = "Non-dune",cex.axis=1.5,tick=FALSE, family="A")
abline(v=3.7,lty=2,lwd=1, col="grey70")  

legend("topleft", legend=c("Dune source", "Hybrid source", "Non-dune source"), col=c(colD,colH,colN), 
       pch=16,cex=1.5,text.font=6)
## ------------------------------------






## LOOK AT HOW DIFFERENT SEED SURVIVAL ESTIMATES AFFECT LAMBDA --------------------------------------
#Order: dune-in-dune, dune-in-non, hyb-in-dune, hyb-in-non, non-in-dune, non-in-non

sdSurv.test <- 0.59       #Try different values for seed survival rate
lam * sdSurv.test
lam.mod

#16% - 58% seed survival will results in ONLY local ecotype having lambda GREATER than 1. 
#Below that, nondune in nondunes has lambda less than 1
#Above that, hybs in dunes has lambda above 1
## --------------------------------------------------------------------------------------------------



## LOOK AT HOW SENSITIVE LAMBDA IS TO EMERGENCE VALUES ---------------------------------------------
pert <- 0.5

## Make storage vectors for all combinations for fitness components & lambda
prEmrg.test <- rep(NA, nrow(prEmrg))

## Calc fitness components, & corresponding lambda for each source/habitat combination
prEmrg.test[1] <- 1/(1+exp(-(smry.E.mm$coefficients[1]))) * pert
prEmrg.test[2] <- 1/(1+exp(-(smry.E.mm$coefficients[1] + smry.E.mm$coefficients[2]))) * pert
prEmrg.test[3] <- 1/(1+exp(-(smry.E.mm$coefficients[1] + smry.E.mm$coefficients[3]))) * pert
prEmrg.test[4] <- 1/(1+exp(-(smry.E.mm$coefficients[1] + smry.E.mm$coefficients[4]))) * pert
prEmrg.test[5] <- 1/(1+exp(-(smry.E.mm$coefficients[1] + smry.E.mm$coefficients[2] + smry.E.mm$coefficients[4]
                            + smry.E.mm$coefficients[5]))) * pert
prEmrg.test[6] <- 1/(1+exp(-(smry.E.mm$coefficients[1] + smry.E.mm$coefficients[3] + smry.E.mm$coefficients[4]
                            + smry.E.mm$coefficients[6]))) * pert

lam.test <- prEmrg.test * prSurv.ord * prSds.ord

## Adjust lambda to include estimates of seed survival
sdSurv <- 0.3       #Set value for seed survival rate (eg. predation, fungal attack, etc.)
lam.testMod <- lam.test * sdSurv
lam.modOrd <- lam.mod[order(index)]
lam.modOrd
lam.testMod
## --------------------------------------------------------------------------------------------------



 

## PARAMETRIC BOOTSTRAP TO ESTIMATE ERROR ON LAMBDA ESTIMATES ---------------------------------------
## Obtain parameter values and covariance matrix for each fitness component
coefsEmrg <- smry.E.mm$coefficients    
covEmrg <- vcov(fitEmrg.mm)

coefsSurv <- smry.S.mm$coefficients    
covSurv <- vcov(fitSurv.mm)

coefsSds <- smry.F.mm$coefficients
covSds <- vcov(fitSds.mm)



## Loop & calculate different values for fitness components & corresponding lambdas 
## from parameters selected from multivariate normal distribtuion
reps <- 10000                                  #Number of repetitions 
matLam <- matrix(NA,reps,nrow(prEmrg))         #Matrices for lambda and fitness component storage 
matEmrg <- matrix(NA,reps,nrow(prEmrg))
matSurv <- matrix(NA,reps,nrow(prEmrg))
matSds <- matrix(NA,reps,nrow(prEmrg))


i <- 1
while (i <= reps) {                                                         #Use while loop so can include conditional if estimated seeds are negative
  
  ## Emergence 
  rand_coefsEmrg <- mvrnorm(n=1, coefsEmrg[,1], covEmrg)                    #Select new param values using mean & covar of params, & mvn distribution
  fitEmrg.mm@beta <- rand_coefsEmrg                                         #Put selected parameters into model fit
  prEmrg.rand <- ggpredict(fitEmrg.mm, c("ecotype", "habitat"))             #Calc predicted emeregence values for all source & habitat combos
  vecEmrg <- prEmrg.rand$predicted
  
  ## Survival
  rand_coefsSurv <- mvrnorm(n=1, coefsSurv[,1], covSurv)       
  fitSurv.mm@beta <- rand_coefsSurv                    
  prSurv.rand <- ggpredict(fitSurv.mm, c("ecotype", "habitat"))
  vecSurv <- prSurv.rand$predicted
  
  ## Seeds produced   
  rand_coefsSds <- mvrnorm(n=1, coefsSds[,1], covSds)                      
  fitSds.mm@beta <- rand_coefsSds
  vecSds <- rep(NA,nrow(coefsSds))
  prSds.rand <- ggpredict(fitSds.mm, c("ecotype", "habitat"))
  vecSds <- prSds.rand$predicted
  
  
  if (all(vecSds > 0)) {                                    #Conditional: proceed only if all seed estimates are positive
    matEmrg[i,] <- vecEmrg                                  #Store emergence values
    matSurv[i,] <- vecSurv                                  #Store survival values
    matSds[i,] <- vecSds                                    #Store seed number values
    
    vecLam <- vecEmrg * vecSurv * vecSds                    #Calculate and store lambda
    matLam[i,] <- vecLam
    
    i <- i + 1                                              #Increment counter 
  }
  
} 
## End loop

## Save output
saveRDS(matLam, file="matLam.rds")
saveRDS(matEmrg, file="matEmrg.rds")
saveRDS(matSurv, file="matSurv.rds")
saveRDS(matSds, file="matSds.rds")



## Multiply lambda by seed survival constant
matLam <- matLam * sdSurv




## CALCULATE DIFFERENTIAL LAMBDA & 95% QUANTILES TO IDENTIFY SIGNIFICANT COMPARISONS 
## Order of source x habitat combos
#dune-in-dunes, dune-in-non.dunes, hyb-in-dunes, hyb-in-non.dunes, non.dune-in-dunes, non.dune-in-non.dunes 

## In dune habitat
DvH <- quantile(log(matLam[,1]) - log(matLam[,3]), probs=c(0.025, 0.975))
DvN <- quantile(log(matLam[,1]) - log(matLam[,5]), probs=c(0.025, 0.975))
HvN <- quantile(log(matLam[,3]) - log(matLam[,5]), probs=c(0.025, 0.975))

## In nondune habitat
NvH.n <- quantile(log(matLam[,6]) - log(matLam[,4]), probs=c(0.025, 0.975))
NvD.n <- quantile(log(matLam[,6]) - log(matLam[,2]), probs=c(0.025, 0.975))
HvD.n <- quantile(log(matLam[,4]) - log(matLam[,2]), probs=c(0.025, 0.975))

## If zero is less than lower quantile or greater than upper quantile, comparison is 'significant' (based on 95% quantiles)




## Plot distributions of log lambda ratios
par(mfrow=c(3,2))  #Plot in 6 panels

## In dune habitat
hist(log(matLam[,1]) - log(matLam[,3]), breaks=40, xlim=c(-6,10),xlab=NA,main="Dune vs Hybrid in dune habitat",
     cex.lab=2,cex.main=2,col=rgb(t(col2rgb("grey")),alpha=180,maxColorValue=255), border=FALSE, family="A",font.main=1)  
abline(v=DvH[1])
abline(v=DvH[2])
abline(v=0, lty=2, col="red")

hist(log(matLam[,1]) - log(matLam[,5]), breaks=80, xlim=c(-6,10),xlab=NA,main="Dune vs Nondune in dune habitat",
     cex.lab=2,cex.main=2,col=rgb(t(col2rgb("grey")),alpha=180,maxColorValue=255), border=FALSE, family="A",font.main=1) 
abline(v=DvN[1])
abline(v=DvN[2])
abline(v=0, lty=2, col="red")


hist(log(matLam[,3]) - log(matLam[,5]), breaks=90, xlim=c(-6,10),xlab="Log lambda ratio",main="Hybrid vs Nondune in dune habitat",
     cex.lab=2,cex.main=2,col=rgb(t(col2rgb("grey")),alpha=180,maxColorValue=255), border=FALSE, family="A",font.main=1) 
abline(v=HvN[1])
abline(v=HvN[2])
abline(v=0, lty=2, col="red")


## In non-dune habitat
hist(log(matLam[,6]) - log(matLam[,4]), breaks=80, xlim=c(-6,10),xlab=NA,main="Nondune vs Hybrid in nondune habitat",
     cex.lab=2,cex.main=2,col=rgb(t(col2rgb("grey")),alpha=180,maxColorValue=255), border=FALSE, family="A",font.main=1)  
abline(v=NvH.n[1])
abline(v=NvH.n[2])
abline(v=0, lty=2, col="red")


hist(log(matLam[,6]) - log(matLam[,2]), breaks=70, xlim=c(-6,10),xlab=NA,main="Nondune vs Dune in nondune habitat",
     cex.lab=2,cex.main=2,col=rgb(t(col2rgb("grey")),alpha=180,maxColorValue=255), border=FALSE, family="A",font.main=1) 
abline(v=NvD.n[1])
abline(v=NvD.n[2])
abline(v=0, lty=2, col="red")


hist(log(matLam[,4]) - log(matLam[,2]), breaks=90, xlim=c(-6,10),xlab="Log lambda ratio",main="Hybrid vs Dune in nondune habitat",
     cex.lab=2,cex.main=2,col=rgb(t(col2rgb("grey")),alpha=180,maxColorValue=255), border=FALSE, family="A",font.main=1) 
abline(v=c(HvD.n[1],HvD.n[2]))
abline(v=0, lty=2, col="red")
## ------------------------------------------------------------------------------------------------------------------







## LTRE -------------------------------------------------------------------------------------------------------------
## Order of source x habitat combos
#dune-in-dunes, dune-in-non.dunes, hyb-in-dunes, hyb-in-non.dunes, non.dune-in-dunes, non.dune-in-non.dunes 

## FOR DUNE AND NONDUNE SOURCE IN DUNE HABITAT ----------------------------------------------------------------------
## Calculate mean FC value across d and nd sources
mean.emrg.d <- mean(c(prEmrg$predicted[1],prEmrg$predicted[5])) #Average of dune and non.dune in dune habitat
mean.surv.d <- mean(c(prSurv$predicted[1],prSurv$predicted[5]))
mean.fec.d <- mean(c(prSds$predicted[1],prSds$predicted[5]))
mean.sdSurv <- 0.3


## Calculate sensitivities
vr.d <- list(emrg=mean.emrg.d, surv=mean.surv.d, fec=mean.fec.d, sdSurv=mean.sdSurv )
lambda.d <- expression(emrg * surv * fec * sdSurv)
x.d <- vitalsens(lambda.d, vr.d)



## Calculate change in fitness components
## local FC - foreign FC 
delta.emrg.d <- prEmrg$predicted[1] - prEmrg$predicted[5]
delta.surv.d <- prSurv$predicted[1] - prSurv$predicted[5] 
delta.fec.d <- prSds$predicted[1] - prSds$predicted[5]


## LTRE 
ltr.emrg.d <- delta.emrg.d * x.d$sensitivity[1]
ltr.surv.d <- delta.surv.d * x.d$sensitivity[2]
ltr.fec.d <- delta.fec.d * x.d$sensitivity[3]
ltr.sdSurv.d <- 0

diff.lam.d <- lam.mod[1] - lam.mod[5]
ltr.adj <- c(ltr.emrg.d, ltr.surv.d, ltr.fec.d, ltr.sdSurv.d, diff.lam.d)



## Calculate 95% quantiles ---------------------------------------------------------
## Loop over FC matrices and re-calc sensitivity and delta FC
store_ltr <- matrix(NA, nrow=nrow(matEmrg), ncol=3)

for (cc in 1:nrow(matEmrg)) {
  
  ## Calc mean FC
  mean.emrg.d <- mean(c(matEmrg[cc,1],matEmrg[cc,5]))
  mean.surv.d <- mean(c(matSurv[cc,1],matSurv[cc,5]))
  mean.fec.d <- mean(c(matSds[cc,1],matSds[cc,5]))
  
  ## Calculate sensitivities
  vr.d <- list(emrg=mean.emrg.d, surv=mean.surv.d, fec=mean.fec.d, sdSurv=mean.sdSurv )
  lambda.d <- expression(emrg * surv * fec * sdSurv)
  x.d <- vitalsens(lambda.d, vr.d)  
  
  ## Calc delta FC
  delta.emrg.d <- matEmrg[cc,1] - matEmrg[cc,5]
  delta.surv.d <- matSurv[cc,1] - matSurv[cc,5] 
  delta.fec.d <- matSds[cc,1] - matSds[cc,5]
  
  
  ## LTRE (product of sensitivity and change in FC)
  ltr.emrg.d <- delta.emrg.d * x.d$sensitivity[1]
  ltr.surv.d <- delta.surv.d * x.d$sensitivity[2]
  ltr.fec.d <- delta.fec.d * x.d$sensitivity[3]
  
  store_ltr[cc,1] <- ltr.emrg.d 
  store_ltr[cc,2] <- ltr.surv.d 
  store_ltr[cc,3] <- ltr.fec.d 
}


## Calculate 95% quantiles
emrg_CI <- quantile(store_ltr[,1], probs=c(0.025, 0.975))
surv_CI <- quantile(store_ltr[,2], probs=c(0.025, 0.975))
fec_CI <- quantile(store_ltr[,3], probs=c(0.025, 0.975))


## Calc quantile differences in lambda 
DvN <- quantile(log(matLam[,1]) - log(matLam[,5]), probs=c(0.025, 0.975))
## ---------------------------------------------------------------------------------------------------





## FOR DUNE AND NONDUNE SOURCE IN  NONDUNE HABITAT ---------------------------------------------------
## Calculate mean FC value across different sources
mean.emrg.nd <- mean(c(prEmrg$predicted[2], prEmrg$predicted[6]))
mean.surv.nd <- mean(c(prSurv$predicted[2], prSurv$predicted[6]))
mean.fec.nd <- mean(c(prSds$predicted[2], prSds$predicted[6]))


## Calculate sensitivities
vr.nd <- list(emrg=mean.emrg.nd, surv=mean.surv.nd, fec=mean.fec.nd, sdSurv=mean.sdSurv )
lambda.nd <- expression(emrg * surv * fec * sdSurv)
x.nd <- vitalsens(lambda.nd, vr.nd)


## Calculate change in fitness components
## Local FC - foreign FC 
delta.emrg.nd <- prEmrg$predicted[6] - prEmrg$predicted[2]
delta.surv.nd <- prSurv$predicted[6] - prSurv$predicted[2] 
delta.fec.nd <- prSds$predicted[6] - prSds$predicted[2]


## LTRE barplot
ltr.emrg.nd <- delta.emrg.nd * x.nd$sensitivity[1]
ltr.surv.nd <- delta.surv.nd * x.nd$sensitivity[2]
ltr.fec.nd <- delta.fec.nd * x.nd$sensitivity[3]
ltr.sdSurv.nd <- 0

diff.lam.nd <- lam.mod[6] - lam.mod[2]
ltr.adj.nd <- c(ltr.emrg.nd, ltr.surv.nd, ltr.fec.nd, ltr.sdSurv.nd, diff.lam.nd)




## Calculate 95% quantiles ---------------------------------------------------------
## Loop over FC matrices and re-calculate sensitivity and delta FC
store_ltr.nd <- matrix(NA, nrow=nrow(matEmrg), ncol=3)

for (cc in 1:nrow(matEmrg)) {
  
  ## Calc mean FC
  mean.emrg.nd <- mean(c(matEmrg[cc,2],matEmrg[cc,6]))
  mean.surv.nd <- mean(c(matSurv[cc,2],matSurv[cc,6]))
  mean.fec.nd <- mean(c(matSds[cc,2],matSds[cc,6]))
  
  ## Calc sensitivities
  vr.nd <- list(emrg=mean.emrg.nd, surv=mean.surv.nd, fec=mean.fec.nd, sdSurv=mean.sdSurv )
  lambda.nd <- expression(emrg * surv * fec * sdSurv)
  x.nd <- vitalsens(lambda.nd, vr.nd)  
  
  ## Calc delta FC
  delta.emrg.nd <- matEmrg[cc,6] - matEmrg[cc,2]
  delta.surv.nd <- matSurv[cc,6] - matSurv[cc,2] 
  delta.fec.nd <- matSds[cc,6] - matSds[cc,2]
  
  
  ## LTRE (product of sensitivity and change in FC)
  ltr.emrg.nd <- delta.emrg.nd * x.nd$sensitivity[1]
  ltr.surv.nd <- delta.surv.nd * x.nd$sensitivity[2]
  ltr.fec.nd <- delta.fec.nd * x.nd$sensitivity[3]
  
  store_ltr.nd[cc,1] <- ltr.emrg.nd 
  store_ltr.nd[cc,2] <- ltr.surv.nd 
  store_ltr.nd[cc,3] <- ltr.fec.nd 
}


## Calculate 95% quantiles
emrg_CI.nd <- quantile(store_ltr.nd[,1], probs=c(0.025, 0.975))
surv_CI.nd <- quantile(store_ltr.nd[,2], probs=c(0.025, 0.975))
fec_CI.nd <- quantile(store_ltr.nd[,3], probs=c(0.025, 0.975))


## Calculate quantiles differences in lambda 
NvD.n <- quantile(log(matLam[,6]) - log(matLam[,2]), probs=c(0.025, 0.975))








## PLOT LTRE -----------------------------------------------------------------------------------------

## Assign colours
colD <- rgb(t(col2rgb("orange")),alpha=215,maxColorValue = 255)
col_plot <- c(colD, colD, colD, colD, "grey")

colN <- rgb(t(col2rgb("royalblue")),alpha=215,maxColorValue = 255)
col_plot.N <- c(colN, colN, colN, colN, "grey")


svg('GSD_fig2.svg', width=17, height=12)
par(mfrow=c(2,1))  #Plot in 2 panels
par(mar=c(2.5, 6, 2, 0),xpd=FALSE)  
windowsFonts(A = windowsFont("Times New Roman"))


## In dune habitat -----------------
bar.centers <- barplot(t(matrix(ltr.adj, ncol=1, byrow=TRUE,
                                dimnames=list(c("Emergence", "Survival","Fecundity", "Seed survival", "Delta lambda")))), 
                       space=c(0,0.1), cex.lab=2.25, border=FALSE, cex.names=1.8,
                       xlab=NA, col=col_plot, beside=TRUE, main="Dune vs non-dune ecotype in dune habitat", 
                       ylab="Contribution to total local minus\nforeign ecotype fitness difference",
                       legend.text = FALSE,args.legend = list(x = "top",cex=2.25), ylim=c(-2,7.5),
                       cex.lab=1.5, cex.main=2, font.main=6, family="A")

abline(h=0)


## Add error bars to barplot
CI_up <- c(emrg_CI[2],surv_CI[2],fec_CI[2],0,DvN[2]) 
CI_lo <- c(emrg_CI[1],surv_CI[1],fec_CI[1],0,DvN[1]) 

arrows(bar.centers, CI_up, bar.centers, CI_lo, lwd = 1.5, angle = 90,
       code = 3, length = 0)





## In non-dune habitat -----------
bar.centers.nd <- barplot(t(matrix(ltr.adj.nd, ncol=1, byrow=TRUE, 
                                   dimnames=list(c("Emergence", "Survival","Fecundity", "Seed survival", "Delta lambda")))), 
                          space=c(0,0.1), cex.lab=2.25, border=FALSE, cex.names=1.8,
                          xlab=NA, col=col_plot.N, beside=TRUE, main="Dune vs non-dune ecotype in non-dune habitat", 
                          ylab="Contribution to total local minus\nforeign ecotype fitness difference", 
                          legend.text=FALSE,args.legend=list(x = "bottomleft"), ylim=c(-30,42),
                          cex.lab=1.5, cex.main=2, font.main=6, family="A") 
abline(h=0)

## Add error bars to barplot
CI_up.nd <- c(emrg_CI.nd[2],surv_CI.nd[2],fec_CI.nd[2],0,NvD.n[2])
CI_lo.nd <- c(emrg_CI.nd[1],surv_CI.nd[1],fec_CI.nd[1],0,NvD.n[1])

arrows(bar.centers.nd, CI_up.nd, bar.centers.nd, CI_lo.nd, lwd = 1.5, angle = 90,
       code = 3, length = 0)
dev.off()
## ------------------------------------------------------------------------------------------------






## PLOT DISTRIBUTIONS OF FITNESS COMPONENT CONTRIBUTIONS TO LAMBDA --------------------------------

col_bar <- rgb(t(col2rgb("grey")), alpha=255, maxColorValue=255)
par(mfrow=c(3,2))  #Plot in 6 panels
par(mar=c(2.5, 6, 2, 2),xpd=FALSE)  


## In dune habitat
hist(store_ltr[,1],breaks=80, main="Contribution of emergence in dunes", xlim=c(-35,40),
     col=col_bar, border=FALSE, family="A", font.main=1, cex.lab=2, cex.main=2)
abline(v=emrg_CI[1])
abline(v=emrg_CI[2])

hist(store_ltr[,2], breaks=30, main="Contribution of survival in dunes", xlim=c(-35,40),
     col=col_bar, border=FALSE, family="A", font.main=1, cex.lab=2, cex.main=2)
abline(v=surv_CI[1])
abline(v=surv_CI[2])

hist(store_ltr[,3],breaks=85, main="Contribution of fecundity in dunes", xlim=c(-35,40),
     col=col_bar, border=FALSE, family="A", font.main=1, cex.lab=2, cex.main=2)
abline(v=fec_CI[1])
abline(v=fec_CI[2])



## In non-dune habitat
hist(store_ltr.nd[,1],breaks=350, main="Contribution of emergence in nondunes", xlim=c(-35,40),
     col=col_bar, border=FALSE, family="A", font.main=1, cex.lab=2, cex.main=2)
abline(v=emrg_CI.nd[1])
abline(v=emrg_CI.nd[2])

hist(store_ltr.nd[,2], breaks=350, main="Contribution of survival in nondunes", xlim=c(-35,40),
     col=col_bar, border=FALSE, family="A", font.main=1, cex.lab=2, cex.main=2)
abline(v=surv_CI.nd[1])
abline(v=surv_CI.nd[2])

hist(store_ltr.nd[,3],breaks=500, main="Contribution of fecundity in nondunes", xlim=c(-35,40),
     col=col_bar, border=FALSE, family="A", font.main=1, cex.lab=2, cex.main=2)
abline(v=fec_CI.nd[1])
abline(v=fec_CI.nd[2])
## ---------------------------------------------------------------------------------------------------------  





