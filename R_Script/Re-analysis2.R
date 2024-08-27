###########################
########################### KOALA INBREEDING AND DISEASE 
########################### Re-analysis post-review
########################### Kasha Strickland

setwd("~/Documents/Papers in prep/Submitted/Romane_diseasekoalainbreeding/March2020_updateanalyses")

library(adehabitatHR)
library(maptools)
library(rgdal)
library(spatstat)
#library(devtools) # only need for install
#install_git("git://github.com/gsk3/taRifx.geo.git")
#library(taRifx.geo)
library(rgeos)
library(raster)
library(dplyr)
library(corpcor)
library(vegan)
library(gdata)
library(brms)
library(MCMCglmm)
library(rstan)

list.files()

###load new g matrix
g2<-read.csv("Rmatrix_new2.csv" ,stringsAsFactors = FALSE)
#ids<-read.csv("final_koalas.csv")
str(g2)
####make the R matrix
summary(g2)
table1<-aggregate(g2$new1,list(g2$new1),length)
table2<-aggregate(g2$new2,list(g2$new2),length)
names(table1)<-c("Name", "Number")
table1<-table1[order(table1$Number,decreasing=TRUE),]
names(table2)<-c("Name", "Number")
table2<-table2[order(table2$Number,decreasing=FALSE),]
#Reorder rows of the data frame to allow filling the matrix correctly
g2b<-g2[order(match(g2$new2, table2$Name)),]
g2b<-g2b[order(match(g2b$new1, table1$Name)),]
nam=unique(c(as.character(table1$Name),as.character(table2$Name)))
#Create an empty matrix and fill it
gm<-matrix(0,442,442)
rownames(gm)<-colnames(gm)<-nam
gm[lower.tri(gm,diag=FALSE)]<-g2b$quellergt
diag(gm)<-1.001
gm[upper.tri(gm)]<-t(gm)[upper.tri(gm)]
g2[g2$new1=="Abby" & g2$new2=="STARKEY",]#check it works
gm[which(colnames(gm)=="STARKEY"),which(rownames(gm)=="Abby") ]##check it works
det(gm)#if zero, matrix is singular

##models require that our G is positive definite. check here
is.positive.definite(gm)
#is.positive.definite(GRM_old)
##ok so it's not - we can convert it to the closest doing this
gm2<-make.positive.definite(gm)
#GRM_old2<-make.positive.definite(GRM_old)
is.positive.definite(gm2)
hist(gm2)
##great - but check that we haven't lost the structure of the data by correlating this new matrix with the old one
#mantel(GRM_new,GRM_new2)# r = 0.971; plot it to see
gmsN1<-as.data.frame(unmatrix(gm)) 
gmsN2<-unmatrix(gm2)
var(gmsN2)##check variance in R values 
plot(gmsN1$`unmatrix(gm)`,gmsN2)##looks good!

unmatrix_without_diag <- function(mat) {
  non_diag_indices <- !diag(nrow(mat))
  as.vector(mat[non_diag_indices])
}

# Apply function to matrices
gmsN1 <- as.data.frame(unmatrix_without_diag(gm))
gmsN2 <- unmatrix_without_diag(gm2)

plot(gmsN1$`unmatrix_without_diag(gm)`,gmsN2)##looks good!

######################################################
####OK now we need to make a home range overlap matrix
##Because I have made it already, just read in the object as follows (to save space and time)
load("UDOI.R")
##script to learn how to make above HRO object 
#full<-read.csv("Final_spatial_subsetted.csv",stringsAsFactors = FALSE)
##subset for at least 5 points
#HRdata<-subset(full,full$relocations>4)
#HRdata<-subset(HRdata,!HRdata$Name=="Anna")##her points have issues and make a really innacurate home range so removing
#str(HRdata)
#HRdata$Name<-factor(HRdata$Name) #<- this is for selecting only animal name, X and Y, from dataframe
#xydata<-data.frame(Name=HRdata$Name,X=HRdata$UTM_Eastin,Y=HRdata$UTM_Northi)
#xydata$Y<-as.numeric(xydata$Y)  # need Name as factor and Easting/Northing as numeric for future
#xydata$X<-as.numeric(xydata$X)
#class(xydata$X) #CHECKS
#class(xydata$Y) #CHECKS
#class(xydata$Name) #CHECKS
#setting grid resolution to 50 m here
#x <- seq(min(xydata[,"X"])-5000,max(xydata[,"X"])+5000,by=50) # where resolution is the pixel size you desire. 
#y <- seq(min(xydata[,"Y"])-5000,max(xydata[,"Y"])+5000,by=50) #This gives a 10m x 10m grid, with a 100m buffer
#xy <- expand.grid(x=x,y=y)
#coordinates(xy) <- ~x+y
#gridded(xy) <- TRUE
#class(xy)
# create spatialpointsDF needed for UD
#hrxydata<-SpatialPointsDataFrame(xydata[,2:3],xydata["Name"])
#plot(hrxydata)
#NB MAKE SURE EACH INDIV HAS AT LEAST 5 POINTS TO FIT A HOME RANGE
# create UD for each indiv with chosen smoothing parameter (h)
#ud_href<-kernelUD(hrxydata[,1],grid=xy)##using href smoothing parameter because after visual checks, looks like it's doing a good job
#image(ud_href)
##calculate home range overlap
#hrud<-kerneloverlaphr(ud_href,method = c("UDOI"),percent=95)

###Figure 1
diag(gm)<-NA
pdf("Figure 1.pdf",width=10)
par(mfrow=c(1,2))
hist(dat2$IR_new,xlab="Internal relatedness",main="",breaks=20)
hist(gm,xlab="Pairwise relatedness",main="",breaks=20)
dev.off()

##OK ready to model!
########1st model - chlamydia at first catch
list.files()
dat<-read.csv("M1.csv",stringsAsFactors = FALSE)
str(dat)
dat2<-na.omit(dat)
str(dat2)
dat2$Name2<-dat2$Name

###Check whether HRO and G correlate too strongly to be put in model together
######################################
##and now hro matrix
HRO1<-hrud[rownames(hrud)%in%dat2$Name,colnames(hrud)%in%dat2$Name]
dim(HRO1)
is.positive.definite(HRO1)
GRM1<-gm2[rownames(gm2)%in%dat2$Name,colnames(gm2)%in%dat2$Name]
dim(GRM1)
is.positive.definite(GRM1)
#GRM2<-make.positive.definite(GRM1)
diag(HRO1)=1
diag(GRM1)=1
##check the correlation between g and hro
HROs<-as.data.frame(unmatrix(HRO1))
GRMs<-as.data.frame(unmatrix(GRM1))
cor(GRMs,HROs)##r2 = 0.3002359
cor.test(GRMs$`unmatrix(GRM1)`,HROs$`unmatrix(HRO1)`)
plot(GRMs$`unmatrix(GRM1)`,HROs$`unmatrix(HRO1)`, xlab ="pairwise relatedness", ylab = "pairwise home range overlap")
var(HROs$`unmatrix(HRO1)`)
var(GRMs$`unmatrix(GRM1)`)

############################################################
############run full model (change = y variable is now chlamydia sick at first capture)
### G = Va; S = shared environment; M = maternal effect

##check covariance of fixed effects
pairs(dat2[c('Chlamydia_sick','Sex','Age','IR_new','Breeding')])
cor(dat2[c('Chlamydia_sick','Sex','Age','IR_new','Breeding')])

##in MCMCglmm
Ginv<-solve(gm2)
Hinv<-solve(hrud)
Ginvsp<-as(Ginv,'CsparseMatrix')
Hinvsp<-as(Hinv,'CsparseMatrix')

priorM1<-list(G=list(G1=list(V=1,nu=1000,alpha.mu=0,alpha.V=1),
                     G2=list(V=1,nu=1000,alpha.mu=0,alpha.V=1)),
                     R=list(V=1,fix=1))

mcM1_GS<-MCMCglmm(Chlamydia_sick~Sex+Age*IR_new+Breeding,
               random=~Name+Name2,
               ginverse = list(Name=Ginvsp,Name2=Hinvsp),
               family="threshold",
               trunc=TRUE,pr=TRUE,
               prior=priorM1,
               data=dat2,
               nitt=1030000,thin=1000,burnin=30000)
plot(mcM1_GS)
autocorr(mcM1_GS$VCV)
summary(mcM1_GS)

varR<-as.data.frame(mcM1_GS$VCV)
names(varR)
h2<-varR$Name/(varR$Name+varR$Name2+varR$units)
h2
hist(h2)
mean(h2)
quantile(h2,c(0.05,0.95))


gm2<-gm2[rownames(gm2) %in% inds,colnames(gm2) %in% inds]

gmsN1<-as.data.frame(unmatrix(gm)) 
gmsN2<-unmatrix(gm2)
plot(gmsN1$`unmatrix(gm)`,gmsN2)
