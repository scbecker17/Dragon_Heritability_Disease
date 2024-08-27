##Working MCMCglmm
library(MCMCglmm); library(poppr); library(tidyverse); library(gdata) 
library(dartR); library(corpcor);library(QGglmm);library(purrr)

name_to_rm <- c("Kim",
                "Vicky",
                "Gertie",
                "Telly",
                "Byron",
                "Brooke",
                "Mitch",
                "Squid",
                "Mojo",
                "Annie",
                "Hagrid",
                "December",
                "Tinly",
                "Tim",
                "Jungle",
                "Scruffy.1",
                "Chanel.1")

#Making GRM
g2<-read.csv("Data/R_estimates_drag.csv" ,stringsAsFactors = FALSE)
table1<-aggregate(g2$ind1.id,list(g2$ind1.id),length)
table2<-aggregate(g2$ind2.id,list(g2$ind2.id),length)
names(table1)<-c("Name", "Number")
table1<-table1[order(table1$Number,decreasing=TRUE),]
names(table2)<-c("Name", "Number")
table2<-table2[order(table2$Number,decreasing=FALSE),]
g2b<-g2[order(match(g2$ind2.id, table2$Name)),]
g2b<-g2b[order(match(g2b$ind1.id, table1$Name)),]
nam=unique(c(as.character(table1$Name),as.character(table2$Name)))
gm<-matrix(0,880,880)
rownames(gm)<-colnames(gm)<-nam
gm[lower.tri(gm,diag=FALSE)]<-g2b$quellergt
diag(gm)<-1.001
gm[upper.tri(gm)]<-t(gm)[upper.tri(gm)]
gm <- gm[!(rownames(gm) %in% name_to_rm),
         !(colnames(gm) %in% name_to_rm)]
indsm <- rownames(gm)
off_diagonals <- gm[upper.tri(gm, diag = FALSE)]
var(off_diagonals)

#Read in phenotypic data
dragdisease <- read.csv("Data/FullData.csv", header = T)
str(dragdisease)

#Filter data
dragdisease <- dragdisease %>% filter(Season > 4, Sex != "Unknown", AgeClass != " ")

# Dragons common to phenotypic data and GRM
inds <- intersect(dragdisease$Name, indsm)
# write.csv(inds, "Output/Names_for_GRM.csv")

# Subset phenotype data  
disdragons <- dragdisease[dragdisease$Name %in% inds,]
length(unique(disdragons$Name))

# Subset GRM and make positive definite
# gm<-gm[rownames(gm) %in% inds,colnames(gm) %in% inds]
ginv<-MASS::ginv(gm) #invert GRM
is.positive.definite(ginv)
ginv2<-make.positive.definite(ginv) #get pos def of inverse 
correlation <- cor(c(ginv), c(ginv2), method = "pearson") #check correlation
plot(as.vector(ginv),as.vector(ginv2)) #visualise correlation
Ginvsp<-as(ginv2,"dgCMatrix") #convert to sparse
dimnames(Ginvsp)<-list(rownames(gm),colnames(gm))
Ginvsp <- Ginvsp[inds, inds]

# gmsN1<-as.data.frame(unmatrix(ginv)) 
# gmsN2<-unmatrix(ginv2)
# plot(gmsN1$`unmatrix(ginv)`,gmsN2)


#Make sure names are same in GRM and phenotypic data
# gp_nam <- (rownames(gm2))
# dis_nam <- (as.character(disdragons$Name))
# 
# unique_to_gp = setdiff(gp_nam, dis_nam)
# unique_to_dis = setdiff(dis_nam, gp_nam)
# 
# cat("Names unique to gp:\n")
# cat(unique_to_gp, sep = "\n")
# 
# cat("\nNames unique to dis:\n")
# cat(unique_to_dis, sep = "\n")

#Priors
prior1 <- list(
  G = list(G1 = list(V=1,nu=1000,alpha.mu=0,alpha.V=1), 
           # G2 =  list(V=1,nu=1000,alpha.mu=0,alpha.V=1), 
           G2 = list(V=1,nu=0.002),
           # G3 = list(V=1,nu=1000,alpha.mu=0,alpha.V=1)),
           G3 = list(V=1,nu=0.002)),
  R=list(V=1, fix=1)
)


colSums(is.na(disdragons))
disdragons2 <- disdragons %>%
  filter(!is.na(InfectedDensity.Annual))
str(disdragons2)
disdragons2$Sex <- as.factor(disdragons2$Sex)
disdragons2$Name <- as.factor(disdragons2$Name)
disdragons2$AgeClass <- as.factor(disdragons2$AgeClass)
disdragons2$Season <- as.factor(disdragons2$Season)
disdragons2$Name2 <- disdragons2$Name
disdragons2$Dragon.Year <- year(disdragons2$Date)
str(disdragons2)

# plot(disdragons2$Dragon.Year, disdragons2$Time.in.pop)
# cor(disdragons2$Dragon.Year, disdragons2$Time.in.pop)

Dragons1 <- 
  disdragons2 %>% 
  group_by(Name, Dragon.Year, Sex, AgeClass, ActiveFungus) %>% 
  filter(Sex != 'Unknown') %>%
  summarise_all(~mean(.x, na.rm = T)) %>% 
  ungroup()

length(unique(Dragons1$Name))

Dragons1 %>% 
  filter(!Sex == "Unknown", !AgeClass == "Unknown") %>% 
  arrange(Dragon.Year) %>% 
  dplyr::select(Name, Dragon.Year, Sex, Time.in.pop, AgeClass,
                ActiveFungus, InfectedDensity.Annual) %>% #na.omit %>% 
  ungroup %>% droplevels -> 
  TestDF1
TestDF1$Name2 <- TestDF1$Name
write.csv(TestDF1, "Data/One_per_season.csv", row.names = F)
TestDF1 <- read.csv("Data/One_per_season.csv", header = T)


#Animel Model
mcM1_GS<-MCMCglmm(ActiveFungus ~ 1, #Sex + Time.in.pop + InfectedDensity.Annual, 
                  random = ~Name + Name2 + Dragon.Year,
                  ginverse = list(Name = Ginvsp),
                  family = "threshold",
                  trunc = TRUE, pr = TRUE,
                  prior = prior1,
                  data = disdragons2,
                  nitt = 1700000, thin = 1000, burnin = 30000)

save(mcM1_GS, file = "Output/mcmc_GRM+Name.rds")
load( "Output/mcmc_GRM+Name+DYear(Random_Only).rds")
save(mcM1_GS, file = "Output/mcmc_GRM+Name+DYear(Random_Only).RData")
plot(mcM1_GS$VCV)
autocorr(mcM1_GS$VCV)
summary(mcM1_GS)
heidel.diag(mcM1_GS$VCV)
effectiveSize(mcM1_GS$VCV)

# latent scale h^2
h2_l <- mcM1_GS[["VCV"]][ , "Name"] / rowSums(mcM1_GS[["VCV"]]) 
mean(h2_l) 
HPDinterval(h2_l)



# observed data scale h^2
pr <-  purrr::pmap_dfr(list(mu = mcM1_GS[["Sol"]][ , "(Intercept)"], 
                            var.a = mcM1_GS[["VCV"]][ , "Name"], 
                            var.p = rowSums(mcM1_GS[["VCV"]])-1), 
                       QGparams, model = "binom1.probit", verbose = FALSE) 
mean(pr[["h2.obs"]]) 
HPDinterval(as.mcmc(pr[["h2.obs"]]))


