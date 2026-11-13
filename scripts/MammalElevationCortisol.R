library(fromhere)
library(tidyverse)
library(WebPower)
library(lmtest)
library(broom)
library(ggplot2)
library(rotl)
library(phylosignal)
library(nlme)
library(ape)
HorBase<-read.csv(fromhere::from_rproj("data/raw/mammalraw.csv"))

##Clean the data
MammalData <- HorBase|>
  filter(Vert_Group=="Mammal",CORT=="cortisol",!Genus%in%c("Balaenoptera","Trichechus","Phoca","Halichoerus","Spermophilus"),!(F_BC_N=="null"&M_BC_N=="null")) |> 
  dplyr::select(Species, F_BC_Mean,M_BC_Mean,Common_name,F_BC_N,M_BC_N,Genus,Elevation)
##Mammal Cortisol Data


SavedMammal<-MammalData
MammalData$F_BC_Mean<-as.numeric(MammalData$F_BC_Mean)
MammalData$F_BC_N<-as.numeric(MammalData$F_BC_N)
MammalData$M_BC_Mean<-as.numeric(MammalData$M_BC_Mean)
MammalData$M_BC_N<-as.numeric(MammalData$M_BC_N)


MammalData[is.na(MammalData)]<-0
Newcolumn<-MammalData |> mutate(Combined_BC_Mean=((F_BC_Mean*F_BC_N)+(M_BC_Mean*M_BC_N))/(F_BC_N+M_BC_N))
Newcolumn$F_BC_Mean=SavedMammal$F_BC_Mean
Newcolumn$M_BC_Mean=SavedMammal$M_BC_Mean

FinalMammal<-Newcolumn


FinalMammal$Elevation<-as.numeric(FinalMammal$Elevation)
FinalMammal$Combined_BC_Mean<-as.numeric(FinalMammal$Combined_BC_Mean)

cort_by_elevation<-FinalMammal |> group_by(Elevation) |> summarise(Cortisol=mean(Combined_BC_Mean))

ggplot(cort_by_elevation)+geom_point(aes(x=Elevation,y=Cortisol))+
  geom_smooth(aes(x=Elevation,y=Cortisol),method=lm,se=FALSE)+
  theme_bw()+
  ylab("Cortisol (ng/mL)")+
  xlab ("Elevation (m)")



#Linear regression power analysis
wp.regression(n = 56, p1 = 17, f2 = 0.48197, alpha = 0.05, power = NULL)


#Spearman's Test
model2<-cor.test(x=cort_by_elevation$Elevation,y=cort_by_elevation$Cortisol,method='spearman')

#Linear Regression Model

model<-lm(Combined_BC_Mean~Elevation,data=FinalMammal)
summary(model)
#Checking linear model parameters
par(mfrow=c(2,2)) 
plot(model)
##plotting model
ggplot(model)+geom_point(aes(x=Elevation,y=Combined_BC_Mean))+
  ylab("Mean Cortisol (ng/mL)")+
  xlab("Elevation (m)")



MammalOutput <- tidy(model, conf.int = TRUE)
MammalOutput
glance(model)

##Creating a tree

species <- c("Connochaetes taurinus", "Ctenomys talarum", "Lynx pardinus", "Macaca mulatta", "Octodon degus", "Ornithorhynchus anatinus", "Ovibos moschatus", "Papio cynocephalus", "Tamias amoenus", "Tragelaphus strepsiceros", "Urocitellus parryii", "Ursus maritimus")
resolved<-tnrs_match_names(species)
tr<-tol_induced_subtree(ott_ids=ott_id(resolved))
plot(tr)
head(tr)
class(tr)


#pgls

pagelslambda<-lambdaTest(FinalMammal$Elevation,tr)


####struggling with tree
install.packages("rotl")
install.packages("ape")
library(rotl)
library(ape)

species <- c("Connochaetes taurinus", "Ctenomys talarum", "Lynx pardinus", "Macaca mulatta", "Octodon degus", "Ornithorhynchus anatinus", "Ovibos moschatus", "Papio cynocephalus", "Tamias amoenus", "Tragelaphus strepsiceros", "Urocitellus parryii", "Ursus maritimus")
taxa <- tnrs_match_names(species)
class(taxa)

tree <- tol_induced_subtree(ott_ids = taxa$ott_id, file="mammal_p")
plot(tree)
tree$tip.label <- strip_ott_ids(tree$tip.label, remove_underscores = FALSE)
tree$tip.label
dat$species <- gsub(" ", "_", dat$species)
setdiff(dat$species, tree$tip.label)   
setdiff(tree$tip.label, dat$species) 
class(tree)

pglsModel <- gls(Combined_BC_Mean ~ Elevation, correlation = corBrownian(phy = tr),
                 data = FinalMammal, method = "ML")
(pglsModel)

