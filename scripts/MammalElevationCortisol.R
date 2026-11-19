#install.packages("fromhere")
#install.packages("lmtest")
#install.packages("phylosignal")
#install.packages("broom.mixed")  
library(broom.mixed)
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

##Normalizing Data
FinalMammal<-FinalMammal |> mutate(logCortisol=log(Combined_BC_Mean))

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


#Linear Regression Model without phylogenetic signal

model<-lm(Combined_BC_Mean~Elevation,data=FinalMammal)
summary(model)
#Checking linear model parameters
par(mfrow=c(2,2)) 
plot(model)
model<-lm(logCortisol~Elevation,data=FinalMammal)
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


hist(FinalMammal$logCortisol)

#glance(model)

##Creating a tree

#species <- c("Connochaetes taurinus", "Ctenomys talarum", "Lynx pardinus", "Macaca mulatta", "Octodon degus", "Ornithorhynchus anatinus", "Ovibos moschatus", "Papio cynocephalus", "Tamias amoenus", "Tragelaphus strepsiceros", "Urocitellus parryii", "Ursus maritimus")
#resolved<-tnrs_match_names(species)
#tr<-tol_induced_subtree(ott_ids=ott_id(resolved))
#my_tree <- tol_induced_subtree(ott_ids = resolved$ott_id)
#plot(tr)
#head(tr)
#class(tr)
#class(my_tree)



####struggling with tree
#install.packages("rotl")
#install.packages("ape")
#library(rotl)
#library(ape)
#library(nlme)

#species <- c("Connochaetes taurinus", "Ctenomys talarum", "Lynx pardinus", "Macaca mulatta", "Octodon degus", "Ornithorhynchus anatinus", "Ovibos moschatus", "Papio cynocephalus", "Tamias amoenus", "Tragelaphus strepsiceros", "Urocitellus parryii", "Ursus maritimus")
#taxa <- tnrs_match_names(species)

#class(tr)

#tree <- tol_induced_subtree(ott_ids = taxa$ott_id, file="mammal_p")
#plot(tree)
#tree$tip.label <- strip_ott_ids(tree$tip.label, remove_underscores = FALSE)
#tree$tip.label
#dat$species <- gsub(" ", "_", dat$species)
#setdiff(dat$species, tree$tip.label)   
#setdiff(tree$tip.label, dat$species) 
#class(tree)
#str(tr)




##retry
#phylogentic tree

library(rotl)
library(ape)
FinalMammal<-FinalMammal |> mutate(GenSpec=paste(Genus,Species,sep="_"))
                                   
species <- c("Connochaetes taurinus", "Ctenomys talarum", "Lynx pardinus", "Macaca mulatta", "Octodon degus", "Ornithorhynchus anatinus", "Ovibos moschatus", "Papio cynocephalus", "Tamias amoenus", "Tragelaphus strepsiceros", "Urocitellus parryii", "Ursus maritimus")
taxa <- tnrs_match_names(species)
taxa


tree <- tol_induced_subtree(ott_ids = taxa$ott_id)

# clean tip labels
tree$tip.label <- strip_ott_ids(tree$tip.label, remove_underscores = FALSE)

# add branch lengths
tree <- multi2di(tree)
tree <- compute.brlen(tree, method = "Grafen")

is.ultrametric(tree)   
#plot(tree)
class(tree)

######ATTENTION TAs AND JOEY, PLEASE READ LINE BELOW THIS.
##please note that the following 6 lines may crash your R. To obtain the results from the repository run "readRDS" line below.
##Running Model
#pglsModel <- gls(Combined_BC_Mean ~ Elevation, correlation = corBrownian(phy = tree),
                 data = FinalMammal, method = "ML")



pglsModel <- gls(
  Combined_BC_Mean ~ Elevation,
  correlation = corPagel(value = 1.099927, phy = tree, form = ~Species),
  data = FinalMammal,
  method = "ML"
)
pglsModel <- gls(Combined_BC_Mean ~ Elevation, correlation = corCAR1(lambda=1.099927),
                 data = FinalMammal, method = "ML")
pglscorMartins <- gls(Combined_BC_Mean ~ Elevation, correlation = corMartins(value=1,phy=tree),
                 data = FinalMammal, method = "ML")
summary(pglsModel)
pglsOutput <- tidy(pglsModel, conf.int = TRUE)
glimpse(pglsOutput)
saveRDS(pglsModel, file = "pgls.output.rds")

pgls.output<-readRDS(fromhere::from_rproj("scripts/pgls.Output.rds")) ##obtaining model output from direct model path.
#Evaluating output and generating confidence intervals.
pglsOutput <- tidy(pgls.output, conf.int = TRUE)
glimpse(pglsOutput)


vcvtree<-vcv.phylo(tree)
class(tree)
treenum<-as.numeric(tree)

lambdaTest(x=FinalMammal$Combined_BC_Mean,vcv=p4tree)
lambdaTest(x=tree,vcv)

str(vcvtree)
p4tree<-as.matrix(tree)

library(phylobase)
phyloSignal(,methods="Lambda")
phylobase::hasEdgeLength(rand_p4_tree)
sapply(phylobase::tdata(rand_p4_tree),class)

p4tree <- as(tree, "phylo4d") # Converted to phylo4d

class(p4tree)
phyloSignal



inherits(p4tree,"phylo")
inherits(p4tree,"phylo4d") # verifies object type

##What actually worked
x<-FinalMammal$Combined_BC_Mean
names(x)<-rownames(FinalMammal)
x<-x[rownames(vcvtree)]
lambdaTest(x,vcvtree)

install.packages("performance")
library(performance)
install.packages("see")
library(see)
check_model(pgls.output)
check_model(model)


norm_res<-residuals(pglmmModel,type="normalized")
qqnorm(norm_res);qqline(norm_res)
AIC(pgls.output)

pglsModel <- gls(
  Combined_BC_Mean ~ Elevation,
  correlation = corPagel(value = 1, phy = tree, form = ~Species),
  data = FinalMammal,
  method = "ML"
)
install.packages("phyr")
library(phyr)
pglmmModel <- pglmm(
  Combined_BC_Mean ~ Elevation + (1 | GenSpec__), 
  data = FinalMammal, 
  cov_ranef = list(GenSpec = tree)
)
summary(pglmmModel)
install.packages("DHARMa")
library(DHARMa)
sim_res<-simulateResiduals(fittedModel=pglmmModel,plot=TRUE)
sim_lm<-simulateResiduals(fittedModel=model,plot=TRUE)


unique(FinalMammal$GenSpec)
print(tree$tip.label)
identical(FinalMammal$GenSpec,tree$tip.label)
which(!(FinalMammal$GenSpec %in% tree$tip.label))
