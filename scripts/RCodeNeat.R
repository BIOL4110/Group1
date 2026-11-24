##Packages
library(phyr)
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
##Importing Dataset
HorBase<-read.csv(fromhere::from_rproj("data/raw/mammalraw.csv"))
##Cleaning Dataset
MammalData <- HorBase|>
  filter(Vert_Group=="Mammal",CORT=="cortisol",!Genus%in%c("Balaenoptera","Trichechus","Phoca","Halichoerus","Spermophilus"),!(F_BC_N=="null"&M_BC_N=="null")) |> 
  dplyr::select(Species, F_BC_Mean,M_BC_Mean,Common_name,F_BC_N,M_BC_N,Genus,Elevation)

###We now need to obtain combined baseline cortisol averages rather than keeping them separated by sex

#Let's Store the data that we've already obtained
SavedMammal<-MammalData
#Convert our data to be able to do operations on
MammalData$F_BC_Mean<-as.numeric(MammalData$F_BC_Mean)
MammalData$F_BC_N<-as.numeric(MammalData$F_BC_N)
MammalData$M_BC_Mean<-as.numeric(MammalData$M_BC_Mean)
MammalData$M_BC_N<-as.numeric(MammalData$M_BC_N)

#Math to form new combined column
#We are renaming this new data frame to "New Column".
MammalData[is.na(MammalData)]<-0
Newcolumn<-MammalData |> mutate(Combined_BC_Mean=((F_BC_Mean*F_BC_N)+(M_BC_Mean*M_BC_N))/(F_BC_N+M_BC_N))
Newcolumn$F_BC_Mean=SavedMammal$F_BC_Mean
Newcolumn$M_BC_Mean=SavedMammal$M_BC_Mean

#Let's use the name "Final Mammal" for clarity as we proceed.
FinalMammal<-Newcolumn
#Create a tree to account for species non-independence
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
class(tree)


#pglmmModel <- pglmm(
#Combined_BC_Mean ~ Elevation + (1 | GenSpec__), 
#data = FinalMammal, 
#cov_ranef = list(GenSpec = tree),
#REML="ML"
#) 

pglmmModel <- pglmm(
  Combined_BC_Mean ~ Elevation + (1 | FactGenSpec__), 
  data = FinalMammal, 
  cov_ranef = list(FactGenSpec = tree),
  REML=FALSE) 

FinalMammal<-FinalMammal |> mutate(logCort=log(Combined_BC_Mean))
pglmmModellog <- pglmm(
  logCort ~ Elevation + (1 | FactGenSpec__), 
  data = FinalMammal, 
  cov_ranef = list(FactGenSpec = tree),
  REML=FALSE) 
summary(pglmmModellog)

summary(pglmmModel)
##Plot
pglmmplot <- data.frame(
  Elevation = seq(min(FinalMammal$Elevation, na.rm = TRUE),
                  max(FinalMammal$Elevation, na.rm = TRUE),
                  length.out = 100)
)
X_new <- model.matrix(~ Elevation, data = pglmmplot)
pglmmplot$fit <- as.vector(X_new %*% pglmmModel$B)


ggplot(FinalMammal, aes(x = Elevation, y = Combined_BC_Mean)) +
  geom_point() +
  geom_line(data = pglmmplot, aes(x = Elevation, y = fit),colour = "blue", linewidth = 1) 

#Okay after hella googling I found out that everything sucks and the NLME is a dumbass R package that doesnt work if you have multiples
#of the same species so all that stuff was useless and we just have to do a phylogenetic mixed model.



##Everything below this line is bs

##Attempting CI & R2
est <- 0.0020147
se  <- 0.0004931
c(lower = est - 1.96 * se, upper = est + 1.96 * se)
#That gives about 0.00105 to 0.00298. To express it as percent change per 100 m:
(exp(100 * c(0.00105, 0.00298)) - 1) * 100
#which is roughly 11% to 35%. Because you only have 12 species, 
#a t-based interval would be slightly wider, so mention that in your write-up.

X <- model.matrix(~ Elevation, data = FinalMammal)
var_fixed  <- var(as.vector(X %*% pglmmModellog$B))
var_random <- 0.8737 + 0.0000003235
var_resid  <- 0.5042

marginal    <- var_fixed / (var_fixed + var_random + var_resid)
conditional <- (var_fixed + var_random) / (var_fixed + var_random + var_resid)
c(marginal = marginal, conditional = conditional)








