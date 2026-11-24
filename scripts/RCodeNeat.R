##Packages
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