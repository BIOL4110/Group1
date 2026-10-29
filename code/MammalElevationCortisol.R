HorBase<-read.csv("HormoneBase_v1.csv")
library(dplyr)
hobase<-HorBase |> filter(Genus=="Lepomis") |> select(F_A_Mean,M_A_Mean,Years,Year_1,Year_final,Species,Elevation)
str(hobase)
library(ggplot2)
library(tidyverse)
str(hobase)
tidydata<-hobase |> pivot_longer(cols=c("F_A_Mean","M_A_Mean"),names_to = "A_Mean_type",values_to="values")
ggplot(tidydata, aes(x=Year_final,y=values),group=Year_final)+geom_point()

##Mammal Cortisol Data
str(MammalData)
print(MammalData)
class(MammalData$F_BC_Mean)
MammalData<-HorBase |> filter(Vert_Group=="Mammal",CORT=="cortisol",!Genus%in%c("Balaenoptera","Trichechus","Phoca","Halichoerus","Spermophilus"))  |> select(Species, F_BC_Mean,M_BC_Mean,Common_name,F_BC_N,M_BC_N,Genus,Elevation)                    

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

library(ggplot2)

FinalMammal$Elevation<-as.numeric(FinalMammal$Elevation)
FinalMammal$Combined_BC_Mean<-as.numeric(FinalMammal$Combined_BC_Mean)

cort_by_elevation<-FinalMammal |> group_by(Elevation) |> summarise(Cortisol=mean(Combined_BC_Mean))

ggplot(cort_by_elevation)+geom_point(aes(x=Elevation,y=Cortisol))+geom_smooth(aes(x=Elevation,y=Cortisol),method=lm,se=FALSE)+theme_bw()

