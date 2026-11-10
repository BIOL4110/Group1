

##Read the CSV file, and load packages
HorBase <- read.csv("data/raw/mammalraw.csv")
library(dplyr)
library(ggplot2)
install.packages("WebPower")
library(WebPower)
library(lmtest)

##Sort Mammal Cortisol Data
MammalData <- HorBase
str(MammalData)
print(MammalData)
class(MammalData$F_BC_Mean)
MammalData<-HorBase |> filter(Vert_Group=="Mammal",CORT=="cortisol",!Genus%in%c("Balaenoptera","Trichechus","Phoca","Halichoerus","Spermophilus"))  |> select(Species, F_BC_Mean,M_BC_Mean,Common_name,F_BC_N,M_BC_N,Genus,Elevation) |>  filter_out(F_BC_Mean=="Null",M_BC_Mean=="Null")                   

MammalData

write.csv(MammalData, "mammalfinal.csv", row.names = FALSE)


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

FinalMammal$F_BC_Mean<-as.numeric(FinalMammal$F_BC_Mean)

cort_by_elevation<-FinalMammal |> group_by(Elevation) |> summarise(Cortisol=mean(Combined_BC_Mean))

ggplot(cort_by_elevation)+geom_point(aes(x=Elevation,y=Cortisol))+geom_smooth(aes(x=Elevation,y=Cortisol),method=lm,se=FALSE)+theme_bw()

###Take out null species

MammalData <- read.csv("data/mammalfinal.csv")
MammalData <- MammalData[-c(10:13), ]
write.csv(MammalData, "data/mammalfinal.csv", row.names = FALSE)
View(MammalData)

MammalData <- read.csv("data/mammalfinal.csv")

MammalData <- read.csv("data/mammalfinal.csv", na.strings = "null")

MammalData <- MammalData |>
  mutate(across(c(F_BC_N, M_BC_N, Elevation), ~ replace(., is.na(.), 0)))

write.csv(MammalData, "data/mammalfinal.csv", row.names = FALSE, na = "null")



###Statistical Tests

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

#Linear regression power analysis
wp.regression(n = 57, p1 = 17, f2 = 0.48197, alpha = 0.05, power = NULL)






