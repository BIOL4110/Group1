hares<-Pedersen_Elevation_Data
str(hares)
library(ggplot2)
hares |> ggplot(aes(Elevation))+geom_bar()
hares |> ggplot(aes(x=Elevation,y=Cortisol,color=Sex))+geom_point(alpha=1)+geom_jitter(alpha=0.75)+xlab("Elevation (m)")+ylab("Cortisol (ng/mL)")+theme_bw()
cor.test(hares$Elevation,hares$Cortisol,method = 'spearman')
#We should run a Kolmogorov-Smirnow test since we have less than n=50
ks.test(hares$Cortisol,'pnorm')

#The results indicated non-normal data

kruskal.test(Cortisol~Elevation,losLength,data=hares)

kruskal.test(Cortisol~Elevation,data=hares)

library(dplyr)
hares |> mutate(logElev=log(Elevation)) |> ggplot(aes(Elevation))+geom_bar()
hares<-mutate(hares,transElev=Elevation^0.5)
ks.test(hares$Loslength,'pnorm')

library(boot)
func<-lm(Cortisol~Elevation,Loslength,data=hares)
