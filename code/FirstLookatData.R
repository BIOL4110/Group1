hares<-Pedersen_Elevation_Data
str(hares)
library(ggplot2)
hares |> ggplot(aes(Elevation))+geom_bar()
hares |> ggplot(aes(x=Elevation,y=Cortisol,color=Sex))+geom_point(alpha=1)+geom_jitter(alpha=0.75)+xlab("Elevation (m)")+ylab("Cortisol (ng/mL)")+theme_bw()

#Spearman test. Not significant results. Unable to include losLength as covariate.
cor.test(hares$Elevation,hares$Cortisol,method = 'spearman')

#We should run a Kolmogorov-Smirnow test since we have less than n=50
ks.test(hares$Cortisol,'pnorm')

#The results indicated non-normal data

kruskal.test(Cortisol~Elevation,losLength,data=hares)

kruskal.test(Cortisol~Elevation,data=hares)


##Attempting to transform data to be normal
library(dplyr)
hares |> mutate(logElev=log(Elevation)) |> ggplot(aes(Elevation))+geom_bar()
hares<-mutate(hares,transElev=Elevation^0.5)
ks.test(hares$Loslength,'pnorm')

#Would like to bootstrap the data. Still need to figure it out.
library(boot)
set.seed(0)
func<-lm(Cortisol~Elevation,Loslength,data=hares) {
  d=hares[i,]
  fit=lm(Cortisol~losLength,sex,data=d)
  return(coef(fit))
}
help("boot")
func2<-lm(hares,i){
  df<-hares[i,]
  c(cor(df[,2],df[,3]),
    median(df[,2]),
    mean(df[,1])
    )
}
boot(hares,lm(Cortisol~Elevation,Loslength,Sex,data=hares),R=10000)

boot(data = hares, statistic = func, R = 10000, formula =)
