#### Invertebrates on at Enclosure Level ####
library(readxl);library(tidyverse)
#bring in treatment data from enclosures
TreatENC<-read_excel("../FEn17/FEn17_data/FEn17OKTreatments.xlsx") 
#bring in the invert data from basket sampling
FEn17Inv12<-read.csv("../FEn17/FEn17_data/FEn17InvMeas.csv", stringsAsFactors = F) %>%
  filter(Label!="rulerxocc.tif") %>% select(-Mean, -Max, -Min) %>%
  mutate(TEid=substring(Label, 6,11),
         Enc=substring(Label, 9,11),
         Week=substring(Label, 6,8))
FEn17Inv09<-read.csv("../FEn17/FEn17_data/FEn17w09.csv", stringsAsFactors = F) %>% 
  filter(Label!="rulerxocc.tif") %>% select(-Mean, -Max, -Min, - SampleSplit)%>%
  mutate(Week="w09",
         TEid=paste(Week, Enc, sep=""))
EnclInv<-rbind(FEn17Inv12, FEn17Inv09) #contains every insect identified from baskets

unique(EnclInv$Taxa)[order(unique(EnclInv$Taxa))] #check taxa for correct spelling

#### Calculate area sampled ####
BaskSampled<-read.csv("../FEn17/FEn17_data/SubSamFEn17OK.csv", stringsAsFactors = F) %>%
  mutate(WeekG=case_when(Week==8~"w09",
                         Week==4~"w04",
                         Week==12 ~"w12")) %>% 
  left_join(TreatENC, by="Enclosure") %>% 
  filter(Enclosure!="WC" & Type.x=="Insect") %>%
  select(Enclosure, Enc2, Basket., WeekG) %>%
  mutate(TEid=paste(WeekG, Enc2, sep="")) %>% filter(!duplicated(TEid))
  
#Some filters didn't have baskets associated with them
#assumed 3 baskets
BaskSampled[is.na(BaskSampled$Basket.),"Basket."]<-3

#macroinvertebrate summary in enclosures
ECounts<-EnclInv %>% group_by(TEid, Enc, Week) %>% 
  summarize(Nsam=n(),
            richness=length(unique(Taxa))) %>% left_join(BaskSampled) %>%
  mutate(AreaSamp=(Basket.*.03315),
         InvertDensity.npm2=Nsam/AreaSamp) %>% 
  full_join(TreatENC) %>% 
  select(-Enclosure, -Enc2, -Basket., -WeekG)
View(ECounts %>% group_by(Type, TreatA) %>% summarize_if(is.numeric, mean))
View(ECounts %>% group_by(Week) %>% summarize_if(is.numeric, mean))
ECounts <-ECounts %>% 
  mutate(Enc2=substr(TEid, 4,6), Week=substr(TEid, 1,3)) %>%
  left_join(MusBiomass) %>% replace(is.na(.),0) %>%
  mutate(MusselBiomass.g.m2=sumBM/0.25)
summary(lmer(InvertDensity.npm2~MusselBiomass.g.m2+Week+(1|Enc2), data=ECounts))

#long data frame of Enc, Sp, Taxa, Count
EnclComMat<-EnclInv %>% group_by(TEid,Enc, Week, Taxa) %>% 
  summarize(N=n()) %>%
  left_join(ECounts) %>% mutate(IndDensity.m2=N/AreaSamp)%>%
  select(Enc, Week,TEid, Taxa, IndDensity.m2) %>% 
  mutate(id=1:n()) %>% filter(Taxa!="MYSTERY") %>%
  spread(Taxa,IndDensity.m2) %>% group_by(TEid) %>%
  summarise_if(is.numeric,funs(mean(.[which(!is.na(.))]))) %>%
  select(-id) 

#need to make it proportional abundance
EnclPCom<-EnclInv %>% group_by(TEid,Enc, Week, Taxa) %>% 
  filter(Taxa!="MYSTERY", Taxa!="MysteryTricoptera", Taxa!="Col.Myst",
         Taxa!="Tri.Pupa") %>% 
  summarize(N=n()) %>%  mutate(id=1:n()) %>% 
  spread(Taxa,N) %>% group_by(TEid)%>%
  summarise_if(is.numeric,funs(mean(.[which(!is.na(.))]))) %>%
  select(-id) %>% replace(is.na(.),0) %>%
  mutate(rowsum=rowSums(.[,-1])) %>% group_by(TEid) %>%
  summarise_if(is.numeric,funs(./rowsum)) %>%
  select(-rowsum) %>%
  select(TEid, everything())

FTaxaTable<-read_excel("Macroinv Power Law Coeffs TBP.XLSX", sheet=2)
match(names(EnclComMat[,-1]), FTaxaTable$Taxa) #checking to make sure all taxa in taxa table

#so we have a species community matrix ShellComMat
#need a trait table of those species
EnclTraits<-FTaxaTable[match(names(EnclComMat[,-1]), FTaxaTable$Taxa),]
EnclTraits$Taxa==names(EnclComMat[,-1]) #need it to all be true
EnclTraits <- EnclTraits %>%
  filter(T.TropP!=0) %>%
  select(Taxa,T.Habit, T.TropP, T.TropS,
         T.dev,T.Lifespan,T.crwl,T.swim,T.MatSize)


#### biomass data ####
mreg<-read_excel("../FEn17/FEn17_data/LENGTH-MASS_CLA_CCV_20161208-reg coefficients.xlsx",
                 sheet=2)
MusselData %>% left_join(TreatENC) %>%
  mutate(BiomassE=case_when(Genus=="AMB"~ 1.57e-6*L^3.21,
                            Genus=='ACT'~ 4.7e-7*L^3.51))%>%
  group_by(Genus) %>%
  summarize(meanBM=mean(BiomassE, na.rm=T))
MusBiomass<-MusselData %>% left_join(TreatENC) %>%
  mutate(BiomassE=case_when(Notes=="NotRecovered"&Genus=="AMB"~3.47,
                            Notes=="NotRecovered"&Genus=="ACT"~7.64,
                            Genus=="AMB"~ 1.57e-6*L^3.21,
                            Genus=='ACT'~ 4.7e-7*L^3.51)) %>%
  group_by(Enc2, Treatment) %>%
  summarize(sumBM=sum(BiomassE, na.rm=T), 
            meanBM=mean(BiomassE,na.rm=T), 
            sdBM=sd(BiomassE, na.rm=T))

#### exploring family abundance ####
head(EnclComMat)
View(EnclPCom %>%
       gather(Taxa, Abundance, -TEid) %>% group_by(Taxa) %>%
       summarize(mean.A=mean(Abundance, na.rm=T)*100))

View(EnclComMat %>%
       gather(Taxa, Abundance, -TEid) %>% group_by(Taxa) %>%
       summarize(mean.A=mean(Abundance, na.rm=T), sum.A=sum(Abundance, na.rm=T)))

EnclGraph<-EnclComMat %>% 
  mutate(Enc2=substr(TEid, 4,6), Week=substr(TEid, 1,3)) %>%
  left_join(TreatENC) %>% left_join(MusBiomass) %>%
  mutate(MusselBiomass.g.m2=sumBM/.25) %>% replace(is.na(.),0)

#Heptageniidae
ggplot(EnclGraph, aes(x=MusselBiomass.g.m2,y=Eph.Heptageniidae))+
  geom_point(size=2, aes(color=TreatA))+geom_smooth(method="lm")
summary(lm(Eph.Heptageniidae~MusselBiomass.g.m2, data=EnclGraph))
#ChironomidaeL - FFG minx
ggplot(EnclGraph, aes(x=MusselBiomass.g.m2,y=Dip.ChironomidaeL))+
  geom_point(size=2, aes(color=TreatA))+geom_smooth(method="lm")
summary(lm(Dip.ChironomidaeL~MusselBiomass.g.m2, data=EnclGraph))
#Elmidae larvae
ggplot(EnclGraph, aes(x=MusselBiomass.g.m2,y=Col.ElmL))+
  geom_point(size=2, aes(color=TreatA))+geom_smooth(method="lm")
summary(lm(Col.ElmL~MusselBiomass.g.m2, data=EnclGraph))
#Polycentropidae
ggplot(EnclGraph, aes(x=MusselBiomass.g.m2,y=Tri.Polycentropidae))+
  geom_point(size=2, aes(color=TreatA))+geom_smooth(method="lm")
summary(lm(Tri.Polycentropidae~MusselBiomass.g.m2, data=EnclGraph))

ggplot(EnclGraph, aes(x=log1p(Dip.ChironomidaeL)))+geom_histogram()

log1p.EnclComMat<-EnclComMat %>% #top 20 taxa
  select(Eph.Heptageniidae, Dip.ChironomidaeL,Col.ElmL,Tri.Polycentropidae,
         Od.Damselflies,Eph.Trcorythidae,Eph.Leptophlebiidae, 
         Dip.ChironomidaeP,Col.Psephenidae,
         Col.ElmA, Eph.Polymitarcyidae,Gas.Ancylidae,Oligochaeta.misc,
         Biv.Corbicula,Hydrozoans.miscH,Od.Dragonflies, Megaloptera.miscM,
         Copepoda.Cyclopoida, Hirudinea.leach, Eph.Myst,-TEid) %>% 
  replace(is.na(.),0) %>%
  mutate_if(is.numeric, .funs=log1p)
encl.pca<-prcomp(log1p.EnclComMat)
plot(encl.pca, type="l")
summary(encl.pca)
print(encl.pca)
encl.pca$rotation[,2]

EnclGraph<-EnclGraph %>% mutate(PCA1=encl.pca$x[,1],
                                  PCA2=encl.pca$x[,2])
ggplot(EnclGraph, aes(x=PCA1, y=PCA2))+
  geom_point(aes(color=TreatA, shape=Week))+
  scale_color_futurama()
ggplot(EnclGraph, aes(x=MusselBiomass.g.m2, y=PCA1))+
  geom_point(aes(color=TreatA, shape=Week))+geom_smooth(method="lm")
PCAlineE<-data.frame(species=rownames(encl.pca$rotation),Scores=encl.pca$rotation[,2])
ggplot(PCAlineE)+
  geom_text(aes(-1,Scores,label=round(Scores,3)))+
  geom_text(aes(1,Scores, label=species))+
  geom_vline(xintercept=0)+
  scale_y_continuous(trans="log")+
  coord_cartesian(xlim=c(-3,3))+
  labs(x=NULL,y=NULL)+theme_classic()
