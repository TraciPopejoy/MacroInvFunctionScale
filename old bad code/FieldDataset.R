#libraries
library(tidyverse); library(ggplot2); library(readxl); library(vegan)

#############     Invert Count in Field Surber Samples    ##############
Count2016<-read_excel("2016 Field Data.xlsx", sheet = "Counts") %>%
  mutate(Reach=paste(Site, Treatment, sep="-")) %>% 
  group_by(Sample, Location, Site, Treatment, Reach, SamplingSeason) %>% 
  summarize_if(is.numeric, sum, na.rm=T)
Count2015<-read_excel("2015 Data Hotspot Field Surber.xlsx") %>% 
  mutate(Reach=paste(Site, Treatment, sep="-")) %>%
  group_by(Sample, Reach, Treatment, SamplingSeason, Site, Location) %>%
  summarize_if(is.numeric, sum, na.rm=T)
CountField<-rbind(Count2016, Count2015)

CountField$InvDensity.npm2<-rowSums(CountField[,7:68], na.rm=T)/0.092903
CountField[is.na(CountField)]<-0
CountField$richness<-specnumber(CountField[,7:68])

fieldSiteKey<-CountField %>% group_by(Sample, Reach, SamplingSeason, Treatment) %>% 
  tally() %>% select(-n) %>% mutate(FSamID=paste(Reach,SamplingSeason, sep="."))
surbernum<-fieldSiteKey %>% group_by(Reach, SamplingSeason) %>% tally()

FieldComMat<-CountField %>% group_by(Reach, SamplingSeason) %>% 
  summarize_if(is.numeric, mean) %>% group_by(Reach, SamplingSeason) %>%
  summarise_if(is.numeric, funs(./0.092903)) %>%
  mutate(FSamID=paste(Reach, SamplingSeason, sep=".")) %>% ungroup()%>%
  select(-Reach, -Other.other, -richness)%>%
  mutate(Treatment=substr(FSamID, 4,5)) %>%
  select(FSamID, Treatment,SamplingSeason,everything())

FieldPCom<-CountField %>% group_by(Reach, SamplingSeason) %>% 
  summarize_if(is.numeric, mean) %>%
  mutate(FSamID=paste(Reach, SamplingSeason, sep="."))%>% ungroup()%>%
  select(-InvDensity.npm2,-Other.other,-Tri.Pupa, -richness,-SamplingSeason,-Reach) %>%
  gather(Taxa,Count, -FSamID) %>% group_by(FSamID, Taxa) %>% 
  filter(Taxa!="MYSTERY", Taxa!="MysteryTricoptera", Taxa!="Col.Myst",
         Taxa!="Tri.Pupa", Taxa!="Dip.Adult", Taxa!="Dip.Cyclorrhaphous",
         Taxa!="Dip.Other", Taxa!="Dip.Orthorrhaphous", Taxa!="Dip.Thaumaleidae",
         Taxa!="Hymenoptera.miscH", Taxa!="Dip.Ascillidae") %>%
  mutate(id=1:n()) %>% 
  spread(Taxa,Count) %>%  select(-id) %>% ungroup() %>%
  mutate(rowsum=rowSums(.[,-1])) %>% group_by(FSamID) %>%
  summarise_if(is.numeric,funs(./rowsum))%>%
  select(-rowsum) %>%   select(FSamID, everything())

FTaxaTable<-read_excel("Macroinv Power Law Coeffs TBP.XLSX", sheet=2)
length(match(names(FieldPCom[,-1]), FTaxaTable$Taxa)) #checking to make sure all taxa in taxa table

#so we have a species community matrix ShellComMat
#need a trait table of those species
FieldTraits<-FTaxaTable[match(names(FieldPCom[,-1]), FTaxaTable$Taxa),]
FieldTraits$Taxa==names(FieldPCom[,-1]) #need it to all be true
FieldTraits <- FieldTraits %>%
  filter(T.TropP!=0) %>%
  select(Taxa,T.Habit, T.TropP, T.TropS,
         T.dev,T.Lifespan,T.crwl,T.swim,T.MatSize)

FieldPComSUR<-CountField %>% group_by(Reach, SamplingSeason) %>% 
  #mutate(FSamID=paste(Reach, SamplingSeason, sep="."))%>% 
  ungroup()%>%
  select(-InvDensity.npm2,-Other.other,-Tri.Pupa, -richness,-SamplingSeason,-Reach,
         -Location,-Site,-Treatment)%>%
  gather(Taxa,Count, -Sample) %>% 
  filter(Taxa!="MYSTERY", Taxa!="MysteryTricoptera", Taxa!="Col.Myst",
         Taxa!="Tri.Pupa", Taxa!="Dip.Adult", Taxa!="Dip.Cyclorrhaphous",
         Taxa!="Dip.Other", Taxa!="Dip.Orthorrhaphous", Taxa!="Dip.Thaumaleidae",
         Taxa!="Hymenoptera.miscH")%>%
  mutate(id=1:n()) %>%  spread(Taxa,Count) %>% group_by(Sample) %>%
  summarise_all(funs(mean(.[which(!is.na(.))]))) %>% select(-id) %>% ungroup()%>%
  mutate(rowsum=rowSums(.[,-1])) %>% group_by(Sample) %>%
  summarise_if(is.numeric,funs(./rowsum))%>%
  select(-rowsum) %>%  select(Sample, everything()) %>%
  filter(Sample!="KT-MR Oct16 S1")

#### mussel biomass estimate ####
FieldBMest<-read.csv("dw_L-W_est_by individualOLS-USETHIS.csv") %>%  
  group_by(Site) %>%
  summarize(Mussel.bm.g=sum(sppsp_ols_mass..g.),
            Quadrat.n=max(Quadrat)) %>%
  mutate(Mussel.g.m2=Mussel.bm.g/(Quadrat.n*0.25),
         Treatment=case_when(substr(Site, 3,4)=="NM"~"NM",
                             T~"MR"),
         Reach=paste(substr(Site,1,2),Treatment, sep="-"))

#### exploring family abundance ####
head(FieldComMat)
FieldPCom %>% 
       gather(Taxa, Abundance, -FSamID) %>% group_by(Taxa) %>%
       summarize(mean.A=mean(Abundance)*100)
FieldComMat %>% select(-InvDensity.npm2) %>%
  gather(Taxa, Abundance, -FSamID) %>% group_by(Taxa) %>%
  summarize(mean.A=mean(Abundance), sum.A=sum(Abundance))

FieldGraph<-F.ComDens %>% left_join(fieldSiteKey[,-1]) %>%
  filter(!duplicated(.)) %>% left_join(FieldBMest)

#ChironomidaeL - FFG minx
ggplot(FieldGraph, aes(x=Mussel.g.m2,y=Dip.ChironomidaeL))+
  geom_smooth(method="lm", color="black")+
  geom_point(size=2, aes(color=SamplingSeason))
#Hydropsychidae Filter Feeder
ggplot(FieldGraph, aes(x=Mussel.g.m2, y=Tri.Hydropsychidae))+
  geom_smooth(method="lm", color="black") +
  geom_point(size=2, aes(color=SamplingSeason))
#Caenidae - C-Gather
ggplot(FieldGraph, aes(x=Mussel.g.m2,y=Eph.Trcorythidae))+
  geom_smooth(method="lm", color="black")+
  geom_point(size=2, aes(color=SamplingSeason))
#Elmidae larvae
ggplot(FieldGraph, aes(x=Mussel.g.m2,y=Col.ElmL))+
  geom_smooth(method="lm", color="black")+
  geom_point(size=2, aes(color=SamplingSeason))
#Oligochaeta
ggplot(FieldGraph, aes(x=Mussel.g.m2,y=Oligochaeta.misc))+
  geom_point(size=2, aes(color=SamplingSeason))+geom_smooth(method="lm")
#Heptageniidae - Grazers
ggplot(FieldGraph, aes(x=Mussel.g.m2,y=Eph.Heptageniidae))+
  geom_point(size=2, aes(color=SamplingSeason))+geom_smooth(method="lm")

ggplot(FieldGraph, aes(x=log1p(Dip.ChironomidaeL)))+geom_histogram()
ggplot(FieldGraph, aes(x=log1p(Tri.Hydropsychidae)))+geom_histogram()

log1p.FieldCom<-F.ComDens %>% 
  select(-FSamID, -Tri.Bracy, -Dip.Chaoboridae, -`FLAT WORMS`, -Entom.Isotomidae,
         -Dip.Culicidae,-Tri.Odontoceridae,-Neu.Sisyridae,-Col.Gyrinidae,
         -Dip.Tupulidae,-Hem.Gerridae, -Col.Hydrophilidae, -Hem.Corixidae,
         -Hydrozoans.miscH, -Lep.Pyralidae, -Isopoda.miscI, -Eph.Myst,
         -Amphipoda.miscA, -Biv.Unionidae,-Dip.Similidae,-Hirudinea.leach,
         -Ostracoda,-Col.Psephenidae, -Tri.Helicophyche,-Tri.Uenoidae,-Copepoda.miscC,
         -Tri.Molannidae,-InvDensity.npm2, -Dip.Ascillidae) %>% #removing taxa <0.1%
  mutate_if(is.numeric, .funs=log1p)
log1p.FieldCom20<-F.ComDens %>% #20 most abundant taxa
  select(Dip.ChironomidaeL, Tri.Hydropsychidae, Eph.Trcorythidae, Col.ElmL,
         Oligochaeta.misc, Eph.Heptageniidae, Biv.Corbicula, Dip.ChironomidaeP,
         Tri.Philopteridae, Tri.Hydroptilidae, Cladocera.miscD, Eph.Polymitarcyidae,
         Eph.Baetidae, Eph.Leptophlebiidae, Tromb.Hydracarina, Tri.Leptoceridae,
         Col.ElmA, Od.Damselflies, Plecoptera.Perlidae, Tri.Polycentropidae) %>% 
  mutate_if(is.numeric, .funs=log1p)
field.pca<-prcomp(log1p.FieldCom20)
field.pcaL<-prcomp(log1p.FieldCom[,6:40])
plot(field.pcaL, type="l")
summary(field.pcaL)
print(field.pcaL)
field.pcaL$rotation[,1]

FieldGraph<-F.ComDens[,1:7] %>% mutate(PCA1=field.pcaL$x[,1],
                                  PCA2=field.pcaL$x[,2]) %>%
  left_join(FieldChlA)
ggplot(FieldGraph, aes(x=WC_CHLA_MG.L, y=PCA1))+
  geom_point(aes(color=Treatment, shape=Season),size=2)+
  scale_color_aaas()+theme_krementz()+
  ylab("PCA1 25% var")+xlab("Water Column ChlA mg/L")

PCAlineF<-data.frame(species=rownames(field.pcaL$rotation),
                     Scores=field.pcaL$rotation[,2])
ggplot(PCAlineF)+
  geom_text(aes(-1,Scores,label=round(Scores,3)))+
  geom_text(aes(1,Scores, label=species))+
  geom_vline(xintercept=0)+
  #scale_y_continuous(trans="log")+
  coord_cartesian(xlim=c(-3,3))+
  labs(x=NULL,y=NULL)+theme_classic()


library(ggsci)
ggplot(FieldGraph, aes(x=PCA1, y=PCA2))+
  geom_point(aes(color=Treatment, shape=Season),size=2)+
  scale_color_aaas()+theme_krementz()+
  xlab("PCA1 25% variance")+ylab("PCA2 14% variance")

### describe invert density ###
FCounts<-CountField %>% ungroup() %>%
  select(Reach, SamplingSeason, InvDensity.npm2, richness, Treatment) %>%
  group_by(Reach, SamplingSeason, Treatment) %>%
  summarize_if(is.numeric, mean, na.rm=T) %>%
  mutate(SeasF=factor(SamplingSeason, 
                      levels=c("Summer2015","Fall2015","Summer2016","Fall2016"))) %>%
  left_join(FieldBMest) %>% replace(is.na(.), 0)
library(lme4);library(lmerTest)
summary(lmer(log1p(InvDensity.npm2)~Mussel.g.m2+SeasF+(1|Reach), data=FCounts))
View(FCounts %>% group_by(Treatment) %>%
  select(-Reach, -SamplingSeason, -Site) %>%summarize_all(mean))

#field sum
FieldSum<-CountField %>% ungroup() %>% 
  summarize_if(is.numeric,sum, na.rm=T) %>% select(-InvDensity.npm2, -richness) %>%
  gather() %>% 
  arrange(desc(value)) %>% mutate(rank=1:n())

ggplot(FieldSum, aes(x=rank, y=value))+geom_bar(stat="identity")+
  xlab("Abundance Rank")+ylab("Total Count")

FieldSumFFG<-CountField %>% ungroup() %>%
  select(-Sample,-Location, -Site,-Treatment,-Reach,-SamplingSeason,
         -richness,-InvDensity.npm2) %>%
  gather(Taxa, value) %>% group_by(Taxa) %>%
  summarize_if(is.numeric, sum, na.rm=T) %>%left_join(FTaxaTable) %>%
  select(Taxa, value, Order, T.TropP) %>%
  mutate(FFG=case_when(T.TropP==1~"C-Gatherer",
                       T.TropP==2~"C-Filterer",
                       T.TropP==3~"Herbivore",
                       T.TropP==4~"Predator",
                       T.TropP==5~"Shredder")) %>%
  arrange(desc(value)) %>% mutate(rank=1:n())
ggplot(FieldSumFFG, aes(x=rank, y=value))+
  geom_bar(stat="identity", aes(fill=FFG))+
  scale_y_continuous(trans="sqrt", breaks=c(100,2000,4000,8000,12000,14000))+
  xlab("Abundance Rank")+ylab("Total Count")

########## OLD CODE #########
#############     Field Invert Biomass Calculation     #############
BiomassReg<-read_excel("Macroinv Power Law Coeffs TBP.xlsx", sheet = 1)
TaxaList<-read_excel("Macroinv Power Law Coeffs TBP.xlsx", sheet = 2)

biomass<-function(Fam, Ord, Length) {
  #tmass<-NULL
  #if not ID'd to family, use order level regressions
  if(is.na(match(Fam, BiomassReg$Family))){ 
    plcoeffs<-BiomassReg[BiomassReg$Order == Ord &
                           !is.na(BiomassReg$Order == Ord),]
  }else{ #pull all the regressions for that family
    plcoeffs<-BiomassReg[BiomassReg$Family==Fam&
                           BiomassReg$Order == Ord&
                           !is.na(BiomassReg$Family==Fam), ] }
  #which regressions were actually built for insects this size
  idx2<- c(Length>=plcoeffs[,19] & Length<=plcoeffs[,20]) 
  idx2[is.na(idx2)]<-T #if no size range listed, use the regression anyways
  d1<-plcoeffs[idx2,] #row by row????
  mass<-d1$a*(Length^d1$b) #power law to determine biomass
  #print(paste(Fam, Ord, sep=" - "))
  mean(mass, na.rm=T)
}

#### 2016 biomass invertebrate calculations ####
Lengths2016<-read_excel("2016 Field Data.xlsx", sheet= 2,
                        col_types = c("text","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",#5
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",#10
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric")) 
(Mtax16<-colnames(Lengths2016)[2:63])

#this takes about two minutes because 17200 observations
melt2016lengths<-Lengths2016 %>% group_by(Sample) %>% 
  gather(Taxa, Length.cm, -Sample)%>% 
  drop_na() %>% left_join(TaxaList) %>% filter(Order!="misc")%>% 
  rowwise()%>% 
  mutate(Length.mm=Length.cm*10,
         BM.mg=mean(biomass(Family, Order, Length.mm))) 

samp2016count<-Count2016 %>% ungroup() %>% select(-Location, -Site, -Treatment, 
                                                  -totalinverts, -richness,-Reach, 
                                                  -Other.other) %>%
  gather(Taxa, Count, -Sample, -SamplingSeason) %>% filter(Count!=0) %>%
  filter(!is.na(Count))

SxT2016Long<-melt2016lengths %>% group_by(Sample,Taxa) %>% 
  summarize(maxBM=max(BM.mg), 
            minBM=min(BM.mg), 
            meanBM=mean(BM.mg, na.rm=T),
            SDBM=sd(BM.mg, na.rm=T)) %>% 
  left_join(samp2016count) %>%
  mutate(tbm.mg.avg=meanBM*Count)
SxT2016Long %>% filter(is.na(meanBM)) #likely samples that are to small/large

BM2016<-SxT2016Long %>% group_by(Sample) %>% 
  summarize(Total.BM.mg=sum(tbm.mg.avg, na.rm=T)) %>% 
  left_join(fieldSiteKey)
CT2016<-Count2016 %>% select(Sample, Reach, SamplingSeason, Site, Location, 
                             totalinverts,richness)
InvSummary2016<-full_join(BM2016, CT2016)%>% #there are some NA, likely spelling mistakes
  group_by(Reach, SamplingSeason, Site, Location, Treatment) %>%
  mutate(BMdens.g.m2=Total.BM.mg*0.001/0.092903,
         invdens.m2=totalinverts/0.092903) %>%
  summarize_if(is.numeric,funs(mean,sd), na.rm=T)
