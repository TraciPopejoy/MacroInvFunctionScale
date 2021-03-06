#libraries
library(tidyverse)

### feeding groups table from LeRoy Poff 2008 ###
ffg<-read_csv("feedingtraits.csv")
FTaxaTable<-read_excel("Macroinv Power Law Coeffs TBP.XLSX", sheet=2)
taxa<-sort(unique(FTaxaTable$Family[FTaxaTable$Family!="misc"]))
FFGunk<-ffg %>% right_join(FTaxaTable[FTaxaTable$Family!="misc",]) %>%
  select(Family, Order, Feed_mode_prim, Feed_mode_sec, Feed_mode_comments) %>%
  distinct()
traitDIST<-FFGunk %>% group_by(Family) %>% filter(!is.na(Feed_mode_sec)) %>%
  summarize(filt=sum(str_count(Feed_mode_sec, "filterer"), na.rm=T),
            gath=sum(str_count(Feed_mode_sec, "gatherer"), na.rm=T),
            other=sum(str_count(Feed_mode_sec, "Other"), na.rm=T),
            parasite=sum(str_count(Feed_mode_sec, "Parasite"), na.rm=T),
            pierc=sum(str_count(Feed_mode_sec, "Piercer"), na.rm=T),
            pred=sum(str_count(Feed_mode_sec, "Predator"), na.rm=T),
            scrap=sum(str_count(Feed_mode_sec, "Scraper"), na.rm=T),
            shred=sum(str_count(Feed_mode_sec, "Shredder"), na.rm=T))
#use this table to supplement other data sources 
#(mainly Le Roy Poff 2006 table)

##### MACROINVERTEBRATE DATA FROM OTHER SCRIPTS #####
source("FieldDataset.R")#Field Data 
head(FieldPCom)
head(FieldTraits)

source("EnclosureDataset.R")#Enclosure Data
head(EnclPCom)
head(EnclTraits)

#Shell Data
head(ShellPCom)
head(ShellTraits)

##### Structural Difference Between Treatments #####
library(vegan)
## alpha diversity ##
EnAlp<-data.frame(Enc=EnclComMat$Enc,TreatA=EnclComMat$TreatA,alpha=specnumber(EnclComMat))
boxplot(alpha~TreatA,EnAlp)

## beta diversity ##
SorEncl<-betadiver(EnclComMat[50:99,-c(1:3)], "sor")
(EnclbetaMod <- betadisper(SorEncl, EnclComMat$TreatA[50:99]))
anova(EnclbetaMod) #difference between groups significant?
permutest(EnclbetaMod, pairwise = TRUE, permutations = 99)
## Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(EnclbetaMod))

## gamma diversity ##
gammaEncl<-data.frame(TreatA=c("AMBL","AMBS","ACTL","ACTS","CTRL"),
            tot.spp=c(sum(colSums(EnclComMat[EnclComMat$TreatA=="AMBL",-c(1:3)])>0),
                      sum(colSums(EnclComMat[EnclComMat$TreatA=="AMBS",-c(1:3)])>0),
                      sum(colSums(EnclComMat[EnclComMat$TreatA=="ACTL",-c(1:3)])>0),
                      sum(colSums(EnclComMat[EnclComMat$TreatA=="ACTS",-c(1:3)])>0),
                      sum(colSums(EnclComMat[EnclComMat$TreatA=="CTRL",-c(1:3)])>0)))

ran.G.E<-round(rnorm(999,mean(gammaEncl$tot.spp), sd(gammaEncl$tot.spp)),0)
hist(ran.G.E)
abline(v=gammaEncl$tot.spp, lty=2, lwd=4, 
       col=c("blue","purple","red","orange","yellow"))
gammaEncl$Zscore<-(gammaEncl$tot.spp - mean(ran.G.E)) / sd(ran.G.E)
gammaEncl

gammaField<-data.frame(Treatment=c("MR F2015","MR F2015","MR F2015","MR F2015",
                                   "NM F2015","NM S2015","NM F2016","NM S2016"),
               tot.spp=c(sum(colSums(FieldComMat[FieldComMat$Treatment=="MR" &
                                                   FieldComMat$SamplingSeason=="Fall2015",
                                                 -c(1:3)])>0),
                         sum(colSums(FieldComMat[FieldComMat$Treatment=="MR" &
                                                   FieldComMat$SamplingSeason=="Summer2015",
                                                 -c(1:3)])>0),
                         sum(colSums(FieldComMat[FieldComMat$Treatment=="MR" &
                                                   FieldComMat$SamplingSeason=="Fall2016",
                                                 -c(1:3)])>0),
                         sum(colSums(FieldComMat[FieldComMat$Treatment=="MR" &
                                                   FieldComMat$SamplingSeason=="Summer2016",
                                                 -c(1:3)])>0),
                         sum(colSums(FieldComMat[FieldComMat$Treatment=="NM" &
                                                   FieldComMat$SamplingSeason=="Fall2015",
                                                 -c(1:3)])>0),
                         sum(colSums(FieldComMat[FieldComMat$Treatment=="NM" &
                                                   FieldComMat$SamplingSeason=="Summer2015",
                                                 -c(1:3)])>0),
                         sum(colSums(FieldComMat[FieldComMat$Treatment=="NM" &
                                                   FieldComMat$SamplingSeason=="Fall2016",
                                                 -c(1:3)])>0),
                         sum(colSums(FieldComMat[FieldComMat$Treatment=="NM" &
                                                   FieldComMat$SamplingSeason=="Summer2016",
                                                 -c(1:3)])>0)))
ran.G.F<-round(rnorm(999,mean(gammaField$tot.spp), sd(gammaField$tot.spp)),0)
hist(ran.G.F)
abline(v=gammaField$tot.spp, lty=2, lwd=4, 
       col=c("blue","purple","red","orange","yellow","goldenrod4","green","forestgreen"))
gammaEncl$Zscore<-(gammaEncl$tot.spp - mean(ran.G.E)) / sd(ran.G.E)
gammaEncl
                   
##### AlphaFD.R from DecomposingFD #####
source("../DecomposingFD/R/AlphaFD.R")
#Field Samples
max(dist(FieldTraits[,-1]))
FTdist<-dist(FieldTraits[,-1])/max(dist(FieldTraits[,-1]))
FTD(FTdist)
#surbers stand as own samples
FAlphaSUR<-FTD.comm(as.matrix(FTdist), as.matrix(FieldPComSUR[,-1]))
FSitesAlphaSUR<-cbind(FAlphaSUR[[1]], FieldPComSUR[,1]) %>% left_join(fieldSiteKey)%>%
  filter(!duplicated(.))
ggplot(FSitesAlphaSUR, aes(x=SamplingSeason, y=qDT, fill=Treatment))+
  geom_boxplot(aes(group=interaction(SamplingSeason,Treatment)))+
  ylab("qDT - Functional Species Diversity")+
  theme(legend.direction="horizontal",legend.position = "bottom")
ggplot(FSitesAlphaSUR, aes(x=SamplingSeason, y=M, fill=Treatment))+
  geom_boxplot(aes(group=interaction(SamplingSeason,Treatment)))+
  ylab("M - magnitude of functional dispersion")+
  theme(legend.direction="horizontal",legend.position = "bottom")
#aggregated surber samples
FAlphaF<-FTD.comm(as.matrix(FTdist), as.matrix(FieldPCom[,-1]))
FSitesAlpha<-cbind(FAlphaF[[1]], FieldPCom[,1]) %>% left_join(fieldSiteKey[,-1])%>%
  filter(!duplicated(.)) %>% left_join(FieldBMest, by="Reach") %>%
  mutate(SamSeaF=factor(SamplingSeason, 
                        levels=c("Summer2015", "Fall2015", "Summer2016", "Fall2016")))
ggplot(FSitesAlpha, aes(x=Mussel.g.m2, y=qDTM, color=SamSeaF))+
  geom_point()+
  ylab("qDTM\nintegrates Fn diversity & dispersion")+
  geom_smooth(aes(group=SamplingSeason), method="lm", level=0, color="black")+
  facet_wrap(~SamSeaF)+theme_bw()
summary(lm(qDTM~Treatment.x+Mussel.g.m2+SamSeaF, data=FSitesAlpha))
ggplot(FSitesAlpha, aes(x=SamSeaF, y=qDTM, fill=Treatment.x))+
  geom_boxplot(aes(group=interaction(SamplingSeason,Treatment.x)))+
  ylab("qDTM\nintegrates Fn diversity & dispersion")+
  theme(legend.direction="horizontal",legend.position = "bottom")
ggplot(FSitesAlpha, aes(x=SamplingSeason, y=qDT, fill=Treatment.x))+
  geom_boxplot(aes(group=interaction(SamplingSeason,Treatment.x)))+
  ylab("qDT - Functional Species Diversity")+
  theme(legend.direction="horizontal",legend.position = "bottom")
ggplot(FSitesAlpha, aes(x=SamplingSeason, y=M, fill=Treatment.x))+
  geom_boxplot(aes(group=interaction(SamplingSeason,Treatment.x)))+
  ylab("M - magnitude of functional dispersion")+
  theme(legend.direction="horizontal",legend.position = "bottom")


#Enclosure Samples
max(dist(EnclTraits[,-1]))
ETdist<-dist(EnclTraits[,-1])/max(dist(EnclTraits[,-1]))
FTD(ETdist)

EnAlphaF<-FTD.comm(as.matrix(ETdist), as.matrix(EnclPCom[,-1]))
EnSitesAlpha<-cbind(EnAlphaF[[1]], EnclPCom[,1]) %>% 
  mutate(Enc2=substr(TEid,4,6),Week=substr(TEid,1,3)) %>%
  left_join(TreatENC) %>% left_join(MusBiomass) %>%
  replace(is.na(.), 0) %>% mutate(MusselDens.g.m2=sumBM/0.25)
ggplot(EnSitesAlpha, aes(x=MusselDens.g.m2, y=qDTM))+
  geom_point(size=2, aes(color=TreatA))+
  #geom_smooth(method="lm", aes(group=TreatA))+
  geom_smooth(method="lm", aes(group=Type, linetype=Type), 
              level=0, color="black")+
  ylab("qDTM\nintegrates Fn diversity & dispersion")+
  xlab("Mussel Biomass g m-2")+
  scale_color_manual(values=CP[c(1,2,5,4,3)])+
  theme_bw()+facet_wrap(~Week)
ggplot(EnSitesAlpha, aes(x=Week, y=qDTM, fill=TreatA))+
  geom_boxplot(aes(group=interaction(Week,TreatA)))+
  ylab("qDTM\nintegrates Fn diversity & dispersion")+
  scale_fill_manual(values=CP[c(1,2,5,4,3)])
ggplot(EnSitesAlpha, aes(x=Week, y=qDT, fill=TreatA))+
  geom_boxplot(aes(group=interaction(Week,TreatA)))+
  ylab("qDT - Functional Species Diversity")+
  scale_fill_manual(values=CP[c(1,2,5,4,3)])
ggplot(EnSitesAlpha, aes(x=Week, y=M, fill=TreatA))+
  geom_boxplot(aes(group=interaction(Week,TreatA)))+
  ylab("M - magnitude of functional dispersion")+
  scale_fill_manual(values=CP[c(1,2,5,4,3)])
 
#Shell Samples
max(dist(ShellTraits[,-1]))
STdist<-dist(ShellTraits[,-1])/max(dist(ShellTraits[,-1]))
FTD(STdist)

ShAlphaF<-FTD.comm(as.matrix(STdist), as.matrix(ShellPCom[,-1]))
ShSitesAlpha<-cbind(ShAlphaF[[1]], ShellPCom[,1]) %>% 
  mutate(Enc2=substr(SamID, 1,3), Shell=substr(SamID,5,7)) %>% left_join(TreatENC) %>%
  mutate(shelTreat=paste(Type,Shell,sep="-"),
         TF=factor(TreatA, levels=c("ACTL","AMBL","ACTS","AMBS", "CTRL"))) %>%
  left_join(MusBioShell) %>%
  mutate(MusselBiomass.g.m2=sumBM/.25)
ggplot(ShSitesAlpha, aes(x=MusselBiomass.g.m2, y=qDTM))+
  geom_point(aes(color=shelTreat), size=2)+
  ylab("qDTM\nintegrates Fn diversity & dispersion")+
  xlab("Mussel Biomass g m-2")+
  geom_smooth(method="lm", aes(group=Type, linetype=Type),
              level=0, color="black")+
  scale_color_manual(values=CG[c(1,5,2,4)])
shMusBioMod<-lm(qDTM~MusselBiomass.g.m2+Type, data=ShSitesAlpha)
ggplot(ShSitesAlpha, aes(x=TF, y=qDTM, fill=shelTreat))+
  geom_boxplot(aes(group=interaction(TreatA,Shell)))+
  ylab("qDTM\nintegrates Fn diversity & dispersion")+
  xlab("Enclosure Treatment")+
  scale_fill_manual(values=CG[1:4])
ggplot(ShSitesAlpha, aes(x=TF, y=qDT, fill=shelTreat))+
  geom_boxplot(aes(group=interaction(TreatA,Shell)))+
  ylab("qDT - Functional Species Diversity")+
  xlab("Enclosure Treatment")+
  scale_fill_manual(values=CG[1:4])
ggplot(ShSitesAlpha, aes(x=TF, y=M, fill=shelTreat))+
  geom_boxplot(aes(group=interaction(TreatA,Shell)))+
  ylab("M - magnitude of functional dispersion")+
  xlab("Enclosure Treatment")+
  scale_fill_manual(values=CG[1:4])

##### BetaFD.R from DecomposingFD #####
source("../DecomposingFD/R/BetaFD.R")
FBetaF<-FTD.beta(as.matrix(FTdist), as.matrix(FieldPCom[,-1]), abund = T)
str(FBetaF)
EBetaF<-FTD.beta(as.matrix(ETdist), as.matrix(EnclPCom[,-1]), abund = T)
str(EBetaF)
SBetaF<-FTD.beta(as.matrix(STdist), as.matrix(ShellPCom[,-1]), abund = T)
str(SBetaF)
# "functional differences among sites are small
# functionally distince communities do not occur along the gradient"

##### GammaFD.R from DecomposingFD #####
source("../DecomposingFD/R/GammaFD.R")
FGam<-FTD.gamma.str(as.matrix(FTdist), as.matrix(FieldPCom[,-1]), abund = T)
EGam<-FTD.gamma.str(as.matrix(ETdist), as.matrix(EnclPCom[,-1]), abund = T)
SGam<-FTD.gamma.str(as.matrix(STdist), as.matrix(ShellPCom[,-1]), abund = T)

#### Models ####
library(lmerTest);library(emmeans)
Emod<-lmer(qDTM~TreatA*Week+(1|Enc2), data=EnSitesAlpha)
anova(Emod)
summary(Emod)

Smod<-lm(qDTM~TreatA*Shell, data=ShSitesAlpha)
anova(Smod)
summary(Smod)

Fmod<-lmer(qEt~Treatment*SamplingSeason+(1|Reach), data=FSitesAlpha)
anova(Fmod)
summary(Fmod)

#### NMDS of Species Present ####

FieldPCom2<-FieldPCom %>% left_join(fieldSiteKey[,-1]) %>%
  filter(!duplicated(.))%>%
  rename(sample.id=FSamID, Time=SamplingSeason)%>%
  select(-Reach) %>%
  select(sample.id,Time,Treatment, everything())
EnclPCom2<-EnclPCom %>% mutate(Enc2=substr(TEid,4,6),Week=substr(TEid,1,3)) %>%
  left_join(TreatENC) %>%
  rename(sample.id=TEid, Time=Week, Treatment=TreatA)%>%
  select(-TreatF,-Enc2,-Enclosure,-Type,-Spp) %>%
  select(sample.id,Time,Treatment, everything())
ShellPCom2<-ShellPCom %>% mutate(Enc2=substr(SamID,1,3),Week="w12") %>%
  left_join(TreatENC) %>%
  rename(sample.id=SamID, Time=Week, Treatment=TreatA)%>%
  select(-TreatF,-Enc2,-Enclosure,-Type,-Spp) %>%
  select(sample.id,Time,Treatment, everything())

MasterInvPCom<-bind_rows(FieldPCom2, EnclPCom2, ShellPCom2) %>%
  replace(is.na(.),0)
traitDis<-FTaxaTable %>% replace(is.na(.),0) %>%
  select(-Trait.Source,-EPA.Traits,-LowestID,-Macro.meio,-Life.Stage) %>%
  filter(T.TropP!=0) %>%
  select(Taxa,Order,T.Habit, T.TropP, T.TropS,
       T.dev,T.Lifespan,T.crwl,T.swim,T.MatSize)
nmdsM<-metaMDS(traitDis[,-c(1:2)])
plot(nmdsM)
data.scores<-as.data.frame(scores(nmdsM))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$taxa<-as.vector(traitDis$Taxa)
data.scores$order<-as.vector(traitDis$Order)
data.scores$FFG<-as.factor(traitDis$T.TropP)
head(data.scores)  #look at the data
#Using the scores function from vegan to
#extract the species scores and convert to a data.frame
species.scores <- as.data.frame(scores(nmdsM, "species"))  
species.scores$traits <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data

ggplot() + 
  #geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,color=FFG), alpha=0.5, size=4) +  # add the species labels
  geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=taxa), size=3) +  # add the species labels
  coord_equal() +
  #scale_color_discrete(labels=c("C-Gatherer","C-Filterer","Herbivore","Predator",
  #                            "Shredder","Parasite"))+
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())

#### Graphs ####
library(colorspace);library(scales)
CP<-diverge_hcl(5, h=c(180,70), c = 100, l = c(50, 90), power = 1)
CP[3]<-"lightgrey"
show_col(CP)

CG<-diverge_hcl(5, h=c(280,180), c = 100, l = c(50, 90), power = 1)
CG[3]<-"lightgrey"
show_col(CG)


head(FCounts)
ECounts<-ECounts %>% 
  mutate(TreFac=factor(TreatA, levels=c("CTRL","ACTL","ACTS","AMBL","AMBS")))
head(ECounts)
SCounts<-SCounts %>%
  filter(!is.na(ShellSpecies)) %>%
  mutate(ShID=paste(ShellSpecies,Type, sep="x"))
head(SCounts)

fdens<-ggplot(FCounts, aes(x=SeasF, y=InvDensity.npm2, fill=Treatment))+
  geom_boxplot(aes(group=interaction(SamplingSeason,Treatment)))+
  scale_fill_manual(name="River Reach\nTreatment",
                    labels=c("Mussel Bed","Control"),
                    values = CP[c(1,3)])+
  scale_y_continuous(name=expression("Macroinvertebrates per m " ^ -2), trans="log1p",
                     breaks = c(1,500, 750,1000,2500,5000,10000,15000))+
  scale_x_discrete(name="Sampling Period", 
                   labels=c("Su 2015", "F 2015", "Su 2016", "F 2016"))+
  theme(legend.justification=c(1,1),legend.position = c(0.27,1),#legend.position = c(0.575,1)
        legend.text=element_text(size=10))
edenss<-ggplot(ECounts, aes(x=Week, y=InvertDensity.npm2, fill=TreFac))+
  geom_boxplot(aes(group=interaction(Week,TreFac)))+
  scale_fill_manual(name="Enclosure\nTreatment",
                    values = CP[c(3,1,2,5,4)])+
  scale_y_continuous(name=expression("Macroinvertebrates per m " ^ -2))+
  scale_x_discrete(name="Sampling Period", 
                   labels=c("Sept. 2017","Oct. 2017"))+
  theme(legend.justification=c(1,1),legend.position = c(0.27,1),#legend.position = c(0.5,1)
        legend.text=element_text(size=10))
sdens<-ggplot(SCounts, aes(x=TreatA, y=InvertDensity.npcm2/0.0001, fill=ShID))+
  geom_boxplot(aes(group=interaction(TreatA,ShID))) +
  scale_fill_manual(name="Shells",
                    values = CP[c(1,2,5,4)],
                    labels=c("Live ACT", "Sham ACT","Live AMB","Sham AMB"))+
  scale_y_continuous(name=expression("Macroinvertebrates per m " ^ -2))+
  scale_x_discrete(name="Enclosure Treatment")+
  theme(legend.justification = c(1,1),legend.position = c(1,1),
        legend.direction="horizontal", legend.text=element_text(size=10))+
  guides(fill=guide_legend(nrow=2))
library(cowplot)
plot_grid(sdens,edenss,fdens, ncol=3,labels = "AUTO")
ggsave("InvertDensity.tiff", dpi=600, width = 10, height = 3.5, units="in")

#### discharge investigation ####
head(FSitesAlpha)
FieldDischarge<-read_excel("DischargeField.xlsx") %>% 
  group_by(Reach, SamplingSeason) %>%
  summarize(Discharge.cms=mean(`Discharge (M3/s)`, na.rm=T))
FSitesDis<-left_join(FSitesAlpha,FieldDischarge)
ggplot(FSitesDis,
       aes(x=Discharge.cms, y=qDTM, color=Treatment.x, shape=SamplingSeason))+
  geom_point()
summary(lm(qDTM~Mussel.g.m2*Discharge.cms, data=FSitesDis))

EncDVw09<-read.csv("../FEn17//FEn17_data/EncPhysDisFEn17OK.csv") %>% 
  rename(Enclosure=`ï..Enclosure`) %>% left_join(TreatENC) %>%
  mutate(Discharge.cms=0.25*V.mps*Depth.m,
         Week="w09") %>%
  select(Enc2,Week,Discharge.cms)
EncDVw12<-read_excel("../FEn17/FEn17_data/videowithabiotics.xlsx") %>% 
  left_join(TreatENC, by=c("Unit"="Enclosure")) %>%
  mutate(Discharge.cms=0.25*Velocity*Depth,
         Week="w12") %>% filter(Month=="Oct") %>%
  select(Enc2,Week,Discharge.cms)
EncSitesDis<-rbind(EncDVw09,EncDVw12) %>% right_join(EnSitesAlpha)
ggplot(EncSitesDis,
       aes(x=Discharge.cms, y=qDTM, color=TreatA, shape=Week))+
  geom_point(size=2)
summary(lm(qDTM~Discharge.cms, data=EncSitesDis))
summary(lm(qDTM~MusselDens.g.m2*Discharge.cms, data=EncSitesDis))

#### debugging code ####
for(k in 1:nrow(FieldPCom)){
  td<-FTaxaTable[match(names(FieldPCom[k,-1]), FTaxaTable$Taxa),] %>%
    filter(T.TropP!=0) %>%
    select(Taxa,T.Habit, T.TropP, T.TropS,
           T.dev,T.Lifespan,T.crwl,T.swim,T.MatSize)
  tdmat<-dist(td[,-1])/max(dist(td[,-1]))
  spmat<-FieldPCom[k,-1]
  testk<-FTD.comm(as.matrix(tdmat), as.matrix(spmat[,-1]))
  print(k)
}
