source("InvDataframes.R")

# Gamma & Alpha Diversity --------

# Field Data
#head(F.ComDens)
ncol(F.ComDens[,-c(1:7,65)]) #gamma diversity?

F.Talpha<-data.frame(F.ComDens[,c(1:7)],
           richness=specnumber(F.ComDens[,-c(1:7, 65)]),
           SimpsonsI=diversity(F.ComDens[,-c(1:7,65)], "simpson")) %>%
  group_by(Reach, SamplingSeason, Treatment, Season, Mussel.g.m2) %>%
  summarize(meanID.npm2=mean(InvDensity.npm2),
            meanR=mean(richness),
            meanSimp=mean(SimpsonsI))
ggplot(F.Talpha, aes(x=Treatment, y=meanSimp))+
  geom_boxplot()+facet_wrap(~Season)

#Enclosure Data
#head(E.ComDens)
ncol(E.ComDens[,-c(1:8)]) #gamma diversity?

E.Talpha<-data.frame(E.ComDens[,c(1:8)],
                     InvDensity.npm2=rowSums(E.ComDens[,-c(1:8)]),
                     richness=specnumber(E.ComDens[,-c(1:8)]),
                     SimpsonsI=diversity(E.ComDens[,-c(1:8)], "simpson"))
ggplot(E.Talpha, aes(x=TreatA, y=SimpsonsI))+
  geom_boxplot()+facet_wrap(~Week)

#Shell Data
#head(E.ComDens)
ncol(S.ComDens[,-c(1:6)]) #gamma diversity?

S.Talpha<-data.frame(S.ComDens[,c(1:7)],
                     InvDensity.npcm2=rowSums(S.ComDens[,-c(1:7)]),
                     richness=specnumber(S.ComDens[,-c(1:7)]),
                     SimpsonsI=diversity(S.ComDens[,-c(1:7)], "simpson"))
ggplot(S.Talpha, aes(x=TreatA, y=SimpsonsI))+
  geom_boxplot()+facet_wrap(~SheType, scales="free_x")

# Beta Diversity ---------

#Field Data
F.Tbeta<-betadiver(F.ComDens[,-c(1:7,65)], method="sor")
(F.Beta.mod <- betadisper(F.Tbeta, paste(F.ComDens$Treatment,F.ComDens$SamplingSeason)))
## Perform test
anova(F.Beta.mod)
## Tukey's Honest Significant Differences
plot(TukeyHSD(F.Beta.mod))
boxplot(F.Beta.mod)

#Enclosure Data
E.Tbeta<-betadiver(E.ComDens[,-c(1:8)], method="sor")
(E.Beta.mod <- betadisper(E.Tbeta, E.ComDens$TreatA))
## Perform test
anova(E.Beta.mod)
## Tukey's Honest Significant Differences
plot(TukeyHSD(E.Beta.mod)) 
boxplot(E.Beta.mod)

#Shell Data
S.Tbeta<-betadiver(S.ComDens[,-c(1:7)], method="sor")
(S.Beta.mod <- betadisper(S.Tbeta, S.ComDens$TreatA))
## Perform test
anova(S.Beta.mod)
## Tukey's Honest Significant Differences
plot(TukeyHSD(S.Beta.mod)) 
boxplot(S.Beta.mod)

# Most Abundant Families -----------
# Field Data
F.ComDens %>% select(-InvDensity.npm2,-Mussel.g.m2,-Other.other) %>%
  gather(Taxa, Abundance, -FSamID,-SamplingSeason,-Season,-Reach,-Treatment) %>% 
  group_by(Taxa) %>%
  summarize(mean.A=mean(Abundance), sd.A=sd(Abundance),sum.A=sum(Abundance)) %>%
  arrange(desc(sum.A))
#results: ChironomidaeL, Hydropsychidae, Caenidae, ElmL, Oligo
# Enclosure Data
E.ComDens %>% select(-Enc2,-Week,-TreatA,-Type,-Spp,-MusselBiomass.g.m2) %>%
  gather(Taxa, Abundance,-TEid,-Enc) %>% 
  group_by(Taxa) %>%
  summarize(mean.A=mean(Abundance), sd.A=sd(Abundance),sum.A=sum(Abundance)) %>%
  arrange(desc(sum.A))
#results: Heptageniidae, ChrionomidaeL, Polycentrop, Damselflies, Caenidae, ElmL
# Field Data
S.ComDens %>% select(-SamID,-Enc2,-ShellSpecies,-TreatA,-Type,-Spp, -SheType) %>%
  gather(Taxa, Abundance) %>% 
  group_by(Taxa) %>%
  summarize(mean.A=mean(Abundance), sd.A=sd(Abundance),sum.A=sum(Abundance)) %>%
  arrange(desc(sum.A))
#results: ChironomidaeL, Polycentropidae, Heptageniidae, ElmL

### So I am interested in:
###  Chironomidae, Heptageniidae, Caenidae, Polycentropidae, Elmidae
### because they are abundant in most of my samples
library(ggsci)
#ChironomidaeL
ggplot(F.ComDens, aes(x=Mussel.g.m2, y=Dip.ChironomidaeL))+
  geom_smooth(method="lm", color="black")+
  geom_point(size=2, aes(color=SamplingSeason))+
  scale_color_futurama()
ggplot(E.ComDens, aes(x=MusselBiomass.g.m2, y=Dip.ChironomidaeL))+
  geom_smooth(method="lm", color="black")+
  geom_point(size=2, aes(color=TreatA))+
  scale_color_futurama()
ggplot(S.ComDens, aes(x=SheType,y=Dip.ChironomidaeL))+
  geom_boxplot()+
  geom_point(size=2, aes(color=TreatA))+
  scale_color_futurama()

#Heptageniidae
ggplot(F.ComDens, aes(x=Mussel.g.m2, y=Eph.Heptageniidae))+
  geom_smooth(method="lm", color="black")+
  geom_point(size=2, aes(color=SamplingSeason))+
  scale_color_futurama()
ggplot(E.ComDens, aes(x=MusselBiomass.g.m2, y=Eph.Heptageniidae))+
  geom_smooth(method="lm", color="black")+
  geom_point(size=2, aes(color=TreatA))+
  scale_color_futurama()
ggplot(S.ComDens, aes(x=SheType,y=Eph.Heptageniidae))+
  geom_boxplot()+
  geom_point(size=2, aes(color=TreatA))+
  scale_color_futurama()

#Caenidae
ggplot(F.ComDens, aes(x=Mussel.g.m2, y=Eph.Trcorythidae))+
  geom_smooth(method="lm", color="black")+
  geom_point(size=2, aes(color=SamplingSeason))+
  scale_color_futurama()
ggplot(E.ComDens, aes(x=MusselBiomass.g.m2, y=Eph.Trcorythidae))+
  geom_smooth(method="lm", color="black")+
  geom_point(size=2, aes(color=TreatA))+
  scale_color_futurama()
ggplot(S.ComDens, aes(x=SheType,y=Eph.Trcorythidae))+
  geom_boxplot()+
  geom_point(size=2, aes(color=TreatA))+
  scale_color_futurama()

#Polycentropidae
ggplot(F.ComDens, aes(x=Mussel.g.m2, y=Tri.Polycentropidae))+
  geom_smooth(method="lm", color="black")+
  geom_point(size=2, aes(color=SamplingSeason))+
  scale_color_futurama()
ggplot(E.ComDens, aes(x=MusselBiomass.g.m2, y=Tri.Polycentropidae))+
  geom_smooth(method="lm", color="black")+
  geom_point(size=2, aes(color=TreatA))+
  scale_color_futurama()
ggplot(S.ComDens, aes(x=SheType, y=Tri.Polycentropidae))+
  geom_boxplot()+
  geom_point(size=2, aes(color=TreatA))+
  scale_color_futurama()

#Elmidae
ggplot(F.ComDens, aes(x=Mussel.g.m2, y=Col.ElmL))+
  geom_smooth(method="lm", color="black")+
  geom_point(size=2, aes(color=SamplingSeason))+
  scale_color_futurama()
ggplot(E.ComDens, aes(x=MusselBiomass.g.m2, y=Col.ElmL))+
  geom_smooth(method="lm", color="black")+
  geom_point(size=2, aes(color=TreatA))+
  scale_color_futurama()
ggplot(S.ComDens, aes(x=SheType,y=Col.ElmL))+
  geom_boxplot()+
  geom_point(size=2, aes(color=TreatA))+
  scale_color_futurama()

# Fuzzy ordination --------------------
#
# Mussel Biomass is already in all of the data sets.
library(fso)
# Discharge ===========================
FieldDischarge<-read_excel("DischargeField.xlsx") %>% 
  group_by(Reach, SamplingSeason) %>%
  summarize(Discharge.cms=mean(`Discharge (M3/s)`, na.rm=T))
####suplementing values find a better way! ###
#based on F16 data from respective reach (or S16 for L3)
FieldDischarge["FSamID"=="GL-MR.Fall2015","Discharge.cms"]<-0.120 
FieldDischarge["FSamID"=="GL-NM.Fall2015","Discharge.cms"]<-0.103 
FieldDischarge["FSamID"=="KT-MR.Fall2015","Discharge.cms"]<-0.109 
FieldDischarge["FSamID"=="L3-MR.Summer2015","Discharge.cms"]<-0.819 
FieldDischarge["FSamID"=="L3-NM.Summer2015","Discharge.cms"]<-0.556 

EncDVw09<-read.csv("../FEn17//FEn17_data/EncPhysDisFEn17OK.csv") %>% 
  rename(Enclosure=`ï..Enclosure`) %>% left_join(TreatENC) %>%
  mutate(Discharge.cms=0.25*V.mps*Depth.m,
         Week="w09") %>%
  dplyr::select(Enc2,Week,Discharge.cms)
EncDVw12<-read_excel("../FEn17/FEn17_data/videowithabiotics.xlsx") %>% 
  left_join(TreatENC, by=c("Unit"="Enclosure")) %>%
  mutate(Discharge.cms=0.25*Velocity*Depth,
         Week="w12") %>% filter(Month=="Oct") %>%
  dplyr::select(Enc2,Week,Discharge.cms)
EncDischarge<-rbind(EncDVw09,EncDVw12)

# Chlorophyll Abundance ===============
FieldChlA<-read_excel("ChlField.xlsx")

# Temperature?? =======================
head(F.BioDens[,-c(1:6)])#those columns just identifiers
F.BioDens$Mussel.g.m2<-F.BioDens$Mussel.g.m2 %>% replace_na(0)
FieldFuzzyData<-F.BioDens %>% left_join(FieldDischarge) %>%
  left_join(FieldChlA)
head(FieldFuzzyData[,-c(1:6,54:56)])
F.Fdist<-dsvdis(FieldFuzzyData[,-c(1:6,54:56)], 'bray')
F.fso<-fso(FieldFuzzyData$Mussel.g.m2, F.Fdist)
summary(F.fso)
F.fso.plot<-FieldFuzzyData %>% mutate(mu=F.fso$mu)
ggplot(F.fso.plot, aes(x=Mussel.g.m2, y=mu))+
  geom_point(aes(color=Season), size=3)+geom_smooth(method="lm")

ordinatedlabels<-LeonRAData[order(LeonRAData$standLRKM),]
unordinatedlabels<-LeonRAData[order(LeonRAData$LRKMlabel),]
greyC<-gray.colors(length(unique(Leongraph$variable)),
                   
                   start = 0.02, end = 0.97, gamma = 2.2, alpha = NULL)
Leonplot<-ggplot(data = Leongraph, aes(x = standLRKM, y = value, fill=var2)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values=greyC, name="Species")+
  scale_x_discrete(labels=rownames(ordinatedlabels))+
  coord_flip()+
  labs(x="Ordinated Sites", y="Relative Abundance") +
  theme(axis.text.x = element_text(size=9,color="black"),
        axis.title.y=element_text(size=20),
        plot.background = element_blank(),
        panel.border=element_blank(),
        panel.grid.major= element_line(colour=NA), 
        panel.grid.minor=element_line(colour=NA),
        title=element_text(size=20),
        panel.background = element_rect(fill = "white"),
        axis.line.x=element_line(colour="black"),
        axis.line.y=element_line(colour="black"))
Leonplot 

##### NOTES:
# so my fuzzy ord worked, mussel biomass can predict community
# need to explore what that acually means in terms of taxonomy
# [ want to make a plot like I did for Con Bio to explore ]
# multivariate fuzzy ord also worked with 
# mussel  biomass, water col chl, disch
# nothing super high, but all significant
# should plot mu mussel biomass and mu chlorophy then add treatment color

F.mfso<-mfso(~FieldFuzzyData$Mussel.g.m2+FieldFuzzyData$Discharge.cms+
               FieldFuzzyData$WC_CHLA_MG.L,
           F.Fdist)
summary(F.mfso)
str(F.mfso)
F.mfso.plot<-data.frame(Reach=FieldFuzzyData$FSamID,
                        MB.mu=F.mfso$mu[,1],
                        Dis.mu=F.mfso$mu[,2],
                        Ch.mu=F.mfso$mu[,3])
ggplot(F.mfso, aes(x=), y=F.mfso$mu[,3])+geom_point()
