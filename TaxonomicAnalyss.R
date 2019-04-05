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

