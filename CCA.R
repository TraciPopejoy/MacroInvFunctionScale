#CCA

### to do: 
#take out weird taxa - adult diptera, isotomidae, Hymenoptera, other.other, 
#isopoda - probably terrestrial fall ins

#Field
#Hypothesis: hydrology/temperature drives community structure
#predict: discharge & season will explain more variation than musses
# mussels should be correlated with discharge though
CCAField.data<-F.ComDens %>% left_join(FieldDischarge) %>% 
  left_join(FieldChlA) %>% dplyr::select(-Hymenoptera.miscH,
                                         -Entom.Isotomidae,
                                         -Other.other,
                                         -Isopoda.miscI) %>%
  mutate(Site=substr(Reach, 1,2))
CCAField.data$Mussel.g.m2<-replace_na(CCAField.data$Mussel.g.m2,0)
CCAField.data[CCAField.data$FSamID=="GL-MR.Fall2015","Discharge.cms"]<-0.0294 #usgs
CCAField.data[CCAField.data$FSamID=="GL-NM.Fall2015","Discharge.cms"]<-0.0274 #usgs
CCAField.data[CCAField.data$FSamID=="KT-MR.Fall2015","Discharge.cms"]<-0.109 
CCAField.data[CCAField.data$FSamID=="L3-MR.Summer2015","Discharge.cms"]<-0.819 
CCAField.data[CCAField.data$FSamID=="L3-NM.Summer2015","Discharge.cms"]<-0.556 

names(CCAField.data[,-c(1:7,66:69)])
cca.F<-cca(CCAField.data[,-c(1:7,66:69)]~
             CCAField.data$Mussel.g.m2+CCAField.data$Discharge.cms+
             CCAField.data$Season)
cca.F
ccaF.plot<-plot(cca.F)
cca.F.sum<-summary(cca.F)
#str(cca.F.sum)
sort(cca.F.sum$biplot[,1])
CCAF.impsp<-c(names(sort(cca.F.sum$species[,1],decreasing=T)[1:4]),
              names(sort(cca.F.sum$species[,1], decreasing=F)[1:4]),
              names(sort(cca.F.sum$species[,2],decreasing=T)[1:4]),
              names(sort(cca.F.sum$species[,2], decreasing=F)[1:4]))#how much variation each CCA contributes to constraint variation explained
cca.F.sum$concont
cca.F.sum$constr.chi/cca.F.sum$tot.chi*100

# Enclosure
#Hypothesis: Mussels biomass drives com structure, but little difference between sp
#Predict: MusselBM explains most variation, followed by chlA abun
CCAEnc.data<-E.ComDens %>% left_join(EncDischarge) %>% 
  left_join(EncChlAraw) %>% filter(!is.na(ChlAdensity)) %>%
  filter(Week=="w12") %>% dplyr::select(-Entom.Isotomidae,
                                        -Isopoda.miscI)
CCAEnc.data[!is.na(CCAEnc.data$ChlAdensity) &
              CCAEnc.data$ChlAdensity<=0,"ChlAdensity"]<-0

names(CCAEnc.data[,-c(1:8,55:56)])
cca.E<-cca(CCAEnc.data[,-c(1:8,55:56)]~CCAEnc.data$MusselBiomass.g.m2+
             CCAEnc.data$Spp+CCAEnc.data$Type+
             CCAEnc.data$Discharge.cms+CCAEnc.data$ChlAdensity)
cca.E
ccaE.plot<-plot(cca.E, scaling=1)
cca.E.sum<-summary(cca.E)
#str(cca.F.sum)
sort(cca.E.sum$biplot[,1]) #loadings on CCA1
CCAE.impsp<-c(names(sort(cca.E.sum$species[,1],decreasing=T)[1:4]),
              names(sort(cca.E.sum$species[,1], decreasing=F)[1:4]),
              names(sort(cca.E.sum$species[,2],decreasing=T)[1:4]),
              names(sort(cca.E.sum$species[,2], decreasing=F)[1:4]))#how much variation each CCA contributes to constraint variation explained
cca.E.sum$concont
#porportion explained by constraints
cca.E.sum$constr.chi/cca.E.sum$tot.chi*100

#Shell
#Hypothesis: Mussel spp/type begins to matter due to differencese 
#            in habitat area and chlA abundance
#Predict: ChlA matters most, and is highly correlated with Live/Sham 
CCAS.data<-S.ComDens %>% 
  left_join(EncDischarge[EncDischarge$Week=="w12",]) %>% 
  left_join(ShellChl) %>% dplyr::select(-Dip.Adult,
                                        -Orthoptera,
                                        -Entom.Isotomidae)
names(CCAS.data[,c(1:7,38:41)])
cca.S<-cca(CCAS.data[,-c(1:7,38:41)]~CCAS.data$ShellSpecies+
             CCAS.data$Type+CCAS.data$Discharge.cms+
             CCAS.data$ChlAdensity)
cca.S
ccaS.plot<-plot(cca.S)
cca.S.sum<-summary(cca.S)
sort(cca.S.sum$biplot[,1]) #loadings on CCA1
sort(cca.S.sum$biplot[,2])
CCAS.impsp<-c(names(sort(cca.S.sum$species[,1],decreasing=T)[1:4]),
              names(sort(cca.S.sum$species[,1], decreasing=F)[1:4]),
              names(sort(cca.S.sum$species[,2],decreasing=T)[1:4]),
              names(sort(cca.S.sum$species[,2], decreasing=F)[1:4]))
#how much variation each CCA contributes to constraint variation explained
cca.S.sum$concont
#porportion explained by constraints
cca.S.sum$constr.chi/cca.S.sum$tot.chi*100 

#install.packages("devtools")
#devtools::install_github("gavinsimpson/ggvegan")
library(ggvegan)
# Field Plot
autoplot(cca.F)
tidy.ccaF<-fortify(cca.F)
tidy.ccaFSp<-tidy.ccaF %>% filter(Score=="species",
                                  Label %in% CCAF.impsp) %>%
  left_join(FTaxaTable, by=c('Label'='Taxa'))
tidy.ccaFSi<-tidy.ccaF %>% filter(Score=="sites")
tidy.ccaFcon<-tidy.ccaF %>% filter(Score!="species",
                                   Score!="sites",
                                   Score!="constraints")%>%
  mutate(LabelG=c("MusselBM","Discharge","Season","Fall","Summer"))

names(tidy.ccaF)
FieldCCA<-ggplot()+
  #geom_point(data=tidy.ccaFSi, aes(x=CCA1,y=CCA2))+
  geom_text(data=tidy.ccaFSp, aes(x=CCA1,y=CCA2, 
                                  label=Tlabel),
            size=2.7) +
  geom_segment(data = tidy.ccaFcon[tidy.ccaFcon$Score=="biplot",],
               aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.0, "cm")))+
  geom_text(data = tidy.ccaFcon[tidy.ccaFcon$Score=="biplot" & 
                                tidy.ccaFcon$LabelG!="Season",],
            aes(x=CCA1, y=CCA2, label=LabelG), color = "red")+
  geom_text(data = tidy.ccaFcon[tidy.ccaFcon$Score=="centroids",],
            aes(x=CCA1, y=CCA2, label=LabelG), color = "navy")+
  scale_x_continuous(name="CCA 1 : 70.2% variance")+
  scale_y_continuous(name="CCA 2 : 24.1% variance")

# Enclosure Plot
autoplot(cca.E)
tidy.ccaE<-fortify(cca.E)
tidy.ccaESp<-tidy.ccaE %>% filter(Score=="species",
                                  Label %in% CCAE.impsp) %>%
  left_join(FTaxaTable, by=c('Label'='Taxa'))
tidy.ccaESi<-tidy.ccaE %>% filter(Score=="sites")
tidy.ccaEcon<-tidy.ccaE %>% filter(Score!="species",
                                   Score!="sites",
                                   Score!="constraints")%>%
  mutate(LabelG=c("MusselBM","AMBdom","Control","Live","Discharge",
                  "ChlA","ACTdom","AMBdom","Control","Live","Sham"))

EncCCA<-ggplot()+
  #geom_point(data=tidy.ccaESi, aes(x=CCA1,y=CCA2))+
  geom_text(data=tidy.ccaESp, aes(x=CCA1,y=CCA2, 
                                  label=Tlabel),
            size=2.7) +
  geom_segment(data = tidy.ccaEcon[tidy.ccaEcon$Score=="biplot",],
               aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.0, "cm")))+
  geom_text(data = tidy.ccaEcon[tidy.ccaEcon$Score=="biplot"& 
                                  tidy.ccaEcon$LabelG!="Live"& 
                                  tidy.ccaEcon$LabelG!="AMBdom"& 
                                  tidy.ccaEcon$LabelG!="Control",],
            aes(x=CCA1, y=CCA2, label=LabelG), color = "red")+
  geom_text(data = tidy.ccaEcon[tidy.ccaEcon$Score=="centroids",],
            aes(x=CCA1, y=CCA2, label=LabelG), color = "navy")+
  scale_x_continuous(name="CCA 1 : 30.9% variance")+
  scale_y_continuous(name="CCA 2 : 19.1% variance")

#Shell Plot
autoplot(cca.S)
tidy.ccaS<-fortify(cca.S)
tidy.ccaSSp<-tidy.ccaS %>% filter(Score=="species",
                                  Label %in% CCAS.impsp) %>%
  left_join(FTaxaTable, by=c('Label'='Taxa'))
tidy.ccaSSi<-tidy.ccaS %>% filter(Score=="sites")
tidy.ccaScon<-tidy.ccaS %>% filter(Score!="species",
                                   Score!="sites",
                                   Score!="constraints") %>%
  mutate(LabelG=c("Amblema","Sham","Discharge","ChlA",
                  "Actinonaias","Amblema","Live","Sham"))

names(tidy.ccaS)
shellCCA<-ggplot()+
  #geom_point(data=tidy.ccaSSi, aes(x=CCA1,y=CCA2))+
  geom_text(data=tidy.ccaSSp, aes(x=CCA1,y=CCA2, 
                                  label=Tlabel),
            size=2.7) +
  geom_segment(data = tidy.ccaScon[tidy.ccaScon$Score=="biplot",],
               aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.0, "cm")))+
  geom_text(data = tidy.ccaScon[tidy.ccaScon$Score=="biplot" &
                                tidy.ccaScon$LabelG=="Discharge" |
                                tidy.ccaScon$LabelG=="ChlA",],
            aes(x=CCA1, y=CCA2, label=LabelG), color = "red",
            position=position_dodge(width=.1))+
  geom_text(data = tidy.ccaScon[tidy.ccaScon$Score=="centroids",],
          aes(x=CCA1, y=CCA2, label=LabelG), color = "navy",
          position=position_dodge(width=.1))+
  scale_x_continuous(name="CCA 1 : 47.2% variance")+
  scale_y_continuous(name="CCA 2 : 25.9% variance")

library(cowplot)
plot_grid(FieldCCA, EncCCA, shellCCA, nrow=1, labels="AUTO")
ggsave("./Figures/CCAplots.tiff", width=12, height=4)
