#CCA

# to do:
# remove FFG from image
# figure out correct interpretation of % variation
# write results

##### Field #####
#Hypothesis: hydrology/temperature drives community structure
#predict: discharge & season will explain more variation than mussels
# mussels should be correlated with discharge though

# which taxa to remove because rare
#47 samples, 3% samples is 1.4
#so removing taxa only found in 2 samples
F.ComDens[,-c(1:7)] %>% gather(taxa,den) %>% filter(den!=0) %>%
  group_by(taxa) %>% tally() %>% arrange(n)

#only including Fall data to avoid psuedoreplication
CCAField.data<-F.ComDens %>% left_join(Fenv.data) %>% 
  filter(SamplingSeason=="Fall2016") %>%
  dplyr::select(-Hymenoptera.miscH,-Dip.Ascillidae,-Dip.Chaoboridae, #only 1
                -Dip.Thaumaleidae,-Hem.Corixidae, -Hem.Gerridae,-Hydrozoans.miscH,
                -Tri.Odontoceridae, #only in 2 samples
                -Biv.Unionidae, #don't care about these
                -Entom.Isotomidae, -Other.other,-Isopoda.miscI) #terrestrials 
CCAField.data$Mussel.g.m2<-replace_na(CCAField.data$Mussel.g.m2,0)
#View(CCAField.data) #check it matches environment

names(CCAField.data[,-c(1:7,58:67)]) #identifying inverts
#partial cca that takes the site affect out (using HUC12num)

# isolate community data and DONT transform bc cca usines chisquare
field.com<-CCAField.data[,-c(1:7, 58:67)] 
# isolate environmental variables and scale to z scores
field.env<-CCAField.data[,c(1:7, 60:67)] %>% mutate_if(is.numeric,scale)

cca.F<-cca(field.com~Mussel.g.m2+Benthic_CHLA_MG.M2+
             Discharge.cms+D50+Condition(HUC12num), field.env, scaling=2)

cca.F
ccaF.plot<-plot(cca.F)
cca.F.sum<-summary(cca.F)
#str(cca.F.sum)
sort(cca.F.sum$biplot[,1])
sort(cca.F.sum$biplot[,2])
#how much variation each CCA contributes to constraint variation explained
CCAF.impsp<-data.frame(CCA1.high=as.character(names(sort(cca.F.sum$species[,1],decreasing=T)[1:4])),
              CCA1.low=as.character(names(sort(cca.F.sum$species[,1], decreasing=F)[1:4])),
              CCA2.high=as.character(names(sort(cca.F.sum$species[,2],decreasing=T)[1:4])),
              CCA2.low=as.character(names(sort(cca.F.sum$species[,2], decreasing=F)[1:4])))
cca.F.sum$concont
cca.F.sum$constr.chi/cca.F.sum$tot.chi*100

vif.cca(cca.F)
# values over 10 mean redundant measures
cca.F0<-cca(field.com~1+Condition(HUC12num), field.env, scaling=2)
anova(cca.F) #if the global model is significant, can do forward model selection
finalF<-ordistep(cca.F0,scope=formula(cca.F), direction="forward")
finalF$anova$padj <- p.adjust (finalF$anova$`Pr(>F)`, method = 'holm', n = 4)
finalF$anova
###### Enclosure #####
#Hypothesis: Mussels biomass drives com structure, but little difference between sp
#Predict: MusselBM explains most variation, followed by chlA abun
CCAEnc.data<-E.ComDens %>% left_join(Eenv.data)%>% 
  filter(!is.na(ChlAdensity)) %>%
  filter(Week=="w12") %>% dplyr::select(-Entom.Isotomidae,
                                        -Isopoda.miscI) %>%
  mutate(TMusBM=ACT+AMB)
CCAEnc.data[!is.na(CCAEnc.data$ChlAdensity) &
              CCAEnc.data$ChlAdensity<=0,"ChlAdensity"]<-0

names(CCAEnc.data[,-c(1:9,56:61)])
# isolate community data and NOT transform bc cca uses chisquare dist
enc.com<-CCAEnc.data[,-c(1:9, 56:61)] 
# isolate environmental variables and scale to z scores
enc.env<-CCAEnc.data[,c(1:9, 56:61)] %>% ungroup()%>% mutate_if(is.numeric,scale)

##### note: have total mussel biomass but it is redundant with treatment factors
cca.E<-cca(enc.com~ACT+AMB+Type+
             Discharge.cms+ChlAdensity+D50, enc.env, scaling=2)
vif.cca(cca.E)
drop1(cca.E, test="perm")#which is least informative constraint
cca.E
ccaE.plot<-plot(cca.E)
cca.E.sum<-summary(cca.E)
#str(cca.F.sum)
sort(cca.E.sum$biplot[,1]) #loadings on CCA1
#how much variation each CCA contributes to constraint variation explained
CCAE.impsp<-data.frame(CCA1.high=names(sort(cca.E.sum$species[,1],decreasing=T)[1:4]),
                       CCA1.low=names(sort(cca.E.sum$species[,1], decreasing=F)[1:4]),
                       CCA2.high=names(sort(cca.E.sum$species[,2],decreasing=T)[1:4]),
                       CCA2.low=names(sort(cca.E.sum$species[,2], decreasing=F)[1:4]))
cca.E.sum$concont
#porportion explained by constraints
cca.E.sum$constr.chi/cca.E.sum$tot.chi*100

anova(cca.E, perm=1000)

##### Shell #####
#Hypothesis: Mussel spp/type begins to matter due to differencese 
#            in habitat area and chlA abundance
#Predict: ChlA matters most, and is highly correlated with Live/Sham 

###### to do - switch fr
CCAS.data<-S.ComDens %>% 
  left_join(EncDischarge[EncDischarge$Week=="w12",]) %>% 
  left_join(ShellChl) %>% dplyr::select(-Dip.Adult,
                                        -Orthoptera,
                                        -Entom.Isotomidae) %>% 
  left_join(ECounts[,-c(4:12)]) %>%
  filter(!is.na(ChlAdensity))
CCAS.data[!is.na(CCAS.data$ChlAdensity) &
              CCAS.data$ChlAdensity<=0,"ChlAdensity"]<-0
names(CCAS.data[,-c(1:7,38:41)])
# isolate community data and NOT transform bc cca uses chisquare
shell.com<-CCAS.data[,-c(1:7, 38:41)]
# isolate environmental variables and scale to z scores
shell.env<-CCAS.data[,c(1:7, 38:41)] %>% mutate_if(is.numeric,scale)

cca.S<-cca(shell.com~ShellSpecies+Type+
             ChlAdensity, shell.env, scaling=2)
cca.S
ccaS.plot<-plot(cca.S)
cca.S.sum<-summary(cca.S)
sort(cca.S.sum$biplot[,1]) #loadings on CCA1
sort(cca.S.sum$biplot[,2])
CCAS.impsp<-data.frame(CCA1.high=names(sort(cca.S.sum$species[,1],decreasing=T)[1:4]),
              CCA1.low=names(sort(cca.S.sum$species[,1], decreasing=F)[1:4]),
              CCA2.high=names(sort(cca.S.sum$species[,2],decreasing=T)[1:4]),
              CCA2.low=names(sort(cca.S.sum$species[,2], decreasing=F)[1:4]))
#how much variation each CCA contributes to constraint variation explained
cca.S.sum$concont
#porportion explained by constraints
cca.S.sum$constr.chi/cca.S.sum$tot.chi*100 
anova(cca.S, perm=1000)
cca.S0<-cca(shell.com~1, shell.env, scaling=2)
finalS<-ordistep(cca.S0,scope=formula(cca.S), direction="forward")
finalS$anova$padj <- p.adjust (finalS$anova$`Pr(>F)`, method = 'holm', n = 3)
finalS$anova
##### Plots #####
#install.packages("devtools")
#devtools::install_github("gavinsimpson/ggvegan")
library(ggvegan)
# Field Plot
#autoplot(cca.F)
tidy.ccaF<-fortify(cca.F)
tidy.ccaFSp<-tidy.ccaF %>% filter(Score=="species") %>%
  left_join(FTaxaTable, by=c('Label'='Taxa'))
tidy.ccaFspespe<-tidy.ccaFSp %>% filter(Label %in% unique(c(t(CCAF.impsp))))
tidy.ccaFSi<-tidy.ccaF %>% filter(Score=="sites") %>% cbind(field.env[,c(2,4)])
tidy.ccaFcon<-tidy.ccaF %>% filter(Score!="species",
                                   Score!="sites",
                                   Score!="constraints")%>%
  mutate(LabelG=c("MusselBM","Benthic Chl.","Discharge","Substrate"))

names(tidy.ccaF)
FieldCCA<-ggplot()+
  geom_point(data=tidy.ccaFSi, aes(x=CCA1,y=CCA2, #label=Reach, 
                                   shape=Treatment))+
  geom_text(data=tidy.ccaFspespe, aes(x=CCA1,y=CCA2, 
                                  label=Tlabel), size=2.8) +
  geom_text(data = tidy.ccaFcon[tidy.ccaFcon$Score=="biplot" & 
                                  tidy.ccaFcon$LabelG!="Season",],
            aes(x=CCA1, y=CCA2, label=LabelG), color = "red")+
  geom_segment(data = tidy.ccaFcon[tidy.ccaFcon$Score=="biplot",],
               aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.0, "cm")), alpha=.7)+
  scale_x_continuous(name="CCA 1 : 33.9% variance")+
  scale_y_continuous(name="CCA 2 : 11.2% variance") +
  scale_color_brewer(palette = "Dark2", name="FFG") +
  theme(legend.position = "none")

# Enclosure Plot
#autoplot(cca.E)
tidy.ccaE<-fortify(cca.E)
tidy.ccaESp<-tidy.ccaE %>% filter(Score=="species",
                                  Label %in% unique(c(t(CCAE.impsp)))) %>%
  left_join(FTaxaTable, by=c('Label'='Taxa'))
tidy.ccaESi<-tidy.ccaE %>% filter(Score=="sites")
tidy.ccaEcon<-tidy.ccaE %>% filter(Score!="species",
                                   Score!="sites",
                                   Score!="constraints")%>%
  mutate(LabelG=c("ActBM","AmbBM","Live","Sham","Discharge","ChlA",
                  "substrate", "Control","Live","Sham"))
    #"ACTS","AMBL","AMBS","CTRL","Discharge","ChlA","Substrate","ACTL",
    #              "ACTS","AMBL","AMBS","CTRL"))


EncCCA<-ggplot()+
  geom_text(data=tidy.ccaESp, aes(x=CCA1,y=CCA2, 
                                  label=Tlabel), size=2.8) +
  geom_segment(data = tidy.ccaEcon[tidy.ccaEcon$Score=="biplot",],
               aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.0, "cm")))+
  geom_text(data = tidy.ccaEcon[tidy.ccaEcon$Score=="biplot"& 
                                  tidy.ccaEcon$LabelG!="Live"& 
                                  tidy.ccaEcon$LabelG!="Sham",],
            aes(x=CCA1, y=CCA2, label=LabelG), color = "red")+
  geom_text(data = tidy.ccaEcon[tidy.ccaEcon$Score=="centroids",],
            aes(x=CCA1, y=CCA2, label=LabelG), color = "navy")+
  scale_x_continuous(name="CCA 1 : 4.0% variance")+
  scale_y_continuous(name="CCA 2 : 2.8% variance")+
  scale_color_brewer(palette = "Dark2", name="FFG", guide=F)

#Shell Plot
#autoplot(cca.S)
tidy.ccaS<-fortify(cca.S)
tidy.ccaSSp<-tidy.ccaS %>% filter(Score=="species",
                                  Label %in% unique(c(t(CCAS.impsp)))) %>%
  left_join(FTaxaTable, by=c('Label'='Taxa'))
tidy.ccaSSi<-tidy.ccaS %>% filter(Score=="sites")
tidy.ccaScon<-tidy.ccaS %>% filter(Score!="species",
                                   Score!="sites",
                                   Score!="constraints") %>%
  mutate(LabelG=c("Amblema","Sham","ChlA",
                  "Actinonaias","Amblema","Live","Sham"))

shellCCA<-ggplot()+
  #geom_point(data=tidy.ccaSSi, aes(x=CCA1,y=CCA2))+
  geom_text(data=tidy.ccaSSp, aes(x=CCA1,y=CCA2, 
                                  label=Tlabel), size=2.8) +
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
  scale_x_continuous(name="CCA 1 : 8.8% variance")+
  scale_y_continuous(name="CCA 2 : 4.9% variance")+
  scale_color_brewer(palette = "Dark2", name="FFG", guide=F)

library(cowplot)
ffgplot<-plot_grid(FieldCCA, EncCCA, shellCCA, nrow=1, labels="AUTO")
legend1<-get_legend(FieldCCA+
                      theme(legend.position = c(.1,.9),
                            legend.direction = "horizontal"))
plot_grid(ffgplot,legend1, nrow=2, rel_heights = c(1,.1))
ggsave("./Figures/CCAplotsG.tiff", width=10, height=4)



head(CCAField.data)
CCAField.data.2<-CCAField.data[complete.cases(CCAField.data),]
field.com<-CCAField.data.2[,-c(1:7, 59:68)]
field.env<-CCAField.data.2[,c(1:7, 59:68)]
Fm1<-cca(field.com~Mussel.g.m2+Year+
           WC_CHLA_MG.L+Benthic_CHLA_MG.M2+ Discharge.cms+D10+D50+
           Condition(HUC12num+Season+Year), field.env[,-c(1:3,7,8,13,17)])
Fm0<-cca(field.com~1, field.env[,-c(1:3,7,8,13,17)])

sink("ccatesting.txt")
field.m<-ordistep(Fm0, scope=formula(Fm1), method="back")
field.m
sink()
vif.cca(field.m)
plot(field.m, display=c("sp","bp"))
field.m$anova


