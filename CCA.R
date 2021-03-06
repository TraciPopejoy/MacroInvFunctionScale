# CCA analysis =================================
source('InvDataframes.R')
source('EnvDataframes.R')

##### Field Data Analysis ----------
#Hypothesis: hydrology/temperature drives community structure
#predict: discharge & season will explain more variation than mussels
# mussels should be correlated with discharge though

# which taxa to remove because rare
#14 samples, 3% samples is 0.42
#so removing taxa only found in 1 samples

# what are my least abundant taxa
tax_rare_field<-F.ComDens[,-c(1:7)] %>% 
  filter(SampSeaF=="Fall2016") %>% #only looking at 2016
  gather(taxa,den) %>% 
  filter(den!=0) %>%
  group_by(taxa) %>% 
  tally() %>% 
  arrange(n) %>%
  filter(n<=1)
empty_col_field<-names(which(colSums(F.ComDens16[,-c(1:7, 70:71)])==0))


#only including Fall data to avoid psuedoreplication
CCAField.data<-F.ComDens %>% 
  left_join(Fenv.data, by=c("Reach", "SamplingSeason")) %>% 
  filter(SamplingSeason=="Fall2016") %>%
  #removing rare taxa
  dplyr::select(-all_of(empty_col_field), #taxa with no observations (empty col)
                -all_of(tax_rare_field$taxa), # taxa only found in 1 site
                -Biv.Unionidae, #don't care about these
                -Entom.Isotomidae, -Other.other,-Isopoda.miscI, -Dip.Other) #terrestrials 
#replacing na in mussel weight with 0
CCAField.data$`Average of STDM (g.m-2)`<-replace_na(CCAField.data$`Average of STDM (g.m-2)`,0) 
#View(CCAField.data) #check it matches environment

names(CCAField.data[,-c(1:7, 53:62)]) #identifying inverts

# isolate community data and DONT transform bc cca uses chisquare
field.com<-CCAField.data[,-c(1:7, 53:62)] 
# isolate environmental variables and scale to z scores
field.env<-CCAField.data[,c(1:7, 53:62)] %>%
  dplyr::select(-Year)

#run the cca with all interesting environmental variables
cca.F<-cca(field.com~`Average of STDM (g.m-2)`+Benthic_CHLA_MG.M2+
             Discharge.cms+Dvar+Condition(HUC12num), field.env, scaling=2)
cca.F #see the results
vif.cca(cca.F) # values over 10 mean redundant environmental variable
ccaF.plot<-plot(cca.F) # default plot
cca.F.sum<-summary(cca.F) #look at the summary of cca results
#how much variation each CCA contributes to constraint variation explained
sort(cca.F.sum$biplot[,1])
sort(cca.F.sum$biplot[,2])  
#what are the species that had the 3 highest and lowest loadings
CCAF.impsp<-data.frame(CCA1.high=as.character(names(sort(cca.F.sum$species[,1],decreasing=T)[1:2])),
              CCA1.low=as.character(names(sort(cca.F.sum$species[,1], decreasing=F)[1:2])),
              CCA2.high=as.character(names(sort(cca.F.sum$species[,2],decreasing=T)[1:2])),
              CCA2.low=as.character(names(sort(cca.F.sum$species[,2], decreasing=F)[1:2])))
cca.F.sum$constr.chi/cca.F.sum$tot.chi*100 # percentage constrained variation

cca.F0<-cca(field.com~1+Condition(HUC12num), field.env, scaling=2) #set a null model for ordistep
anova(cca.F) #if the global model is significant, can do forward model selection
finalF<-ordistep(cca.F0,scope=formula(cca.F), direction="both") #complete forward and backward model selection
finalF$anova$padj <- p.adjust (finalF$anova$`Pr(>F)`, method = 'holm', n = 4) #determine holm adjusted p value for multiple tests
finalF$anova #resulting best model with p values

###### Enclosure Data Analysis -------------------------
#Hypothesis: Mussels biomass drives com structure, but little difference between sp
#Predict: MusselBM explains most variation, followed by chlA abun

# which taxa to remove because rare
#50 samples, 4% samples is 2
#so removing taxa only found in 2 samples
#what are my rare species
tax_rare_enc<-E.ComDens[,-c(1,2,4:8)] %>% filter(Week=="w12") %>%
  gather(taxa,den) %>% filter(den!=0) %>%
  group_by(taxa) %>% tally() %>% arrange(n) %>%
  filter(n<=2)
empty_col_enc<-names(which(colSums(E.ComDens12[,-c(1:8)])==0))

#build dataset with species and environment data while removing rare species
CCAEnc.data<-E.ComDens %>% left_join(Eenv.data)%>% 
  left_join(MusBioG)%>%
  dplyr::select(-all_of(empty_col_enc), #empty columns
                -all_of(tax_rare_enc$taxa), #rare taxa
                -Entom.Isotomidae,-Isopoda.miscI, -Dip.Other, #terrestrials 
                -Treatment) %>% 
  filter(!is.na(ChlAdensity)) %>% #remove instances were I don't have chlorophyll a
  filter(Week=="w12") #only considering week 12


CCAEnc.data[is.na(CCAEnc.data$ACT),44:45]<-0 #control treatments have 0 mussel biomass
CCAEnc.data <- CCAEnc.data %>%  mutate(TMusBM=ACT+AMB) #calculating total mussel biomass
CCAEnc.data[!is.na(CCAEnc.data$ChlAdensity) &
              CCAEnc.data$ChlAdensity<=0,"ChlAdensity"]<-0 #replacing negative chla values with 0

names(CCAEnc.data[,-c(1:8,39:46)]) #names of invertebrates
# isolate community data and NOT transform bc cca uses chisquare dist
enc.com<-CCAEnc.data[,-c(1:8,39:46)] 
# isolate environmental variables and scale to z scores
enc.env<-CCAEnc.data[,c(1:8,39:46)] %>% ungroup()%>%
  #left_join(EnclosureRaster[,c(1,3,4)], by=c("Enc2"="enc"))%>% 
  mutate_if(is.numeric,scale) %>% 
  #changing type into binary to reduce variable redundancy
  mutate(Live=case_when(Type=="Control"~"No",
                        Type=="Sham"~"No",
                        Type=="Live"~"Yes")) 

# run the cca
cca.E<-cca(enc.com~ACT+AMB+Live+
             Discharge.cms+ChlAdensity+Dvar, enc.env, scaling=2)
#not including conditional spatial vector because only explains 3% variance
vif.cca(cca.E) #check for redundant environmental variables
cca.E
ccaE.plot<-plot(cca.E)
cca.E.sum<-summary(cca.E)
#how much variation each CCA contributes to constraint variation explained
sort(cca.E.sum$biplot[,1]) 
#identify important (high or low loadings) species
CCAE.impsp<-data.frame(CCA1.high=names(sort(cca.E.sum$species[,1],decreasing=T)[1:2]),
                       CCA1.low=names(sort(cca.E.sum$species[,1], decreasing=F)[1:2]),
                       CCA2.high=names(sort(cca.E.sum$species[,2],decreasing=T)[1:2]),
                       CCA2.low=names(sort(cca.E.sum$species[,2], decreasing=F)[1:2]))
cca.E.sum$constr.chi/cca.E.sum$tot.chi*100 #percentage contrained inertia

anova(cca.E, perm=1000) #run the anova-like permutations; its not significant

### checking is benthic chlorophyll was similair among treatments
ggplot(CCAEnc.data) + geom_boxplot(aes(x=TreatA, y=ChlAdensity))
ggplot(CCAEnc.data) + geom_boxplot(aes(x=TreatA, y=Dvar))
ggplot(CCAEnc.data) + geom_boxplot(aes(x=TreatA, y=Discharge.cms))


##### Shell Data Analysis --------------------------------
#Hypothesis: Mussel spp/type begins to matter due to differencese 
#            in habitat area and chlA abundance
#Predict: ChlA matters most, and is highly correlated with Live/Sham 

#5% of species is 1.6
# what are my least abundant species
S.ComDens[,-c(1:7)] %>%
  gather(taxa,den) %>% filter(den!=0) %>%
  group_by(taxa) %>% tally() %>% arrange(n)
  
#build dataset without rare species and with environmental data
CCAS.data<-S.ComDens %>% 
  left_join(EncDischarge[EncDischarge$Week=="w12",]) %>% 
  left_join(ShellChl) %>% dplyr::select(-Dip.Adult,
                                        -Orthoptera,
                                        -Entom.Isotomidae) %>% 
  filter(!is.na(ChlAdensity)) #remove observations with missing chlorophyll data
#replacing negative chlorophyll with 0
CCAS.data[!is.na(CCAS.data$ChlAdensity) &
              CCAS.data$ChlAdensity<=0,"ChlAdensity"]<-0

names(CCAS.data[,-c(1:7,38:41)]) #invertebrate taxa names
# isolate community data and NOT transform bc cca uses chisquare
shell.com<-CCAS.data[,-c(1:7, 38:41)]
# isolate environmental variables and scale to z scores
shell.env<-CCAS.data[,c(1:7, 38:41)] %>% mutate_if(is.numeric,scale)

#run cca for the shell data
cca.S<-cca(shell.com~ShellSpecies+Type+
             ChlAdensity, shell.env, scaling=2)
cca.S
ccaS.plot<-plot(cca.S)
cca.S.sum<-summary(cca.S)
sort(cca.S.sum$biplot[,1]) #loadings on CCA1
sort(cca.S.sum$biplot[,2])
#identify high and low loading species on CCA1 and CCA2
CCAS.impsp<-data.frame(CCA1.high=names(sort(cca.S.sum$species[,1],decreasing=T)[1:2]),
              CCA1.low=names(sort(cca.S.sum$species[,1], decreasing=F)[1:2]),
              CCA2.high=names(sort(cca.S.sum$species[,2],decreasing=T)[1:2]),
              CCA2.low=names(sort(cca.S.sum$species[,2], decreasing=F)[1:2]))
#how much variation each CCA contributes to constraint variation explained
cca.S.sum$concont
cca.S.sum$constr.chi/cca.S.sum$tot.chi*100 #porportion explained by constraints

anova(cca.S, perm=1000) #run the anova to look at significance
cca.S0<-cca(shell.com~1, shell.env, scaling=2) #null model for ordistep
#complete forward and backward selection for shell anova
finalS<-ordistep(cca.S0,scope=formula(cca.S), direction="both")
#holm adjusted p values for multiple tests
finalS$anova$padj <- p.adjust (finalS$anova$`Pr(>F)`, method = 'holm', n = 3)
finalS$anova

### checking if chlorophyll was different among shell types
View(CCAS.data)
chlShell<-lm(ChlAdensity ~ Type, data=CCAS.data)
summary(chlShell)
TukeyHSD(chlShell)
ggplot(CCAS.data) +geom_boxplot(aes(x=Type, y=ChlAdensity))

##### CCA Plots ========================================
#install.packages("devtools")
#devtools::install_github("gavinsimpson/ggvegan")
library(ggvegan); library(ggrepel); library(cowplot)

cow<-theme_cowplot()
theme_set(cow)
# Field Plot ----------------------------------
tidy.ccaF<-fortify(cca.F) #make me a single long dataframe
tidy.ccaFSp<-tidy.ccaF %>% filter(Score=="species") %>% #species dataframe
  left_join(FTaxaTable, by=c('Label'='Taxa'))
#important species dataframe
tidy.ccaFspespe<-tidy.ccaFSp %>% filter(Label %in% unique(c(t(CCAF.impsp)))) 
#sites dataframe
tidy.ccaFSi<-tidy.ccaF %>% filter(Score=="sites") %>% cbind(field.env[,c(2,4)])
#environmental dataframe
tidy.ccaFcon<-tidy.ccaF %>% filter(Score!="species",
                                   Score!="sites",
                                   Score!="constraints")%>%
  mutate(Tlabel=factor(c("mussel biomass","chl. a","discharge","substrate"),
                       levels=c("mussel biomass","chl. a","discharge","substrate"),
                       labels=c(expression(paste("mussel biomass")),
                                expression('chlorophyll '*italic(a)),
                                "discharge",
                                expression(D[var])))) #make pretty labels
frep<-rbind(tidy.ccaFspespe[c(1:8,12)],tidy.ccaFcon) #dataframe for labels

FieldCCA<-ggplot()+
  #plot all species
  geom_point(data=tidy.ccaFSp, aes(x=CCA1, y=CCA2), color="grey",
             alpha=.5)+
  #plot important species as black points
  geom_point(data=tidy.ccaFspespe, aes(x=CCA1,y=CCA2)) +
  #plot the environ. variables
  geom_segment(data = tidy.ccaFcon[tidy.ccaFcon$Score=="biplot",],
               aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.0, "cm")), alpha=.7)+
  scale_x_continuous(name="CCA 1 : 35.7%")+
  scale_y_continuous(name="CCA 2 : 10.6%") +
  scale_color_manual(values=c("red","black"), guide=F) +
  scale_size_manual(values=c(3.5,2.8), guide=F) #changing label sizes
FieldCCA

# Enclosure Plot -------------------------------
tidy.ccaE<-fortify(cca.E) #long dataframe of CCA results
#dataframe of all species
tidy.ccaESp<-tidy.ccaE %>% filter(Score=="species") %>%
  left_join(FTaxaTable, by=c('Label'='Taxa'))
#dataframe of important species
tidy.ccaEspespe<-tidy.ccaE %>% 
  filter(Score=="species",
         Label %in% unique(c(t(CCAE.impsp)))) %>%
  left_join(FTaxaTable, by=c('Label'='Taxa'))
#dataframe of environmental contraints
tidy.ccaEcon<-tidy.ccaE %>% filter(Score!="species",
                                   Score!="sites",
                                   Score!="constraints")%>%
  mutate(Tlabel=factor(c("ActBM","AmbBM","Alive","Discharge","ChlA",
                         "substrate","Not Alive","Alive"),
                       levels=c("ActBM","AmbBM","Alive","Discharge","ChlA",
                                "substrate","Not Alive","Alive"),
                       labels=c(expression(paste(italic("A.ligamentina")," biomass")),
                                expression(paste(italic("A.plicata")," biomass")),
                                "alive",
                                "discharge",
                                expression('chlorophyll '*italic(a)),
                                expression(D[var]),
                                expression(paste("not alive")),
                                "alive")))
erep<-rbind(tidy.ccaEspespe[c(1:8,12)],tidy.ccaEcon) %>% #label dataframe
  filter(Score!="biplot" | Tlabel!="alive")

EncCCA<-ggplot()+
  #plot all species grey
  geom_point(data=tidy.ccaESp, aes(x=CCA1, y=CCA2), color="grey",
             alpha=.5)+
  #plot imp. species black
  geom_point(data=tidy.ccaEspespe, aes(x=CCA1,y=CCA2)) +
  #add environmental vectors
  geom_segment(data = tidy.ccaEcon[tidy.ccaEcon$Score=="biplot",],
               aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.0, "cm")))+
  #label environmental and imp. species
  scale_x_continuous(name="CCA 1 : 5.3%")+
  scale_y_continuous(name="CCA 2 : 3.0%")+
  scale_color_manual(values=c("red","navy","black")) +
  scale_size_manual(values=c(4,4,2.8))+
  theme(legend.position = "none")
EncCCA

#Shell Plot ---------------------------------------
tidy.ccaS<-fortify(cca.S) #long data fram of CCA results
#data frame of species loadings
tidy.ccaSSp<-tidy.ccaS %>% filter(Score=="species") %>%
  left_join(FTaxaTable, by=c('Label'='Taxa'))
#dataframe of imp. species loadings
tidy.ccaSspespe<-tidy.ccaS %>% filter(Score=="species",
                                  Label %in% c("Col.ElmA","Ostracoda",
                                               "Eph.LeptA","Tri.Helicophyche",
                                               "Tri.Hydroptilidae","Gas.Planorbidae",
                                               "Dip.ChironomidaeL","Tri.Lepidostomatidae",
                                               "Dip.Tipulidae","Hydrozoans.miscH")) %>%
  left_join(FTaxaTable, by=c('Label'='Taxa'))
#dataframe of environmental variables
tidy.ccaScon<-tidy.ccaS %>% filter(Score!="species",
                                   Score!="sites",
                                   Score!="constraints") %>%
  mutate(Tlabel=factor(c("Amblema","Sham","ChlA",
                         "Actinonaias","Amblema","Live","Sham"),
                levels=c("Amblema","Sham","ChlA",
                         "Actinonaias","Amblema","Live","Sham"),
                labels=c(expression(italic("A. plicata")),
                         "sham",
                         expression('chlorophyll '*italic(a)),
                         expression(italic("A. ligamentina")),
                         expression(italic("A. plicata")),
                         "alive","not alive")))
srep<-rbind(tidy.ccaSspespe[,c(1:8,12)], tidy.ccaScon) %>% #dataframe for labels 
  filter(Score!="biplot" | Label=="ChlAdensity")

shellCCA<-ggplot()+
  geom_point(data=tidy.ccaSSp, aes(x=CCA1, y=CCA2), color="grey",
             alpha=.5)+
  geom_point(data=tidy.ccaSspespe, aes(x=CCA1,y=CCA2)) +
  geom_segment(data = tidy.ccaScon[tidy.ccaScon$Score=="biplot",],
               aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.0, "cm")))+
  scale_x_continuous(name="CCA 1 : 6.9%")+
  scale_y_continuous(name="CCA 2 : 4.8%")+
  scale_color_manual(values=c("red","navy","black"), guide=F) +
  scale_size_manual(values=c(4,4,2.8), guide=F)
shellCCA
#looking at all species for the discussion
#comparing to other manuscripts, which taxa did we find
ggplot()+
  geom_point(data=tidy.ccaSSp, aes(x=CCA1, y=CCA2), color="grey")+
  geom_point(data=tidy.ccaSspespe, aes(x=CCA1,y=CCA2)) +
  geom_segment(data = tidy.ccaScon[tidy.ccaScon$Score=="biplot",],
               aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.0, "cm")))+
  geom_text_repel(data=srep[srep$Score!="species",],
                  aes(x=CCA1,y=CCA2, label=Tlabel,
                             color=Score))+
  geom_text_repel(data=tidy.ccaSSp, aes(x=CCA1,y=CCA2, label=Tlabel))+
  scale_color_manual(values=c("red","blue"), guide=F)


# CCA with all text plot ---------------------------------
FCCA_atext<-FieldCCA+
  geom_text_repel(data=frep, aes(x=CCA1,y=CCA2, #label everything
                               label=Tlabel, color=Score, size=Score),
                parse=T, box.padding =.3)

ECCA_atest<-EncCCA + geom_text_repel(data=erep, aes(x=CCA1,y=CCA2, label=Tlabel,
                               color=Score, size=Score),
                parse=T, direction="both")

SCCA_atext<-shellCCA + 
  geom_text_repel(data=srep, aes(x=CCA1,y=CCA2, label=Tlabel,
                                 size=Score, color=Score),
                          parse=T)

cca_atext<-plot_grid(FCCA_atext, ECCA_atest, SCCA_atext, nrow=1, labels="AUTO") #plot the boxes together
legend1<-get_legend(ECCA_atest+ #pull the legend out and make it better
                      theme(legend.position = c(.3,.8),
                            legend.direction = "horizontal"))
plot_grid(cca_atext,legend1, nrow=2, rel_heights = c(1,.1)) #print both boxes and label
ggsave("./Figures/CCAplotsAllText.tiff", width=11, height=4)

# CCA with just important sp text plots ----------
FCCA_sp<-FieldCCA+
  geom_text_repel(data=frep[frep$Score=="species",], 
                  aes(x=CCA1,y=CCA2, #label the imp. species
                      label=Tlabel),size=3, parse=T)

ECCA_sp<-EncCCA + 
  geom_text_repel(data=erep[erep$Score=="species",], 
                  aes(x=CCA1,y=CCA2, label=Tlabel),size=3, parse=T)

SCCA_sp<-shellCCA + 
  geom_text_repel(data=srep[srep$Score=="species",], 
                  aes(x=CCA1,y=CCA2, label=Tlabel), size=3,parse=T)

plot_grid(FCCA_sp, ECCA_sp, SCCA_sp, nrow=1, 
                  labels="AUTO") #plot the boxes together
ggsave("./Figures/CCAplotsspText2.tiff", width=11, height=4)

