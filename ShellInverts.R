#### Invertebrates on Shells ####

library(readxl);library(tidyverse)
#bring in treatment data from enclosures
TreatENC<-read_excel("../FEn17/FEn17_data/FEn17OKTreatments.xlsx") 
#bring in mussel data from enclosures
MusselData<-read_csv("../FEn17/FEn17_data/MusselBMExpFEn17OK.csv")
names(MusselData)[1]<-"Enclosure"
#calculate shell surface area by enclosure & species
MDat<-MusselData %>% full_join(TreatENC)%>% 
  select(Enc2,Genus, TreatA, Type,Spp,L, H,W) %>% 
  mutate(ShellSpecies=recode(Genus, "AMBL"="AMB", "ACT"="ACT"),
         ShellSurArea.mm2=2*(L*H)+2*(L*W)+2*(H*W),
         SamID=paste(Enc2, ShellSpecies, sep=".")) %>%
  group_by(Enc2,TreatA,SamID) %>% 
  summarize(TShellSurArea.cm2=sum(ShellSurArea.mm2, na.rm=T)*.01)
#the ones with 0 I didn't measure height on them b/c not needed

#bring in macroinvertebrates lengths on shell data
shellInv<-read_excel("../FEn17/FEn17_data/ShellInv.xlsx")
unique(shellInv$Taxa)[order(unique(shellInv$Taxa))] #check taxa for correct spelling
#long list of macroinvertebrate lengths
SInv<-shellInv %>% filter(Taxa!="Spider") %>% 
  select(Enc2, ShellSpecies, SamType, Taxa, Length.cm)

#macroinvertebrate summary on shells
SCounts<-SInv %>% group_by(Enc2, ShellSpecies) %>% 
  summarize(Nsam=n(),
            richness=length(unique(Taxa)))%>%
  mutate(SamID=paste(Enc2,ShellSpecies, sep="."))%>% full_join(MDat) %>%
  mutate(InvertDensity.npcm2=Nsam/TShellSurArea.cm2) %>% full_join(TreatENC)
#long data frame of Enc, Sp, Taxa, Count
ShellComMat<-SInv %>% group_by(Enc2, ShellSpecies, Taxa) %>% 
  mutate(SamID=paste(Enc2,ShellSpecies, sep=".")) %>% summarize(N=n()) %>%
  left_join(SCounts) %>% mutate(IndDensity.cm2=N/TShellSurArea.cm2,
                                SamID=paste(Enc2, ShellSpecies, sep=".")) %>%
  select(Enc2, ShellSpecies,SamID, Taxa, IndDensity.cm2) %>% 
  mutate(id=1:n()) %>% 
  spread(Taxa,IndDensity.cm2) %>% group_by(SamID) %>%
  summarise_if(is.numeric,funs(mean(.[which(!is.na(.))]))) %>%
  select(-id)

##need proportional abundance
ShellPCom<-SInv %>% group_by(Enc2, ShellSpecies, Taxa) %>% 
  filter(Taxa!="MYSTERY", Taxa!="MysteryTricoptera", Taxa!="Col.Myst",
         Taxa!="Tri.Pupa", Taxa!="Dip.Adult", Taxa!="Orthoptera") %>% 
  summarize(N=n()) %>%  mutate(id=1:n(),
                               SamID=paste(Enc2,ShellSpecies, sep=".")) %>% 
  spread(Taxa,N) %>% group_by(SamID)%>%
  summarise_if(is.numeric,funs(mean(.[which(!is.na(.))]))) %>%
  select(-id) %>% replace(is.na(.),0) %>%
  mutate(rowsum=rowSums(.[,-1])) %>% ungroup()%>% group_by(SamID)%>%
  summarise_if(is.numeric,funs(./rowsum)) %>%
  select(-rowsum) %>% 
  select(SamID, everything())

FTaxaTable<-read_excel("Macroinv Power Law Coeffs TBP.XLSX", sheet=2)
match(names(ShellComMat[,-1]), FTaxaTable$Taxa) #checking to make sure all taxa in taxa table


#so we have a species community matrix ShellComMat
#need a trait table of those species
ShellTraits<-FTaxaTable[match(names(ShellComMat[,-1]), FTaxaTable$Taxa),]
ShellTraits$Taxa==names(ShellComMat[,-1]) #need it to all be true
ShellTraits <- ShellTraits %>%
  filter(T.TropP!=0) %>%
  select(Taxa,T.Habit, T.TropP, T.TropS,
         T.dev,T.Lifespan,T.crwl,T.swim,T.MatSize)


#### nmds ####
head(SCountsT)
ShellCom<-SCountsT %>% ungroup() %>% group_by(SamID) %>% spread(Taxa, RA) %>%
  select(-richness, -TShellSurArea.cm2, -InvertDensity.npcm2,-Enclosure,
         -TreatA, -TreatF, -Type, -Spp,-Enc2, -ShellSpecies, -Nsam, -N) %>%
  summarise_all(funs(mean(., na.rm = TRUE)))
ShellCom<-ShellCom[rowSums(ShellCom[,-1], na.rm=T)==100,] %>% 
  mutate(Enc2=substr(SamID, 1,3),
         ShellSpp=substr(SamID, 5,8))
ShellCom[,-c(1,35:37)] <- replace(ShellCom[,-c(1,35:37)], is.na(ShellCom[,-c(1,35:37)]), 0)
library(vegan)
nmdsShell<-metaMDS(ShellCom[,-c(1,35:37)])
plot(nmdsShell)
data.scoresShell <- as.data.frame(scores(nmdsShell))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scoresShell$site<-as.vector(ShellCom$SamID)
data.scoresShell$ShellSpp <- as.vector(ShellCom$ShellSpp)  #  add the grp variable created earlier
data.scoresShell$Enc2 <- as.vector(ShellCom$Enc2)  #  add the grp variable created earlier
data.scoresShell<- data.scoresShell %>% left_join(Treat)
head(data.scoresShell)  #look at the data

species.scores <- as.data.frame(scores(nmdsShell, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data

grp.amb <- data.scoresShell[data.scoresShell$ShellSpp == "AMB", ][chull(data.scoresShell[data.scoresShell$ShellSpp == 
                                                                      "AMB", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.act <- data.scoresShell[data.scoresShell$ShellSpp == "ACT", ][chull(data.scoresShell[data.scoresShell$ShellSpp == 
                                                                                           "ACT", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
spp.hull.data <- rbind(grp.amb, grp.act)  #combine grp.a and grp.b

grp.ambl <- data.scoresShell[data.scoresShell$TreatA == "AMBL", ][chull(data.scoresShell[data.scoresShell$TreatA == 
                                                                                           "AMBL", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.actl <- data.scoresShell[data.scoresShell$TreatA == "ACTL", ][chull(data.scoresShell[data.scoresShell$TreatA == 
                                                                                           "ACTL", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
grp.ambs <- data.scoresShell[data.scoresShell$TreatA == "AMBS", ][chull(data.scoresShell[data.scoresShell$TreatA == 
                                                                                            "AMBS", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.acts <- data.scoresShell[data.scoresShell$TreatA == "ACTS", ][chull(data.scoresShell[data.scoresShell$TreatA == 
                                                                                            "ACTS", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
treat.hull.data <- rbind(grp.ambl, grp.actl, grp.ambs, grp.acts)  #combine grp.a and grp.b


ggplot() + 
  geom_polygon(data=treat.hull.data,aes(x=NMDS1,y=NMDS2,fill=TreatA,group=TreatA),alpha=0.40) + # add the convex hulls
  geom_polygon(data=treat.hull.data,aes(x=NMDS1,y=NMDS2,color=ShellSpp,group=ShellSpp),alpha=0.20) + # add the convex hulls
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.7, size=2) +  # add the species labels
  geom_point(data=data.scoresShell,aes(x=NMDS1,y=NMDS2,shape=ShellSpp),size=4) + # add the point markers
  scale_fill_manual(values=CP[c(1,2,5,4)])+
  coord_equal() +
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

##### biomass #####
MusBioShell<-MusselData %>% left_join(TreatENC) %>%
  mutate(BiomassE=case_when(Genus=="AMB"~ 1.57e-6*L^3.21,
                            Genus=='ACT'~ 4.7e-7*L^3.51)) %>%
  group_by(Enc2, Genus, Treatment) %>%
  summarize(sumBM=sum(BiomassE, na.rm=T), 
            meanBM=mean(BiomassE,na.rm=T), 
            sdBM=sd(BiomassE, na.rm=T))
