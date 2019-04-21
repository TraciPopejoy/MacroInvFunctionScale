
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
  spread(Taxa,IndDensity.cm2) %>% group_by(SamID)%>%
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


#### exploring family abundance ####
head(ShellComMat)
View(ShellPCom %>%
       gather(Taxa, Abundance, -SamID) %>% group_by(Taxa) %>%
       summarize(mean.A=mean(Abundance, na.rm=T)*100))

View(ShellComMat %>%
       gather(Taxa, Abundance, -SamID) %>% group_by(Taxa) %>%
       summarize(mean.A=mean(Abundance, na.rm=T), sum.A=sum(Abundance, na.rm=T)))

##### biomass #####
MusBioShell<-MusselData %>% left_join(TreatENC) %>%
  mutate(BiomassE=case_when(Genus=="AMB"~ 1.57e-6*L^3.21,
                            Genus=='ACT'~ 4.7e-7*L^3.51)) %>%
  group_by(Enc2, Genus, Treatment) %>%
  summarize(sumBM=sum(BiomassE, na.rm=T), 
            meanBM=mean(BiomassE,na.rm=T), 
            sdBM=sd(BiomassE, na.rm=T))

ShellGraph<-ShellComMat %>% 
  mutate(Enc2=substr(SamID, 1,3), Genus=substr(SamID, 5,7)) %>%
  left_join(TreatENC) %>% left_join(MusBioShell) %>%
  mutate(MusselBiomass.g.m2=sumBM/.25,
         Treat.Shell=paste(Type,Genus, sep="."))

#ChironomidaeL - FFG minx
ggplot(ShellGraph, aes(x=MusselBiomass.g.m2,y=Dip.ChironomidaeL))+
  geom_smooth(method="lm", color="black")+
  geom_point(size=2, aes(color=Treat.Shell))
#Polycentropidae - Predator
ggplot(ShellGraph, aes(x=MusselBiomass.g.m2,y=Tri.Polycentropidae))+
  geom_smooth(method="lm", color="black")+
  geom_point(size=2, aes(color=Treat.Shell))
#Heptageniidae
ggplot(ShellGraph, aes(x=MusselBiomass.g.m2,y=Eph.Heptageniidae))+
  geom_smooth(method="lm", color="black")+
  geom_point(size=2, aes(color=Treat.Shell))
#Elmidae larvae
ggplot(ShellGraph, aes(x=MusselBiomass.g.m2,y=Col.ElmL))+
  geom_smooth(method="lm", color="black")+
  geom_point(size=2, aes(color=Treat.Shell))

ggplot(ShellGraph, aes(x=log1p(Dip.ChironomidaeL)))+geom_histogram()
ggplot(ShellGraph, aes(x=log1p(Tri.Polycentropidae)))+geom_histogram()

log1p.ShellCom<-ShellComMat %>% #taxa >0.1%
  select(-SamID) %>% replace(is.na(.),0) %>%
  mutate_if(is.numeric, .funs=log1p)
shell.pca<-prcomp(log1p.ShellCom)
plot(shell.pca, type="l")
summary(shell.pca)
print(shell.pca)

ShellGraph<-ShellGraph %>% mutate(PCA1=shell.pca$x[,1],
                                  PCA2=shell.pca$x[,2])
library(ggsci)
ggplot(ShellGraph, aes(x=PCA1, y=PCA2))+
  geom_point(aes(color=Treat.Shell))+
  scale_color_futurama()
ggplot(ShellGraph, aes(x=MusselBiomass.g.m2, y=PCA2))+
  geom_point(aes(color=Type))+geom_smooth(method="lm")
summary(lm(PCA1~MusselBiomass.g.m2*Type, data=ShellGraph))
PCAline<-data.frame(species=rownames(shell.pca$rotation),Scores=shell.pca$rotation[,2])
ggplot(PCAline)+
  geom_text(aes(-1,Scores,label=round(Scores,3)))+
  geom_text(aes(1,Scores, label=species))+
  geom_vline(xintercept=0)+
  scale_y_continuous(trans="log")+
  coord_cartesian(xlim=c(-3,3))+
  labs(x=NULL,y=NULL)+theme_classic()

#### nmds ####
head(SCounts)
ShellCom<-SCounts %>% ungroup() %>% group_by(SamID) %>% spread(Taxa, RA) %>%
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


###Shell Sum
ShellSum<-SInv %>% group_by(Enc2, ShellSpecies, Taxa) %>% 
  mutate(SamID=paste(Enc2,ShellSpecies, sep=".")) %>% summarize(N=n()) %>%
  group_by(Taxa) %>% summarize_if(is.numeric,sum, na.rm=T)%>%
  arrange(desc(N)) %>% mutate(rank=1:n())
ggplot(ShellSum, aes(x=rank, y=N))+geom_bar(stat="identity")+
  xlab("Abundance Rank")+ylab("Total Count")
  

                                       