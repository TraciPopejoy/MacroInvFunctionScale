

#long data frame of Enc, Sp, Taxa, Count
EnclComMat<-EnclInv %>% group_by(TEid,Enc, Week, Taxa) %>% 
  summarize(N=n()) %>%
  left_join(ECounts) %>% mutate(IndDensity.m2=N/AreaSamp)%>%
  select(Enc,TreatA, Week,TEid, Taxa, IndDensity.m2) %>% 
  mutate(id=1:n()) %>% filter(Taxa!="MYSTERY") %>%
  spread(Taxa,IndDensity.m2) %>% group_by(TEid, Enc, TreatA) %>%
  summarise_if(is.numeric,funs(mean(.[which(!is.na(.))]))) %>%
  select(-id) %>% replace(is.na(.),0)

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




#### exploring family abundance ####
head(EnclComMat)
EnclPCom %>%
 gather(Taxa, Abundance, -TEid) %>% group_by(Taxa) %>%
 summarize(mean.A=mean(Abundance, na.rm=T)*100)

EnclComMat %>%
 gather(Taxa, Abundance, -TEid) %>% group_by(Taxa) %>%
 summarize(mean.A=mean(Abundance, na.rm=T), sum.A=sum(Abundance, na.rm=T))

EnclGraph<-EnclComMat %>% 
  mutate(Enc2=substr(TEid, 4,6), Week=substr(TEid, 1,3)) %>%
  left_join(TreatENC) %>% left_join(MusBiomass) %>%
  mutate(MusselBiomass.g.m2=sumBM/.25) %>% replace(is.na(.),0)

#Heptageniidae
ggplot(EnclGraph, aes(x=MusselBiomass.g.m2,y=Eph.Heptageniidae))+
  geom_smooth(method="lm",color="black")+
  geom_point(size=2, aes(color=TreatA))
summary(lm(Eph.Heptageniidae~MusselBiomass.g.m2, data=EnclGraph))
#ChironomidaeL - FFG minx
ggplot(EnclGraph, aes(x=MusselBiomass.g.m2,y=Dip.ChironomidaeL))+
  geom_smooth(method="lm",color="black")+
  geom_point(size=2, aes(color=TreatA))
summary(lm(Dip.ChironomidaeL~MusselBiomass.g.m2, data=EnclGraph))
#Elmidae larvae
ggplot(EnclGraph, aes(x=MusselBiomass.g.m2,y=Col.ElmL))+
  geom_smooth(method="lm",color="black")+
  geom_point(size=2, aes(color=TreatA))
summary(lm(Col.ElmL~MusselBiomass.g.m2, data=EnclGraph))
#Polycentropidae
ggplot(EnclGraph, aes(x=MusselBiomass.g.m2,y=Tri.Polycentropidae))+
  geom_smooth(method="lm",color="black")+
  geom_point(size=2, aes(color=TreatA))
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

#Enclosure Sum

EnclSum<-EnclInv %>% group_by(Taxa) %>% 
  summarize(N=n()) %>%
  arrange(desc(N)) %>% mutate(rank=1:n())
ggplot(EnclSum, aes(x=rank, y=N))+geom_bar(stat="identity")+
  xlab("Abundance Rank")+ylab("Total Count")
