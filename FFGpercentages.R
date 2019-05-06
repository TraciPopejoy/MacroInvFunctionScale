head(FieldPCom)
head(EnclPCom)
head(ShellPCom)

FieldPerc<-F.ComPer %>% gather(Taxa, Per, -FSamID,-Reach,-SamplingSeason,
                               -Treatment,-Season,-Mussel.g.m2) %>% 
  left_join(FTaxaTable) %>%
  select(FSamID, Taxa, Per, Order, T.TropP) %>%
  group_by(FSamID, T.TropP) %>% summarise(TotalPer=sum(Per)) %>%
  left_join(fieldSiteKey) %>% ungroup() %>% 
  select(-Sample) %>%filter(!duplicated(.)) %>%
  mutate(FFG=case_when(T.TropP==1~"C-Gatherer",
                       T.TropP==2~"C-Filterer",
                       T.TropP==3~"Herbivore",
                       T.TropP==4~"Predator",
                       T.TropP==5~"Shredder",
                       T.TropP==6~"Parasite")) %>%
  left_join(FieldBMest, by="Reach")
FieldPerc %>% group_by(FSamID) %>% summarize(sum(TotalPer))
library(ggsci)

ggplot(FieldPerc,aes(x=Treatment.x,y=TotalPer, color=FFG, fill=FFG, shape=FFG))+
  stat_summary(position=position_dodge(width=.9), size=1.1)+
  scale_y_continuous(name="FFG Relative Abundance",
                     breaks=c(0,.20,.40,.60,.80,1.00),
                     labels=c("0%","20%","40%","60%","80%","100%"))+
  scale_color_manual(values = c("darkred","#c5693c","#a4a934","#6796ca",
                                "#8867d0","#c060a6"))+
    scale_shape_manual(values=c(15:18,15:17))+
  facet_wrap(~SamplingSeason)+theme_krementz()

ggplot(FieldPerc, aes(x=Treatment.x, y=TotalPer, fill=FFG))+
  geom_bar(stat="identity")+facet_wrap(~SamplingSeason)
ggplot(FieldPerc[FieldPerc$FFG=="C-Filterer",], 
       aes(x=Mussel.g.m2, y=TotalPer, color=Treatment.x, shape=SamplingSeason))+
 geom_point(size=3)+scale_color_startrek()
  #stat_summary(fun.y = mean, geom = "line",position=position_dodge(width=1.75))+
  #stat_summary(position=position_dodge(width=.5))

FieldPerc<-FieldPCom %>% gather(Taxa, Per, -FSamID) %>% left_join(FTaxaTable) %>%
  select(FSamID, Taxa, Per, Order, T.TropP) %>%
  group_by(FSamID, T.TropP) %>% summarise(TotalPer=sum(Per)) %>%
  left_join(fieldSiteKey) %>% ungroup() %>% 
  select(-Sample) %>%filter(!duplicated(.)) %>%
  mutate(FFG=case_when(T.TropP==1~"C-Gatherer",
                       T.TropP==2~"C-Filterer",
                       T.TropP==3~"Herbivore",
                       T.TropP==4~"Predator",
                       T.TropP==5~"Shredder",
                       T.TropP==6~"Parasite")) %>%
  left_join(FieldBMest, by="Reach")
FieldPerc
library(ggsci)
ggplot(FieldPerc, aes(x=Treatment, y=TotalPer, fill=FFG))+
  geom_bar(stat="identity")+facet_wrap(~SamplingSeason)
ggplot(FieldPerc[FieldPerc$FFG=="C-Filterer",], 
       aes(x=Mussel.g.m2, y=TotalPer, color=Treatment.x, shape=SamplingSeason))+
  geom_point(size=3)+scale_color_startrek()
#stat_summary(fun.y = mean, geom = "line",position=position_dodge(width=1.75))+
#stat_summary(position=position_dodge(width=.5))

###### Shell
ShellPerc<-S.ComPer %>% gather(Taxa, Per, -SamID,-Enc2,
                               -ShellSpecies,-TreatA,-Type,-Spp) %>% 
  left_join(FTaxaTable) %>%
  group_by(SamID, T.TropP) %>% summarise(TotalPer=sum(Per)) %>%
  mutate(Enc2=substr(SamID,1,3), Genus=substr(SamID, 5,7))%>%
  left_join(TreatENC) %>%
  mutate(FFG=case_when(T.TropP==1~"C-Gatherer",
                       T.TropP==2~"C-Filterer",
                       T.TropP==3~"Herbivore",
                       T.TropP==4~"Predator",
                       T.TropP==5~"Shredder",
                       T.TropP==6~"Parasite"))
ShellPerc
library(ggsci)
ggplot(ShellPerc,aes(x=Genus,y=TotalPer, color=FFG, fill=FFG, shape=FFG))+
  stat_summary(position=position_dodge(width=.9), size=1.1)+
  scale_y_continuous(name="FFG Relative Abundance",
                     breaks=c(0,.20,.40,.60,.80,1.00),
                     labels=c("0%","20%","40%","60%","80%","100%"))+
  scale_color_manual(values = c("darkred","#c5693c","#a4a934","#6796ca",
                                "#8867d0","#c060a6"))+
  scale_shape_manual(values=c(15:18,15:17))+
  facet_wrap(~Type)+theme_krementz()


ggplot(ShellPerc, aes(x=TreatA, y=TotalPer, fill=FFG))+
  geom_bar(stat="identity")+facet_wrap(~Genus+Type, scales="free_x")
ggplot(ShellPerc[ShellPerc$FFG=="C-Filterer",], 
       aes(x=sumBM, y=TotalPer, color=Genus, shape=Type))+
  geom_point(size=3)+scale_color_startrek()
#stat_summary(fun.y = mean, geom = "line",position=position_dodge(width=1.75))+
#stat_summary(position=position_dodge(width=.5))
