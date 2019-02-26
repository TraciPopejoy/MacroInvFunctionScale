head(FieldPCom)
head(EnclPCom)
head(ShellPCom)

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
ShellPerc<-ShellPCom %>% gather(Taxa, Per, -SamID) %>% left_join(FTaxaTable) %>%
  select(SamID, Taxa, Per, Order, T.TropP) %>%
  group_by(SamID, T.TropP) %>% summarise(TotalPer=sum(Per)) %>%
  mutate(Enc2=substr(SamID,1,3), Genus=substr(SamID, 5,7))%>%
  left_join(TreatENC) %>%
  mutate(FFG=case_when(T.TropP==1~"C-Gatherer",
                       T.TropP==2~"C-Filterer",
                       T.TropP==3~"Herbivore",
                       T.TropP==4~"Predator",
                       T.TropP==5~"Shredder",
                       T.TropP==6~"Parasite")) %>%
  left_join(MusBioShell)
ShellPerc
library(ggsci)
ggplot(ShellPerc, aes(x=Treatment, y=TotalPer, fill=FFG))+
  geom_bar(stat="identity")+facet_wrap(~Genus+Type, scales="free_x")
ggplot(ShellPerc[ShellPerc$FFG=="C-Filterer",], 
       aes(x=sumBM, y=TotalPer, color=Genus, shape=Type))+
  geom_point(size=3)+scale_color_startrek()
#stat_summary(fun.y = mean, geom = "line",position=position_dodge(width=1.75))+
#stat_summary(position=position_dodge(width=.5))
