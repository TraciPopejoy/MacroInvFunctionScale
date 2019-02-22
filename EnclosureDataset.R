#### Invertebrates on at Enclosure Level ####

library(readxl);library(tidyverse)
#bring in treatment data from enclosures
TreatENC<-read_excel("../FEn17/FEn17_data/FEn17OKTreatments.xlsx") 
#bring in the invert data from basket sampling
FEn17Inv12<-read.csv("../FEn17/FEn17_data/FEn17InvMeas.csv", stringsAsFactors = F) %>%
  filter(Label!="rulerxocc.tif") %>% select(-Mean, -Max, -Min) %>%
  mutate(TEid=substring(Label, 6,11),
         Enc=substring(Label, 9,11),
         Week=substring(Label, 6,8))
FEn17Inv09<-read.csv("../FEn17/FEn17_data/FEn17w09.csv", stringsAsFactors = F) %>% 
  filter(Label!="rulerxocc.tif") %>% select(-Mean, -Max, -Min, - SampleSplit)%>%
  mutate(Week="w09",
         TEid=paste(Week, Enc, sep=""))
EnclInv<-rbind(FEn17Inv12, FEn17Inv09) #contains every insect identified from baskets

unique(EnclInv$Taxa)[order(unique(EnclInv$Taxa))] #check taxa for correct spelling

#### Calculate area sampled ####
BaskSampled<-read.csv("../FEn17/FEn17_data/SubSamFEn17OK.csv", stringsAsFactors = F) %>%
  mutate(WeekG=case_when(Week==8~"w09",
                         Week==4~"w04",
                         Week==12 ~"w12")) %>% 
  left_join(TreatENC, by="Enclosure") %>% 
  filter(Enclosure!="WC" & Type.x=="Insect") %>%
  select(Enclosure, Enc2, Basket., WeekG) %>%
  mutate(TEid=paste(WeekG, Enc2, sep="")) %>% filter(!duplicated(TEid))
  
#Some filters didn't have baskets associated with them
#assumed 3 baskets
BaskSampled[is.na(BaskSampled$Basket.),"Basket."]<-3

#macroinvertebrate summary in enclosures
ECounts<-EnclInv %>% group_by(TEid, Enc, Week) %>% 
  summarize(Nsam=n(),
            richness=length(unique(Taxa))) %>% left_join(BaskSampled) %>%
  mutate(AreaSamp=(Basket.*.03315),
         InvertDensity.npm2=Nsam/AreaSamp) %>% 
  full_join(TreatENC) %>% 
  select(-Enclosure, -Enc2, -Basket., -WeekG)
#long data frame of Enc, Sp, Taxa, Count
EnclComMat<-EnclInv %>% group_by(TEid,Enc, Week, Taxa) %>% 
  summarize(N=n()) %>%
  left_join(ECounts) %>% mutate(IndDensity.m2=N/AreaSamp)%>%
  select(Enc, Week,TEid, Taxa, IndDensity.m2) %>% 
  mutate(id=1:n()) %>% filter(Taxa!="MYSTERY") %>%
  spread(Taxa,IndDensity.m2) %>% group_by(TEid) %>%
  summarise_if(is.numeric,funs(mean(.[which(!is.na(.))]))) %>%
  select(-id) 

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
  select(Taxa,T.Habit, T.TropP, T.TropS)
