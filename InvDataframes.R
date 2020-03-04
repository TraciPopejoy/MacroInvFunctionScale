# script to get invertebrate data into correct form 
# have four versions - density, percent abundance, biomass, percent biomass
# dataframes named: ScaleInitial.DataType.Quantifications
#  so F.Com.Dens is the field data community matrix in density form (n/m2)
# data frames  have associated treatment data and mussel biomass estimates

# TREATMENT TABLES: fieldSiteKey, TreatENC
# TAXA TABLES: FTaxaTable
# INV DATAFRAMES: 
# 		Field - F.ComDens,F.ComPer, F.BioDens, F.BioPer
#       Enclosure - E.ComDens, E.ComPer
#       Shell - S.ComDens, S.ComPer

library(tidyverse); library(ggplot2); library(readxl); library(vegan)

# Biomass Function ====================================================
BiomassReg<-read_excel("Macroinv Power Law Coeffs TBP.xlsx", sheet = 1) #sheet that holds regression coefficients
BMTaxaTable<-read_excel("Macroinv Power Law Coeffs TBP.xlsx", sheet = 4) #taxa table to match taxa to correct regression
#function to calculate biomass (in mg) given a length
biomass<-function(Fam, Ord, Length) {
  #tmass<-NULL
  #if not ID'd to family, use order level regressions
  if(is.na(match(Fam, BiomassReg$Family))){ 
    plcoeffs<-BiomassReg[BiomassReg$Order == Ord &
                           !is.na(BiomassReg$Order == Ord),]
  }else{ #pull all the regressions for that family
    plcoeffs<-BiomassReg[BiomassReg$Family==Fam&
                           BiomassReg$Order == Ord&
                           !is.na(BiomassReg$Family==Fam), ] }
  #which regressions were actually built for insects this size
  #idx2<- c(Length>=plcoeffs[,19] & Length<=plcoeffs[,20]) 
  #idx2[is.na(idx2)]<-T #if no size range listed, use the regression anyways
  d1<-plcoeffs#[idx2,] #row by row????
  mass<-d1$a*(Length^d1$b) #power law to determine biomass
  #print(paste(Fam, Ord, sep=" - "))
  mean(mass, na.rm=T) #returns the mean biomass estimate from all applied formulas
}

# Field Data ------------------------------
# count data for the surber samples; different tables for 2015, 2016
Count2016<-read_excel("./data/2016 Field Data.xlsx", sheet = "Counts") %>%
  mutate(Reach=paste(Site, Treatment, sep="-")) %>% 
  group_by(Sample, Location, Site, Treatment, Reach, SamplingSeason) %>% 
  summarize_if(is.numeric, sum, na.rm=T)
Count2015<-read_excel("./data/2015 Data Hotspot Field Surber.xlsx") %>% 
  mutate(Reach=paste(Site, Treatment, sep="-")) %>%
  group_by(Sample, Reach, Treatment, SamplingSeason, Site, Location) %>%
  summarize_if(is.numeric, sum, na.rm=T)
CountField<-rbind(Count2016, Count2015) #joining 2015 and 2016 so CountField has raw counts from each surber sample

#calculate total invertebrate density
CountField$InvDensity.npm2<-rowSums(CountField[,7:68], na.rm=T)/0.092903
CountField[is.na(CountField)]<-0 #replace na with 0 for vegan
CountField$richness<-specnumber(CountField[,7:68]) #calculate taxonomic richness
#good raw count dataframe for Field Data

# Sample ID key for building treatment table
fieldSiteKeyA<-CountField %>% group_by(Sample, Reach, SamplingSeason, Treatment) %>% 
  tally() %>% dplyr::select(-n) %>% 
  mutate(FSamID=paste(Reach,SamplingSeason, sep="."),
         Season=case_when(substr(SamplingSeason, 1,6)=="Summer"~"Summer",
                          substr(SamplingSeason, 1,6)=="Fall20"~"Fall"))
# field mussel biomass estimate
FieldBMest<-read_xlsx("./data/Recalculated data for TPD.xlsx", sheet=3) %>%
  select(-`OLD mass estimate g/ M2 with corbicula)`)
fieldSiteKey<-fieldSiteKeyA %>% left_join(FieldBMest, by=c("Reach","Treatment")) %>%
  mutate(SampSeaF=factor(SamplingSeason, levels=c("Summer2015","Fall2015","Summer2016","Fall2016")),
         Year=case_when(SamplingSeason=="Summer2015" |
                        SamplingSeason=="Fall2015"~2015,
                        SamplingSeason=="Summer2016" |
                        SamplingSeason=="Fall2016"~2016))
fieldSiteKey[fieldSiteKey$Reach=="K7-NM",c(8:10)]<-0
#fieldSiteKey is a good 'treatment' dataframe

#field data community matrix in density; units are n/m2
F.ComDens<-CountField %>% 
  group_by(Reach, SamplingSeason) %>% 
  summarize_if(is.numeric, mean) %>% 
  group_by(Reach, SamplingSeason) %>%
  summarise_if(is.numeric, list(~./0.092903)) %>%
  mutate(FSamID=paste(Reach, SamplingSeason, sep="."),
         InvDensity.npm2R=InvDensity.npm2*0.092903) %>% group_by(FSamID)%>%
  left_join(fieldSiteKey) %>% ungroup() %>% 
  dplyr::select(-Sample, -InvDensity.npm2) %>%
  filter(!duplicated(.))%>%
  dplyr::select(-`N quadrats`, -`StdDev of STDM (g.m-2)2`,
                -Site, -richness) %>%
  dplyr::select(FSamID,Reach,SamplingSeason, Treatment,Season,`Average of STDM (g.m-2)`,
                InvDensity.npm2R,everything()) %>%
  replace(is.na(.),0)
# field data community matrix in relative abundance; units are %
F.ComPer<-CountField %>% group_by(Reach, SamplingSeason) %>% 
  summarize_if(is.numeric, mean) %>%
  mutate(FSamID=paste(Reach, SamplingSeason, sep="."))%>% ungroup()%>%
  select(-InvDensity.npm2,-Other.other,-Tri.Pupa,-richness,-SamplingSeason,-Reach) %>%
  gather(Taxa,Count, -FSamID) %>% group_by(FSamID, Taxa) %>% 
  filter(Taxa!="MYSTERY", Taxa!="MysteryTricoptera", Taxa!="Col.Myst",
         Taxa!="Tri.Pupa", Taxa!="Dip.Adult", Taxa!="Dip.Cyclorrhaphous",
         Taxa!="Dip.Other", Taxa!="Dip.Orthorrhaphous", Taxa!="Dip.Thaumaleidae",
         Taxa!="Hymenoptera.miscH", Taxa!="Dip.Ascillidae") %>%
  mutate(id=1:n()) %>% 
  spread(Taxa,Count) %>%  dplyr::select(-id) %>% ungroup() %>%
  mutate(rowsum=rowSums(.[,-1])) %>% group_by(FSamID) %>%
  summarise_if(is.numeric,list(~./rowsum))%>%
  left_join(fieldSiteKey) %>% ungroup() %>% select(-Sample, -rowsum) %>%
  filter(!duplicated(.))%>%
  dplyr::select(-`N quadrats`, -`StdDev of STDM (g.m-2)2`, -Site) %>%
  dplyr::select(FSamID,Reach,SamplingSeason, Treatment,Season,`Average of STDM (g.m-2)`,
                everything())

#like biomass taxa table, but 'correct' taxonomic information
#biomass taxa table has 'misc' for some orders to filter (b/c no regressions for that order)
FTaxaTable<-read_excel("Macroinv Power Law Coeffs TBP.XLSX", sheet=2)
#checking to make sure all taxa in taxa table
which(is.na(match(names(F.ComPer[,-c(1:6,61,62)]), FTaxaTable$Taxa)))

# Field Biomass estimate ===============================================
# lengths measured from the surber data
# all organisms measured up to 100 or 10%, whichever is larger
Lengths2015<-read_excel("./data/2015 Data Hotspot Field Surber.xlsx", sheet= 3,
                        col_types = c("text","text","text","date","text",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "text")) %>%
  filter(substr(Sample,1,11) != "CALIBRATION")
#(Mtax15<-colnames(Lengths2015)[6:50])

#this takes about two minutes because 12300 observations
melt2015lengths<-Lengths2015 %>% 
  group_by(Sample, Reach, SamplingSeason) %>%
  dplyr::select(-`Surber #`, -`Date-collected`, -NOTES, -`average scale`, -Unknown) %>% 
  gather(Taxa, Length.mm, -Sample, -Reach, -SamplingSeason) %>% #long dataset
  drop_na() %>% #remove incomplete cases
  left_join(BMTaxaTable) %>% 
  filter(Order!="misc") %>% #no regressions for these
  rowwise() %>% #throw each row into the mutate
  mutate(BM.mg=mean(biomass(Family, Order, Length.mm))) #calculate biomass at each row

#because not all were measured, taking an average biomass to multiply by the count
samp2015count<-Count2015 %>% ungroup() %>% dplyr::select(-Location, -Site, -Treatment, 
                                                  -Reach, -Other.other) %>%
  gather(Taxa, Count, -Sample, -SamplingSeason) %>% filter(Count!=0)
# 2015 estimate of biomass based on mean biomass * count for each sample and taxa
SxT2015Long<-melt2015lengths %>% group_by(Sample,Taxa) %>% 
  summarize(maxBM=max(BM.mg), 
            minBM=min(BM.mg), 
            meanBM=mean(BM.mg, na.rm=T),
            SDBM=sd(BM.mg, na.rm=T)) %>% 
  left_join(samp2015count) %>%
  mutate(tbm.mg.avg=meanBM*Count)

# 2016 biomass invertebrate calculations
Lengths2016<-read_excel("./data/2016 Field Data.xlsx", sheet= 2,
                        col_types = c("text","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",#5
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",#10
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric",
                                      "numeric","numeric","numeric","numeric")) 

#this takes about two minutes because 17200 observations
melt2016lengths<-Lengths2016 %>% group_by(Sample) %>% 
  gather(Taxa, Length.cm, -Sample)%>% 
  drop_na() %>% #remove na from dataset
  left_join(BMTaxaTable) %>% 
  filter(Order!="misc")%>% #no regression
  rowwise()%>% #throw each row into the mutate
  mutate(Length.mm=Length.cm*10, #convert to mm
         BM.mg=mean(biomass(Family, Order, Length.mm)))  #calculate biomass for each length

samp2016count<-Count2016 %>% ungroup() %>% 
  dplyr::select(-Location, -Site, -Treatment, -Reach,-Other.other) %>%
  gather(Taxa, Count, -Sample, -SamplingSeason) %>% filter(Count!=0) %>%
  filter(!is.na(Count))

# 2016 estimate of biomass based on mean biomass * count for each sample and taxa
SxT2016Long<-melt2016lengths %>% group_by(Sample,Taxa) %>% 
  summarize(maxBM=max(BM.mg), 
            minBM=min(BM.mg), 
            meanBM=mean(BM.mg, na.rm=T),
            SDBM=sd(BM.mg, na.rm=T)) %>% 
  left_join(samp2016count) %>%
  mutate(tbm.mg.avg=meanBM*Count)
# SxT2016Long %>% filter(is.na(meanBM)) #likely samples that are to small/large

# field biomass density units is mg/m2
F.BioDens<-rbind(SxT2015Long, SxT2016Long) %>% 
  select(Sample, Taxa, tbm.mg.avg) %>%
  spread(Taxa, tbm.mg.avg) %>% left_join(fieldSiteKey[,1:5]) %>%
  group_by(FSamID,Reach, SamplingSeason) %>%
  replace(is.na(.),0) %>%
  summarize_if(is.numeric, mean)%>%
  summarize_if(is.numeric, list(~./0.092903)) %>%
  left_join(fieldSiteKey)%>% select(-Sample)%>% ungroup() %>%
  filter(!duplicated(.)) %>%
  select(-Site,-`StdDev of STDM (g.m-2)2`, -`N quadrats`)%>%
  select(FSamID,Reach,SamplingSeason, Treatment,Season,`Average of STDM (g.m-2)`, everything())

F.BioPer<-rbind(SxT2015Long, SxT2016Long) %>% 
  select(Sample, Taxa, tbm.mg.avg) %>%
  spread(Taxa, tbm.mg.avg) %>% left_join(fieldSiteKey[,1:5]) %>%
  group_by(FSamID,Reach, SamplingSeason) %>%
  replace(is.na(.),0) %>%
  summarize_if(is.numeric, mean)%>% ungroup() %>%
  mutate(rowsum=rowSums(.[,-c(1:3)])) %>% group_by(FSamID)%>%
  summarise_if(is.numeric,list(~./rowsum)) %>% 
  left_join(fieldSiteKey) %>% select(-rowsum, -Sample) %>%
  filter(!duplicated(.)) %>%
  select(-Site,-`StdDev of STDM (g.m-2)2`, -`N quadrats`)%>%
  select(FSamID,Reach,SamplingSeason, Treatment,Season,`Average of STDM (g.m-2)`, everything())

# Enclosure Treatment / Biomass data ----------------------
#treatment data frame
TreatENC<-read_excel("../FEn17/FEn17_data/FEn17OKTreatments.xlsx") 
# mussel regressions
# using the tutorial from Van ee et al. 2019
# saved the values to the following file

# bring in mussel data from enclosures
MusselData<-read_csv("../FEn17/FEn17_data/MusselBMExpFEn17OK.csv")
names(MusselData)[1]<-"Enclosure"
#calculate shell surface area by enclosure & species
### Nov 20 2019 - wrapped foil around shells to build correction #
SAcheck<-read_excel("./data/SurfaceAreaCheck.xlsx", sheet = "Sheet1") %>% #raw data
  mutate(ShellSurArea.mm2=2*(Length.mm*Height.mm)+2*(Length.mm*Width.mm)+2*(Height.mm*Width.mm)) #estimated surface area
 
#plot relationship between foil surface area and estimated surface area
ggplot(SAcheck, aes(x=SA, y=ShellSurArea.mm2, color=Genus))+
  geom_smooth(method="lm")+geom_point()+
  xlab("Foil Estimated Surface Area")+
  ylab("Quadrate Estimate Surface Area")
saMod<-lm(ShellSurArea.mm2~SA*Genus, data=SAcheck)
summary(saMod)
lm(ShellSurArea.mm2~SA:Genus, data=SAcheck) #species specific regression

# calculating surface area and applying species specific correction
MDat<-MusselData %>% full_join(TreatENC)%>% 
  dplyr::select(Enc2,Genus, TreatA, Type,Spp,L, H,W) %>% 
  mutate(ShellSpecies=recode(Genus, "AMBL"="AMB", "ACT"="ACT"),
         ShellSurArea.mm2=2*(L*H)+2*(L*W)+2*(H*W),
         corSA.mm2=case_when(Genus=="AMB"~ShellSurArea.mm2/1.750,
                             Genus=="ACT"~ShellSurArea.mm2/1.595),
         SamID=paste(Enc2, ShellSpecies, sep="."))%>%
  group_by(Enc2,TreatA,SamID) %>% 
  summarize(TShellSurArea.cm2=sum(ShellSurArea.mm2, na.rm=T)*.01)
#the ones with 0 I didn't measure height on them b/c not needed

## mussel biomass in each enclosure using tutorial from van Ee et al. (2020)
#building a taxa table to match with known regressions
tax.table<-data.frame(Family="Unionidae",
                      Tribe=c("Amblemini","Lampsilini"),
                      Genus=c("AMB","ACT"),
                      Genus2=c("Amblema","Actinonaias"),
                      Species=c("plicata","ligamentina"))
#getting raw data and covering for loss mussels (so NA with length etc. )
MusBiomassE<-MusselData %>% left_join(TreatENC) %>%
  mutate(Length_MM=case_when(is.na(L)&Genus=="AMB"~92.4,
                             is.na(L)&Genus=="ACT"~111.6,
                             Genus=="AMB"~ L,
                             Genus=='ACT'~ L),
         BiomassE=case_when(Notes=="NotRecovered"&Genus=="AMB"~3.47,
                            Notes=="NotRecovered"&Genus=="ACT"~7.64,
                            Genus=="AMB"~ 1.57e-6*L^3.21,
                            Genus=='ACT'~ 4.7e-7*L^3.51)) %>%
  dplyr::select(-W,-H,-TreatF,-Enclosure,-Exposure,-Notes) %>% left_join(tax.table)

# use MusBiomassE with length tutorial; results are saved as csv
MusBio<-read.csv("./data/EncMussBio.csv") %>%
  group_by(Enc2, Treatment, Genus) %>%
  dplyr::summarize(sumBM=sum(STDM_g, na.rm=T), 
            meanBM=mean(STDM_g,na.rm=T), 
            sdBM=sd(STDM_g, na.rm=T)) %>% 
  mutate(MusselBiomass.g.m2=sumBM/0.25)
MusBioG<-MusBio %>% dplyr::select(-sumBM, -meanBM, -sdBM)%>%
  group_by(Enc2) %>%
  spread(Genus, MusselBiomass.g.m2) 
MusBioS<-MusBio %>% ungroup() %>% group_by(Enc2) %>%
  dplyr::summarise(sumBMall.g.m2=sum(sumBM)/0.25)

# Enclosure Invertebrates ---------
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
# unique(EnclInv$Taxa)[order(unique(EnclInv$Taxa))] #check taxa for correct spelling

# Calculate area sampled 
BaskSamp<-read.csv("../FEn17/FEn17_data/SubSamFEn17OK.csv", stringsAsFactors = F) %>%
  mutate(WeekG=case_when(Week==8~"w09",
                         Week==4~"w04",
                         Week==12 ~"w12")) %>% 
  left_join(TreatENC, by="Enclosure") %>% 
  filter(Enclosure!="WC" & Type.x=="Insect") %>%
  select(Enclosure, Enc2, Basket., WeekG) %>%
  mutate(TEid=paste(WeekG, Enc2, sep="")) %>% filter(!duplicated(TEid))
#Some filters didn't have baskets associated with them
#assumed 3 baskets
BaskSamp[is.na(BaskSamp$Basket.),"Basket."]<-3

# macroinvertebrate summary in enclosures
ECounts<-EnclInv %>% group_by(TEid, Enc, Week) %>% 
  summarize(Nsam=n(),
            richness=length(unique(Taxa))) %>% 
  left_join(BaskSamp) %>%
  mutate(AreaSamp=(Basket.*.03315),
         InvertDensity.npm2=Nsam/AreaSamp) %>% 
  full_join(TreatENC) %>% 
  left_join(MusBioS) %>%
  replace(is.na(.),0)

# enclosure data community matrix in density; units are n/m2
E.ComDens<-EnclInv %>% group_by(TEid,Enc, Week, Taxa) %>% 
  summarize(N=n()) %>%
  left_join(ECounts) %>% mutate(IndDensity.m2=N/AreaSamp)%>%
  select(Enc,TreatA, Week,TEid, Taxa, IndDensity.m2) %>% 
  mutate(id=1:n()) %>% filter(Taxa!="MYSTERY") %>%
  spread(Taxa,IndDensity.m2) %>% group_by(TEid, Enc, TreatA) %>%
  summarise_if(is.numeric,list(~mean(.[which(!is.na(.))]))) %>%
  select(-id) %>% replace(is.na(.),0) %>%
  left_join(ECounts) %>% 
  select(-Nsam,-richness,-Enclosure,-Enc,-WeekG,-Basket., -AreaSamp,
         -InvertDensity.npm2, -TreatF)%>%
  select(TEid,Enc2,Week,TreatA,Type,Spp,sumBMall.g.m2, everything())

# enclosure data community matrix in relative abundance; units are %
E.ComPer<-EnclInv %>% group_by(TEid,Enc, Week, Taxa) %>% 
  filter(Taxa!="MYSTERY", Taxa!="MysteryTricoptera", Taxa!="Col.Myst",
         Taxa!="Tri.Pupa") %>% 
  summarize(N=n()) %>%  mutate(id=1:n()) %>% 
  spread(Taxa,N) %>% group_by(TEid)%>%
  summarise_if(is.numeric,list(~mean(.[which(!is.na(.))]))) %>%
  select(-id) %>% replace(is.na(.),0) %>%
  mutate(rowsum=rowSums(.[,-1])) %>% group_by(TEid) %>%
  summarise_if(is.numeric,list(~./rowsum)) %>%
  left_join(ECounts) %>% 
  select(-rowsum, -Nsam,-richness,-Enclosure,-Enc,-WeekG,-Basket., -AreaSamp,
         -InvertDensity.npm2, -TreatF)%>%
  select(TEid,Enc2,Week,TreatA,Type,Spp,sumBMall.g.m2, everything())

# Enclosure Biomass ===========
#head(EnclInv)
EBiomass<-EnclInv %>% left_join(BMTaxaTable) %>% filter(Order!="misc")%>% 
  rowwise()%>% 
  mutate(Length.mm=Length.cm*10,
         BM.mg=mean(biomass(Family, Order, Length.mm)))

# enclosure data biomass matrix in density; units are mg/m2
E.BioDens<- EBiomass %>% group_by(TEid, Taxa) %>%
  summarize(TBM.mg=sum(BM.mg)) %>%
  left_join(ECounts) %>% mutate(BM.mg.m2=TBM.mg/AreaSamp) %>%
  select(TEid, Taxa, BM.mg.m2) %>% group_by(TEid) %>%
  spread(Taxa, BM.mg.m2) %>%
  left_join(ECounts) %>% 
  select(-Nsam,-richness,-Enclosure,-Enc,-WeekG,-Basket., -AreaSamp,
         -InvertDensity.npm2, -TreatF)%>%
  select(TEid,Enc2,Week,TreatA,Type,Spp,sumBMall.g.m2, everything())

# enclosure data biomass matrix in relative abundance; units are % mg
E.BioPer<- EBiomass %>% group_by(TEid, Taxa) %>%
  summarize(TBM.mg=sum(BM.mg)) %>% group_by(TEid) %>%
  replace(is.na(.),0) %>% spread(Taxa,TBM.mg) %>% ungroup() %>%
  mutate(rowsum=rowSums(.[,-1], na.rm=T))%>% group_by(TEid)%>%
  summarise_if(is.numeric,list(~./rowsum)) %>% 
  left_join(ECounts) %>% 
  select(-rowsum, -Nsam,-richness,-Enclosure,-Enc,-WeekG,-Basket., -AreaSamp,
         -InvertDensity.npm2, -TreatF)%>%
  select(TEid,Enc2,Week,TreatA,Type,Spp,sumBMall.g.m2, everything())

# Shell Invertebrates -----------
#shell invertebrates sampled by collecting mussels or sham 
#and putting them in separate buckets based on their speciese
shellInv<-read_excel("../FEn17/FEn17_data/ShellInv.xlsx")
#long list of macroinvertebrate lengths
SInv<-shellInv %>% filter(Taxa!="Spider") %>% 
  select(Enc2, ShellSpecies, SamType, Taxa, Length.cm)

#macroinvertebrate summary on shells
SCounts<-SInv %>% group_by(Enc2, ShellSpecies) %>% 
  summarize(Nsam=n(),
            richness=length(unique(Taxa)))%>%
  mutate(SamID=paste(Enc2,ShellSpecies, sep="."))%>% full_join(MDat) %>%
  mutate(InvertDensity.npcm2=Nsam/TShellSurArea.cm2) %>% full_join(TreatENC) %>% 
  filter(!is.na(InvertDensity.npcm2))

# shell data community matrix in density; units are n/cm2
S.ComDens<-SInv %>% group_by(Enc2, ShellSpecies, Taxa) %>% 
  mutate(SamID=paste(Enc2,ShellSpecies, sep=".")) %>% summarize(N=n()) %>%
  left_join(SCounts) %>% mutate(IndDensity.cm2=N/TShellSurArea.cm2,
                                SamID=paste(Enc2, ShellSpecies, sep=".")) %>%
  select(Enc2, ShellSpecies,SamID, Taxa, IndDensity.cm2) %>% 
  mutate(id=1:n()) %>% 
  spread(Taxa,IndDensity.cm2) %>% group_by(SamID)%>%
  summarise_if(is.numeric,list(~mean(.[which(!is.na(.))]))) %>%
  replace(is.na(.),0)%>%
  left_join(SCounts)%>% 
  select(-id,-Nsam,-richness,-Enclosure,-TShellSurArea.cm2,
         -InvertDensity.npcm2,-TreatF)%>%
  mutate(SheType=paste(ShellSpecies,Type))%>%
  select(SamID,Enc2,ShellSpecies,TreatA,Type,Spp,SheType, everything())

# shell data community matrix in relative abundance; units are %
S.ComPer<-SInv %>% group_by(Enc2, ShellSpecies, Taxa) %>% 
  filter(Taxa!="MYSTERY", Taxa!="MysteryTricoptera", Taxa!="Col.Myst",
         Taxa!="Tri.Pupa", Taxa!="Dip.Adult", Taxa!="Orthoptera") %>% 
  summarize(N=n()) %>%  mutate(id=1:n(),
                               SamID=paste(Enc2,ShellSpecies, sep=".")) %>% 
  spread(Taxa,N) %>% group_by(SamID)%>%
  summarise_if(is.numeric,list(~mean(.[which(!is.na(.))]))) %>%
  select(-id) %>% replace(is.na(.),0) %>%
  mutate(rowsum=rowSums(.[,-1])) %>% ungroup()%>% group_by(SamID)%>%
  summarise_if(is.numeric,list(~./rowsum)) %>%
  left_join(SCounts)%>% 
  select(-rowsum, -Nsam,-richness,-Enclosure,-TShellSurArea.cm2,
         -InvertDensity.npcm2,-TreatF)%>%
  select(SamID,Enc2,ShellSpecies,TreatA,Type,Spp, everything())

# Shell Biomass ==============
head(SInv)
SBiomass<-SInv %>% group_by(Enc2, ShellSpecies, Taxa) %>% 
  mutate(SamID=paste(Enc2,ShellSpecies, sep=".")) %>%
  left_join(BMTaxaTable) %>% filter(Order!="misc")%>% 
  rowwise()%>% 
  mutate(Length.mm=Length.cm*10,
         BM.mg=mean(biomass(Family, Order, Length.mm)))

# shell data biomass matrix in density; units are mg/ m2
S.BioDens<- SBiomass %>% group_by(SamID, Taxa) %>%
  summarize(TBM.mg=sum(BM.mg))%>%
  left_join(SCounts) %>% mutate(BM.mg.cm2=TBM.mg/TShellSurArea.cm2)%>%
  select(SamID, Taxa, BM.mg.cm2) %>% group_by(SamID) %>%
  spread(Taxa, BM.mg.cm2) %>%
  left_join(SCounts)%>% 
  select(-Nsam,-richness,-Enclosure,-TShellSurArea.cm2,-InvertDensity.npcm2,-TreatF)%>%
  select(SamID,Enc2,ShellSpecies,TreatA,Type,Spp, everything())

# shell data biomass matrix in relative abundance; units are % mg
S.BioPer<- SBiomass %>% group_by(SamID, Taxa) %>%
  summarize(TBM.mg=sum(BM.mg)) %>% group_by(SamID) %>%
  replace(is.na(.),0) %>% spread(Taxa,TBM.mg) %>% ungroup() %>%
  mutate(rowsum=rowSums(.[,-1], na.rm=T))%>% group_by(SamID)%>%
  summarise_if(is.numeric,list(~./rowsum)) %>% 
  left_join(SCounts) %>% 
  select(-rowsum,-Nsam,-richness,-Enclosure,-TShellSurArea.cm2,
         -InvertDensity.npcm2,-TreatF)%>%
  select(SamID,Enc2,ShellSpecies,TreatA,Type,Spp, everything())
