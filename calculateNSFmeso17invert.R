library(readxl);library(tidyverse)
counts<-read_xlsx('~/Desktop/NSF2017 Invert Janell raw.xlsx', sheet=2)
head(counts)

len17<-read_xlsx('~/Desktop/NSF2017 Invert Janell raw.xlsx', sheet=1) %>%
  select(Habitat, Week, Family, Length) %>% rename(Taxa=Family) %>%
  mutate(SampleID=paste(Habitat, paste("wk",Week), sep="-"),
         Length.mm=Length*10)
head(len17)
unique(len17$Taxa)

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

raw_biomass<-len17 %>% filter(Taxa!="Other") %>%
  left_join(BMTaxaTable) %>%
  rowwise() %>% #throw each row into the mutate
  mutate(BM.mg=mean(biomass(Family, Order, Length.mm)))

sum_biomass<-raw_biomass %>% group_by(Habitat, Week, SampleID, Family) %>%
  summarize(mean(BM.mg), sd(BM.mg), mean(Length.mm))
