#libraries
library(tidyverse)

### feeding groups table from LeRoy Poff 2008 ###
ffg<-read_csv("feedingtraits.csv")
FTaxaTable<-read_excel("Macroinv Power Law Coeffs TBP.XLSX", sheet=2)
taxa<-sort(unique(FTaxaTable$Family[FTaxaTable$Family!="misc"]))
FFGunk<-ffg %>% right_join(FTaxaTable[FTaxaTable$Family!="misc",]) %>%
  select(Family, Order, Feed_mode_prim, Feed_mode_sec, Feed_mode_comments) %>%
  distinct()
traitDIST<-FFGunk %>% group_by(Family) %>% filter(!is.na(Feed_mode_sec)) %>%
  summarize(filt=sum(str_count(Feed_mode_sec, "filterer"), na.rm=T),
            gath=sum(str_count(Feed_mode_sec, "gatherer"), na.rm=T),
            other=sum(str_count(Feed_mode_sec, "Other"), na.rm=T),
            parasite=sum(str_count(Feed_mode_sec, "Parasite"), na.rm=T),
            pierc=sum(str_count(Feed_mode_sec, "Piercer"), na.rm=T),
            pred=sum(str_count(Feed_mode_sec, "Predator"), na.rm=T),
            scrap=sum(str_count(Feed_mode_sec, "Scraper"), na.rm=T),
            shred=sum(str_count(Feed_mode_sec, "Shredder"), na.rm=T))
#use this table to supplement other data sources (mainly Le Roy Poff 2006 table)

##### MACROINVERTEBRATE DATA FROM OTHER SCRIPTS #####
#Field Data
head(FieldComMat)
head(FieldTraits)

#Enclosure Data

#Shell Data
head(ShellComMat)
head(ShellTraits)

##### AlphaFD.R from DecomposingFD #####
FTD()
## matrix or dist containing species functional distances (tdmat)
## tdmat should be rescaled to 0<=dij<=1
## Hill coefficient q
## weights are proportions of each species in each plot
FTD.comm()
## wrapper for the above function across multiple communities
## requires distance matrix w/ all species, scaled to 0-1
## community data matrix (communities x species)
## with species in same order as distance matrix
## if abund=T, values in community data matrix are treated as abundances
## otherwise, all are converted to presence-absence
## if match.names=T, the code will match species names across the
## trait distance matrix and comm data matrix and rearrange the latter
