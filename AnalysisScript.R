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
#use this table to supplement other data sources 
#(mainly Le Roy Poff 2006 table)

##### MACROINVERTEBRATE DATA FROM OTHER SCRIPTS #####
#Field Data
head(FieldComMat)
head(FieldTraits)

#Enclosure Data
head(EnclComMat)
head(EnclTraits)

#Shell Data
head(ShellComMat)
head(ShellTraits)

##### AlphaFD.R from DecomposingFD #####
source("../DecomposingFD/R/AlphaFD.R")

max(dist(EnclTraits[,-1]))
ETdist<-dist(EnclTraits[,-1])/max(dist(EnclTraits[,-1]))
FTD(ETdist)
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

#### Graphs ####
library(colorspace)
CP<-diverge_hcl(5, h=c(180,70), c = 100, l = c(50, 90), power = 1)
CP[3]<-"lightgrey"
show_col(CP)

FCounts<-meanCount2016 %>% 
  select(Reach, SamplingSeason, InvDensity.npm2, richness, Treatment) %>%
  mutate(SeasF=factor(SamplingSeason, levels=c("Summer2016","Fall2016")))
head(FCounts)
ECounts<-ECounts %>% 
  mutate(TreFac=factor(TreatA, levels=c("CTRL","ACTL","ACTS","AMBL","AMBS")))
head(ECounts)
SCounts<-SCounts %>%
  filter(!is.na(ShellSpecies)) %>%
  mutate(ShID=paste(ShellSpecies,Type, sep="x"))
head(SCounts)

fdens<-ggplot(FCounts, aes(x=SeasF, y=InvDensity.npm2, fill=Treatment))+
  geom_boxplot(aes(group=interaction(SamplingSeason,Treatment)))+
  scale_fill_manual(name="River Reach\nTreatment",
                    labels=c("Mussel Bed","Control"),
                    values = CP[c(1,3)])+
  scale_y_continuous(name=expression("Macroinvertebrates per m " ^ -2))+
  scale_x_discrete(name="Sampling Period", 
                   labels=c("Summer 2016","Fall 2016"))+
  theme(legend.justification=c(1,1),legend.position = c(0.575,1))
edenss<-ggplot(ECounts, aes(x=Week, y=InvertDensity.npm2, fill=TreFac))+
  geom_boxplot(aes(group=interaction(Week,TreFac)))+
  scale_fill_manual(name="Enclosure\nTreatment",
                    values = CP[c(3,1,2,5,4)])+
  scale_y_continuous(name=expression("Macroinvertebrates per m " ^ -2))+
  scale_x_discrete(name="Sampling Period", 
                   labels=c("Sept. 2017","Oct. 2017"))+
  theme(legend.justification=c(1,1),legend.position = c(0.5,1))
sdens<-ggplot(SCounts, aes(x=TreatA, y=InvertDensity.npcm2/0.0001, fill=ShID))+
  geom_boxplot(aes(group=interaction(TreatA,ShID))) +
  scale_fill_manual(name="Shells",
                    values = CP[c(1,2,5,4)],
                    labels=c("Live ACT", "Sham ACT","Live AMB","Sham AMB"))+
  scale_y_continuous(name=expression("Macroinvertebrates per m " ^ -2))+
  scale_x_discrete(name="Enclosure Treatment")+
  theme(legend.justification=c(1,1),legend.position = c(1,1))
,
        legend.direction="horizontal")+guides(fill=guide_legend(nrow=2))
library(cowplot)
plot_grid(sdens,edenss,fdens, ncol=3,labels = "AUTO")
ggsave("InvertDensity.tiff", dpi=600, width = 10, height = 3.5, units="in")

