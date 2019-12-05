library(tidyverse) #libraries
source('InvDataframes.R') #bring in the invertebrate data
#invertebrates were subsampled for lengths (100 or 10%, whichever higher)
#so we need to increase the data for some taxa
#building a 'key' to determine which taxa should be increased
#relying on counts to be more accurate than lengths
# determining data deficiencies ============
key2 <- Lengths2016 %>% bind_rows(Lengths2015[,-c(2:5, 56:58)]) %>%
  gather(Taxa, Length, -Sample) %>%
  filter(!is.na(Length)) %>% group_by(Sample, Taxa) %>% tally()

key1<-CountField %>% ungroup() %>%
  dplyr::select(-Location, -Site, -Treatment,-Reach, -SamplingSeason,
                -InvDensity.npm2, -richness) %>%
  gather(Taxa, Count, -Sample) %>% filter(Count!=0)

key<-full_join(key1,key2) %>% mutate(DIF=Count-n,
                            perDif=DIF/Count*100) %>%
  arrange(desc(DIF)) %>% #filter(!is.na(Count), !is.na(n)) %>%
  #taxa want to exclude
  filter(Taxa!="Dip.Other", Taxa!="Other.other", Taxa!="Hirudinea.leach", #rare
         Taxa!="Hydrozoans.miscH", #rare
         Taxa!="Cladocera.miscD", Taxa!="Copepoda.miscC", Taxa!="Ostracoda") #regressions are not good
check<-key %>% filter(is.na(DIF))
write.csv(check, file="invertLengthsFix.csv") #taxa I need to check or need to drop

# replicating subsampled data =================
set.seed(051520) #because I'm using sample() need to have random number set
#taxa I actually want to increase using the following code
keyF<-key %>% mutate(mult=Count %/%n, samp=Count%%n) %>% filter(DIF!=0) 

Lengths2015.cm<-Lengths2015[,-c(2:5, 56:58)] %>% 
  group_by(Sample) %>%
  transmute_at(vars(-group_cols()), function(x){x/10})

#function to replicate the data
repData<-function(Taxa, Samp, Lengths) {
 #pull key for that taxa
 keyI<-keyF[keyF$Taxa==Taxa & keyF$Sample==Samp &
            !is.na(keyF$Sample) & !is.na(keyF$Taxa), ] #includes mult & samp
 # mult is the number of full times the lengths can go into the count
 # will use rep() to multiply this number to get the full 'lengths'
 # samp is the remained after dividing count by lengths
 # will use sample() to subsample the lengths measured to fill out the count
 newLengths<-c(rep(Lengths,as.numeric(keyI[,7])), base::sample(Lengths,as.numeric(keyI[,8])))
 tibble(Sample=Samp,Taxa=Taxa, Lengths.n=newLengths) #give me a dataframe
}
repDF<-NULL #for repeated data
for(k in 1:nrow(keyF)){
  #pull out the lengths measured for that sample and taxa
  sub<-Lengths2016 %>%  bind_rows(Lengths2015.cm) %>% #all lengths
    gather(Taxa, Length, -Sample) %>% filter(!is.na(Length)) %>%
    filter(Sample==as.character(keyF[k,"Sample"]) & Taxa==as.character(keyF[k, "Taxa"]))
  #apply the function to get repeated data
  testing<-repData(unique(sub$Taxa), unique(sub$Sample), sub$Length)
  #print(paste(k,paste(keyF[k,"Sample"],keyF[k,"Taxa"], sep="-")))
  repDF<-rbind(repDF, testing) #bind it to the bigger data frame
}

repDF %>% group_by(Sample, Taxa) %>% 
  filter(!is.na(Lengths.n)) %>% tally() %>% arrange(desc(n))
#checking how well it actually increased the data compared to count
t3<-Lengths2016 %>% bind_rows(Lengths2015[,-c(2:5, 56:58)]) %>%
  gather(Taxa, Lengths.n, -Sample) %>%
  filter(!(Sample %in% repDF$Sample)) %>%
  bind_rows(repDF) %>% filter(!is.na(Lengths.n)) %>% group_by(Sample,Taxa) %>%
  tally(name="reppedN")
testfin<-full_join(keyF,t3) %>% mutate(check=Count-reppedN) %>% arrange(check) 
mean(testfin$check, na.rm=T) #have a length for each Count!
nrow(repDF)

gather(CountField, Taxa, count, -InvDensity.npm2, -richness, -Sample,
       -Location,-Site, -Treatment, -Reach, -SamplingSeason) %>%
  select(Sample,Taxa,count) %>% ungroup() %>%
  filter(Taxa!="Dip.Other", Taxa!="Other.other", Taxa!="Hirudinea.leach", #rare
         Taxa!="Hydrozoans.miscH", #rare
         Taxa!="Cladocera.miscD", Taxa!="Copepoda.miscC", Taxa!="Ostracoda") %>% 
  summarize(sum(count))


# building the long dataset ================
#need treatment, length.mm, biomass.g
SizeSpecData<-Lengths2016 %>% 
  bind_rows(Lengths2015.cm) %>% gather(Taxa, Lengths.n, -Sample) %>%
  anti_join(repDF[,-3]) %>% #remove samples that are in repDF
  bind_rows(repDF) %>% filter(!is.na(Lengths.n)) %>% ungroup() %>%
  left_join(BMTaxaTable) %>% filter(Order!="misc") %>% 
  rowwise() %>% 
  #calculate biomass from length data
  #will get NAN for taxa without regressions and for lengths outside size range
  mutate(Length.mm=Lengths.n*10, 
         BM.mg=mean(biomass(Family, Order, Length.mm)),
         BM.g=BM.mg/1000) %>%
  left_join(fieldSiteKey) %>%
  filter(!is.na(BM.mg), BM.mg!=0) #removing NAN from no regression
head(SizeSpecData)
#total count is 41077
#2755 lost because order has no BM regression
#2131 lost because out of regression size range
nrow(SizeSpecData) #38949
# need to account for different sample sizes; need to average at reach level

ggplot(SizeSpecData, aes(x=Treatment, y=Length.mm/10, fill=SampSeaF))+
  geom_violin()+scale_y_log10()
ggplot(SizeSpecData, aes(x=Treatment, y=BM.g, fill=SampSeaF))+
  geom_violin()+scale_y_log10()

# Size Spectrum Analysis - breaks ===========
#code modified from https://github.com/Jpomz/mining-size-sprectra-Freshwater-Biology-accepted
#citation: Pomeranz, Warburton and Harding. 2018. Anthropogenic mining alters macroinvertebrate size spectra in streams. Freshwater Biology.
#### have to load plyr, which means with biomass() won't work
# plyr interferes with rowwise I think
library(plyr)
breaks <- 2^seq(-25,1)
# make sure that range of breaks is greater than range of data
test1 <- min(breaks) < min(SizeSpecData$BM.g) & max(breaks) > max(SizeSpecData$BM.g)

if(test1 == TRUE){
  # applies bin_and_center() to the dataset, after subsetting by "site" variable
  binned <- ddply(SizeSpecData, .(FSamID), #reach is the aggregation level
                  bin_and_center,
                  var = "BM.g",
                  breaks = breaks)
}else{
  "Breaks does not cover range of data"
}

# save results
write.csv(binned, "binned_size_spectraNov5.csv",
          row.names = FALSE)

# 3_statistical_analyses
library(MuMIn);library(quantreg);library(ggplot2)

# read in binned size spectra data
binned <- read.csv("binned_size_spectraNov5.csv",
                   stringsAsFactors = FALSE)
# make a data frame with the appropriate info: Mussel biomass, year, season, huc12
SizeRegData<-binned %>% left_join(unique(fieldSiteKey[,-1])) %>% 
  left_join(FieldSpData@data) %>% 
  dplyr::select(-SamplingSeason,-Quadrat.n,-`USGS gage Identifier`,
                -SiteID,-HUC12,-HUC8)

# global models ####
# full quadratic
global.quadratic <- lm(log_count_corrected~ 
                         log_mids_center * Mussel.g.m2 +
                         I(log_mids_center^2) + 
                         I(log_mids_center^2):Mussel.g.m2, #+Year+Season,
                       data=SizeRegData,
                       na.action = "na.fail")
# full linear model
global.linear <- lm(log_count_corrected ~ 
                      log_mids_center+Mussel.g.m2,#+Year+HUC12num+Season,
                    data=SizeRegData,
                    na.action = "na.fail")
# compare quadratic and linear models
AIC(global.quadratic, global.linear)
# quadratic much better than linear
# move forward with quad model

# systematically test simplified quadratic models using MuMIn::dredge 
dredge.models.quad <- dredge(global.quadratic,
                        beta = "none",
                        extra = "R^2")

dredge.models.lin <- dredge(global.linear,
                             beta = "none",
                             extra = "R^2")

# table 2 for MS #### 
# table with results for all simplified models
sink("table_2Nov5.txt")
dredge.models.quad
sink()

# pick best model based on AICc using MuMIn::get.models
top.models <- get.models(dredge.models.quad,
                         subset = delta < 2)
# single top model selected
sink("top_models.txt")
top.models
sink()
# save top.model for making figure in script 4
saveRDS(top.models[[1]],
        "top_quadratic_model.RDS")


# Mrange ####
# Mrange ####
# calculate range of M values
mrange <- SizeRegData %>% 
  group_by(Reach,FSamID, Mussel.bm.g) %>%
  summarise(mrange = max(log_mids_center) - min(log_mids_center))
# save Mrange as csv
write.csv(mrange, "mrange_data.csv",
          row.names = FALSE)

# mrange ~ gradient linear model
m.mod <- lm(mrange ~ pca1, data = mrange)
summary(m.mod)

# Quantile regression
M.quant <- rq(log_mids~pca1, data = binned, tau = c(0.05, 0.95))
summary(M.quant)


