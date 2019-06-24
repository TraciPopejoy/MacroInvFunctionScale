library(tidyverse) #libraries
source('InvDataframes.R') #bring in the invertebrate data
#invertebrates were subsampled for lengths (100 or 10%, whichever higher)
#so we need to increase the data for some taxa
#building a 'key' to determine which taxa should be increased
#relying on counts to be more accurate than lengths
# determining data deficiencies ============
key2 <- Lengths2016 %>% bind_rows(Lengths2015[,-c(2:5, 50:52)]) %>%
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
         Taxa!="Cladocera.miscD", Taxa!="Copepoda.miscC", Taxa!="Ostracoda") #regressions are not good

View(key %>% filter(is.na(DIF))) #taxa I need to check or need to drop

# replicating subsampled data =================
set.seed(051520) #because I'm using sample() need to have random number set
#taxa I actually want to increase using the following code
keyF<-key %>% mutate(mult=Count %/%n, samp=Count%%n) %>% filter(DIF>=1) 

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
  sub<-Lengths2016 %>%  bind_rows(Lengths2015[,-c(2:5, 50:52)]) %>% #all lengths
    gather(Taxa, Length.cm, -Sample) %>% 
    filter(Sample==as.character(keyF[k,"Sample"]) & Taxa==as.character(keyF[k, "Taxa"]))
  #apply the function to get repeated data
  testing<-repData(unique(sub$Taxa), unique(sub$Sample), sub$Length.cm)
  #print(paste(k,paste(keyF[k,"Sample"],keyF[k,"Taxa"], sep="-")))
  repDF<-rbind(repDF, testing) #bind it to the bigger data frame
}
#its putting a lot of NA in for some reason
repDF %>% group_by(Sample, Taxa) %>% 
  filter(!is.na(Lengths.n)) %>% tally() %>% arrange(desc(n))
#checking how well it actually increased the data compared to count
t3<-Lengths2016 %>% bind_rows(Lengths2015[,-c(2:5, 50:52)]) %>%
  gather(Taxa, Lengths.n, -Sample) %>%
  filter(!(Sample %in% repDF$Sample)) %>%
  bind_rows(repDF) %>% filter(!is.na(Lengths.n)) %>% group_by(Sample,Taxa) %>%
  tally(name="reppedN")
testfin<-full_join(keyF,t3) %>% mutate(check=Count-reppedN) %>% arrange(check) 
mean(testfin$check, na.rm=T) #good enough

#so two problems so far 1. NA happening during for loop & 2. sample() not pulling enough samples

# building the long dataset ================
#need treatment, length.mm, biomass.g
SizeSpecData<-Lengths2016 %>% bind_rows(Lengths2015[,-c(2:5, 50:52)] %>% 
                                          group_by(Sample) %>%
                                          transmute_at(vars(-group_cols()),function(x){x/10})) %>%
  gather(Taxa, Lengths.n, -Sample) %>% #get data in one long frame
  filter(!(Sample %in% fullD$Sample)) %>%
  bind_rows(fullD) %>% filter(!is.na(Lengths.n)) %>% ungroup() %>%
  left_join(BMTaxaTable) %>% filter(Order!="misc") %>% 
  rowwise() %>% 
  #calculate biomass from length data
  #will get NAN for taxa without regressions and for lengths outside size range
  mutate(Length.mm=Lengths.n*10, 
         BM.mg=mean(biomass(Family, Order, Length.mm)),
         BM.g=BM.mg/1000) %>%
  left_join(fieldSiteKey) %>%
  filter(!is.na(BM.mg)) #removing NAN from no regression
head(SizeSpecData)
summary(SizeSpecData$Length.mm) #checking it worked
# need to account for different sample sizes; need to average at reach level

ggplot(SizeSpecData, aes(x=Treatment, y=Length.mm, fill=SampSeaF))+
  geom_violin()+scale_y_log10()
ggplot(SizeSpecData, aes(x=Treatment, y=BM.g, fill=SampSeaF))+
  geom_violin()+scale_y_log10()

# Size Spectrum Analysis - breaks ===========
#code modified from https://github.com/Jpomz/mining-size-sprectra-Freshwater-Biology-accepted
#citation: Pomeranz, Warburton and Harding. 2018. Anthropogenic mining alters macroinvertebrate size spectra in streams. Freshwater Biology.
#### have to load plyr, which means with biomass() won't work
# plyr interferes with rowwise I think
library(plyr)
breaks <- 2^seq(-25,2)
# make sure that range of breaks is greater than range of data
test1 <- min(breaks) < min(SizeSpecData$BM.g) & max(breaks) > max(SizeSpecData$BM.g)

if(test1 == TRUE){
  # applies bin_and_center() to the dataset, after subsetting by "site" variable
  binned <- ddply(SizeSpecData, .(Sample),
                  bin_and_center,
                  var = "BM.g",
                  breaks = breaks)
}else{
  "Breaks does not cover range of data"
}

binSamp<-binned %>% left_join(fieldSiteKey, by="Sample")
# save results
write.csv(binned, "binned_size_spectra.csv",
          row.names = FALSE)
ggplot(binSamp, aes(x=log_mids, y=log_count_corrected))+
  geom_point(size=.5)+
  geom_smooth(formula=y~x^2, level=.5) + facet_wrap(~Treatment+Season)+
  scale_x_continuous(name="Log Mass")+
  scale_y_continuous(name="Log Count")+
  theme_bw()
