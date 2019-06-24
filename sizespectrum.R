library(tidyverse) #libraries
source('InvDataframes.R') #bring in the invertebrate data
melt2016lengthsN<- melt2016lengths %>% mutate(Length.mm = Length.cm*10,
                                              Reach=substr(Sample, 1,5),
                                              SamplingSeason=substr(Sample, 7,11)) %>%
  dplyr::select(Sample, Reach, SamplingSeason, Taxa, Length.mm, Family, Order, Class,BM.mg)
t2 <- Lengths2016 %>% gather(Taxa, Length, -Sample) %>%
  filter(!is.na(Length)) %>% group_by(Sample, Taxa) %>% tally()

t1<-Count2016 %>% ungroup() %>%
  dplyr::select(-Location, -Site, -Treatment,-Reach, -SamplingSeason) %>%
  gather(Taxa, Count, -Sample) %>% filter(Count!=0)

#not all macroinvertebrates were measured for length, but most were. 
#here are my trouble problems
key<-full_join(t1,t2) %>% mutate(DIF=Count-n,
                            perDif=DIF/Count*100) %>%
  arrange(desc(DIF)) %>% filter(!is.na(Count), !is.na(n)) %>%
  #taxa want to exclude
  filter(Taxa!="Dip.Other", Taxa!="Other.other", Taxa!="Hirudinea.leach",
         Taxa!="Cladocera.miscD", Taxa!="Copepoda.miscC", Taxa!="Ostracoda") #regressions are not good


melt2016lengthsN %>% group_by(Sample,Taxa)
key %>% filter(is.na(DIF))
key<-key %>% mutate(mult=Count %/%n,
                                  samp=Count%%n)

repData<-function(Taxa, Samp, Lengths) {
 #pull key for that taxa
 keyI<-key[key$Taxa==Taxa & 
           key$Sample==Samp &
           !is.na(key$Sample==Samp), ]
 newLengths<-c(rep(Lengths,as.numeric(keyI[,7])), sample(Lengths,as.numeric(keyI[,8])))
 tibble(Sample=Samp,Taxa=Taxa, Lengths.n=newLengths)
}
fullD<-NULL
for(k in 1:nrow(key)){
  sub<-Lengths2016 %>% gather(Taxa, Length.cm, -Sample) %>% 
    filter(Sample==as.character(key[k,"Sample"]) & Taxa==as.character(key[k, "Taxa"]))
  testing<-repData(unique(sub$Taxa), unique(sub$Sample), sub$Length.cm)
  #print(paste(k,paste(key[k,"Sample"],key[k,"Taxa"], sep="-")))
  fullD<-rbind(fullD, testing)
}
#its putting a lot of NA in for some reason
fullD %>% group_by(Sample, Taxa) %>% 
  filter(!is.na(Lengths.n)) %>% tally() %>% arrange(desc(n))

t3<-melt2016lengths %>% filter(!(Sample %in% fullD$Sample)) %>%
  bind_rows(fullD) %>% filter(!is.na(Lengths.n)) %>% group_by(Sample,Taxa) %>%
  tally(name="reppedN")
testfin<-full_join(key,t3) %>% mutate(check=Count-reppedN) %>% arrange(check) 
mean(testfin$check) #good enough


test1<-melt2016lengthsN %>% filter(Taxa=="Tri.Hydropsychidae", Sample=="LY-MR Oct16 S2")
hist(test1$Length.mm)
mean(test1$Length.mm);sd(test1$Length.mm)
test2<-rnorm(10,mean(test1$Length.mm), sd(test1$Length.mm))
replicate(4, test2)
test3<-rlnorm(100,log10(mean(test1$Length.mm)), log10(sd(test1$Length.mm)))
library(fitdistrplus)
temp1<-fitdist(test1$Length.mm, distr = "gamma", method = "mle")
summary(temp1)
hist(rgamma(100, temp1[1]$estimate[1], temp1[1]$estimate[2]))
hist(test1$Length.mm)
test3<-rep(test1$Length.mm, length.out=495)
hist(test3)
mean(test3);sd(test3)
hist(test1$Length.mm)
hist(c(test1$Length.mm,test2))
  
size.spec.test<- rbind(melt2015lengths, melt2016lengthsN) %>% 
  mutate(Treatment=substr(Sample, 4,5)) %>%
  filter(!is.na(BM.mg)) %>% mutate(BM.g=BM.mg/1000)
head(size.spec.test)
# need to account for samples where count !=measured
# need to account for different sample sizes

ggplot(size.spec.test, aes(x=Treatment, y=Length.mm, fill=SamplingSeason))+
  geom_violin()
#removing NaN results for BM.mg removes a lot of large things from NM

# make sure that range of breaks is
# greater than range of data
test1 <- min(breaks) < min(size.spec.test$BM.g) & max(breaks) > max(size.spec.test$BM.g)

if(test1 == TRUE){
  # applies bin_and_center() to the dataset, after subsetting by "site" variable
  binned <- ddply(size.spec.test, .(Sample),
                  bin_and_center,
                  var = "BM.g",
                  breaks = breaks)
}else{
  "Breaks does not cover range of data"
}

binSamp<-binned %>% left_join(fieldSiteKey, by="Sample")
# save results
write.csv(binned, "results/binned_size_spectra.csv",
          row.names = FALSE)
ggplot(binSamp, aes(x=log_mids, y=log_count_corrected, color=Reach))+
  geom_point(size=.5)+
  geom_smooth(formula=y~x^2) + facet_wrap(~Treatment)
