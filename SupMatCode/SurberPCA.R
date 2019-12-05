# Using a PCA to look at variation within the macroinvertebrate communities 
# BETWEEN Surber samples from the field surveys
#load invdataframes.R and envdataframes.R 
source('InvDataframes.R'); source('EnvDataframes.R')

#check variation among surber samples at the field level
SurvComMatF16<-samp2016count %>% 
  spread(Taxa, Count) %>% 
  filter(SamplingSeason == "Fall2016")%>% # only considering Fall 2016
  dplyr::select(-SamplingSeason, # pointless since only one sampling season
                -Hem.Gerridae,  -Dip.Chaoboridae,#only in 1
                -Dip.Orthorrhaphous, -Dip.Thaumaleidae, -Lep.Pyralidae, #only in 1
                -Biv.Unionidae) %>% #don't care about these
  replace(is.na(.),0) %>% #replace na with 0 for vegan
  mutate_if(is.numeric,log1p) %>%
  mutate(Reach=substr(Sample, 1,5),
         Site=substr(Sample, 1,2),
         Treat=substr(Sample, 4,5))%>%
  select(Reach, Site, Treat, Sample, everything())

#checking the dataframe has what I need
names(SurvComMatF16)
PCAreach<-prcomp(SurvComMatF16[,-c(1:4)]) #run a pca only on species data
PCAgraphD<-fortify(PCAreach) %>% cbind(SurvComMatF16[,1:4]) #makes a good plotting dataframe

library(ggrepel)
#nice plot of PCA to include in supplementary material
ggplot()+
  geom_point(data=PCAgraphD, aes(x=PC1, y=PC2, shape=Treat), size=2)+ #points for each sample
  geom_text_repel(data=PCAgraphD[PCAgraphD$PC1>=4 | #label the extream points
                                 PCAgraphD$PC2<=-3 |
                                 PCAgraphD$PC2>=3,], 
                  aes(x=PC1, y=PC2, label=Site),size=3)+
  ylab("PC2 (15.2%)")+xlab("PC1 (35.6%)")+theme_bw()
ggsave("./Figures/S2PCAsurbers.tiff",width=3.5, height=3)

nrow(SurvComMatF16)
5*7*2
#conclusion: LY MR has high variation, 
# as does GL NM, K2 NM and KS MR
