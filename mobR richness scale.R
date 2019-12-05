library(devtools)
library(mobr)
source('InvDataframes.R') #invertebrate data
source('EnvDataframes.R') #environment data

SurvComMatF162<-samp2016count %>% spread(Taxa, Count) %>% 
  filter(SamplingSeason == "Fall2016")%>%
  dplyr::select(-SamplingSeason, 
                -Hem.Gerridae,
                -Dip.Chaoboridae,#only in 1
                -Dip.Orthorrhaphous, -Dip.Thaumaleidae, -Lep.Pyralidae, #only in 1
                #-Neu.Sisyridae,-Col.Hydrophilidae,-Dip.Tupulidae, #only in 2
                #-Tri.Odontoceridae, -`FLAT WORMS`, #only in 2 samples
                -Biv.Unionidae) %>% #don't care about these
  replace(is.na(.),0) %>% #mutate_if(is.numeric,log1p) %>%
  mutate(Reach=substr(Sample, 1,5),
         Site=substr(Sample, 1,2),
         Treat=substr(Sample, 4,5))%>%
  select(Reach, Site, Treat, Sample, everything())

#using field.com and field.env from CCA.R
fmob<-FieldSpData[-15,]
fspat<-cbind(fmob, fmob@coords)
fspat2<- left_join(SurvComMatF162[,1:4], fspat@data)
inv_mob_2<-make_mob_in(SurvComMatF162[,-c(1:4)], fspat2,
                       coord_names=c("Longitude","Latitude"))
inv_mob_2
#spatial, sample-based rarefaction sSBR
plot_rarefaction(inv_mob_2,'Treat','spat')

#individual rarefaction (only reflect SAD)
plot_rarefaction(inv_mob_2, 'Treat','indiv', pooled=F, leg_loc='topright')
plot_rarefaction(inv_mob_2, 'Treat','indiv', pooled=T)

#directly examine SAD
plot_abu(inv_mob_2, 'Treat', type='rad', pooled=T, log='x')

inv_stats<-get_mob_stats(inv_mob_2, group_var='Treat', n_perm=500)

plot(inv_stats, 'S')
plot(inv_stats, 'N')
plot(inv_stats, 'S_n')
plot(inv_stats, 'S_PIE')

inv_deltaS <- get_delta_stats(inv_mob_2, 'Treat',
                              ref_group = 'NM',
                              type='discrete', log_scale=T,
                              n_perm=200)
plot(inv_deltaS, 'MR','NM', display='rarefaction')
plot(inv_deltaS, 'MR','NM', display='delta S')

plot(inv_deltaS, 'MR','NM', display='ddelta S')

par(mfrow=c(1,2))
overlap_effects(inv_deltaS, 'MR', display='raw', leg_loc='bottomleft')
overlap_effects(inv_deltaS, 'MR', display='stacked', 
                prop=T, leg_loc=NA)

### enclosure ###
#using enc.com and enc.env from CCA.R, 
# Enclosure Raster in lab meeting checks.R
emob<-left_join(enc.env, EnclosureRaster[,c(1,3,4,11)])
inv_mob_enc<-make_mob_in(enc.com, emob,
                       coord_names=c("xc","yc"))
#spatial, sample-based rarefaction sSBR
plot_rarefaction(inv_mob_enc,'TreatA','spat',
                 leg_loc = 'bottomright')

#individual rarefaction (only reflect SAD)
plot_rarefaction(inv_mob_enc, 'TreatA','indiv', pooled=F, 
                 leg_loc='bottomright')
plot_rarefaction(inv_mob_enc, 'TreatA','indiv', pooled=T,
                 leg_loc='bottomright')

#directly examine SAD
plot_abu(inv_mob_enc, 'TreatA', type='rad', pooled=T, log='x')

inv_stats_enc<-get_mob_stats(inv_mob_enc, group_var='TreatA', n_perm=500)

plot(inv_stats_enc, 'S')
plot(inv_stats_enc, 'N')
plot(inv_stats_enc, 'S_n')
plot(inv_stats_enc, 'S_PIE')

inv_deltaS <- get_delta_stats(inv_mob_enc, 'TreatA',
                              ref_group = 'CTRL',
                              type='discrete', log_scale=T,
                              n_perm=20)
plot(inv_deltaS, 'MR','NM', display='rarefaction')
plot(inv_deltaS, 'ACTL','AMBL', display='delta S')

plot(inv_deltaS, 'MR','NM', display='ddelta S')

par(mfrow=c(1,2))
overlap_effects(inv_deltaS, 'MR', display='raw', leg_loc='bottomleft')
overlap_effects(inv_deltaS, 'MR', display='stacked', 
                prop=T, leg_loc=NA)

