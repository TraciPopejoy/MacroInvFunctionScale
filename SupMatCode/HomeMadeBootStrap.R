# Bootstrap Surber Samples from field surveys ======================
# goal: sample one surber at reach level,
# 		run a CCA and pull relevant metrics
# 		make plots that describe variation in CCA results from variation at Surber level

# taxa match averaged field com data but each row is a surber sample
F.ComDens.surbers<-samp2016count %>%
  filter(!(Taxa %in% c(empty_col_field, tax_rare_field$taxa,
                       'Biv.Unionidae', #don't care about these
                       'Entom.Isotomidae', 'Other.other',
                       'Isopoda.miscI', 'Dip.Other'))) 

# For Loop for bootstrap -------------------------------------
set.seed(5432) #want a reproducible result
sumdataT<-NULL #create empty dataframe for CCA results (Constrained % etc.)
ANO<-NULL #create a different dataframe for ordistep reulsts
spDat1<-data.frame(rowname=unique(F.ComDens.surbers$Taxa)) #start a dataframe for sp on CCA1
spDat2<-data.frame(rowname=unique(F.ComDens.surbers$Taxa)) #start a dataframe for sp on CCA2
for(k in 1:1000) {
	#create a dataframe that has one surber sample for reach, randomly chosen
  surbsamp<-F.ComDens.surbers %>% 
  spread(Taxa, Count) %>% 
    filter(SamplingSeason == "Fall2016")%>%
    dplyr::select(-SamplingSeason) %>% #don't care about these
    replace(is.na(.),0) %>%
    mutate(Reach=substr(Sample, 1,5),
           Site=substr(Sample, 1,2),
           Treat=substr(Sample, 4,5))%>%
    dplyr::select(Reach, Site, Treat, Sample, everything()) %>% 
	group_by(Reach) %>%
    sample_n(1) #this does the random choice part
  # join this new dataframe with environment to pair species with environment	
  sampdf<-surbsamp %>% left_join(field.env, by="Reach")
  # run the CCA on this new dataframe
  cca.samp<-cca(sampdf[,-c(1:4,50:64)]~`Average of STDM (g.m-2)`+
                  Benthic_CHLA_MG.M2+ Discharge.cms+Dvar+
                  Condition(HUC12num), sampdf, scaling=2)
  ss<-summary(cca.samp) #save summary of CCA results
  #identify and sort which summary stats I want to keep (and keep track of which run it came from)
  sumdata<-data.frame(run=k, 
                      uniq=sampdf[1,4],
                      t(ss$biplot[,1]),
                      t(ss$biplot[,2]),
                      ss$constr.chi,
                      ss$unconst.chi,
                      ss$partial.chi,
                      ss$tot.chi,
                      perCCA1=ss$concont$importance[2,1],
                      perCCA2=ss$concont$importance[2,2])
  sumdataT<-rbind(sumdataT, sumdata) #join sumdata to sumdataT to generate total dataframe
  #pull the species loadings for CCA1 and CCA2 so they can be plotted
  #join them to a larger dataframe
  spDat1<-left_join(spDat1, rownames_to_column(as.data.frame(ss$species[,1])), by="rowname") 
  spDat2<-left_join(spDat2, rownames_to_column(as.data.frame(ss$species[,2])), by="rowname")
  
  #run the null global model to determine which variables to add during ordistep
  cca.samp0<-cca(sampdf[,-c(1:4,50:64)]~1+Condition(HUC12num), sampdf, scaling=2)
  an<-anova(cca.samp) #if the global model is significant, can do forward model selection
  #run ordistep to see which variables are significant and add information
  if(an$`Pr(>F)`[1] < 0.05 ) {
    stepx<-ordistep(cca.samp0,scope=formula(cca.samp), direction="both", 
                    trace=F)
  }else{stepx<-NULL} #if global wasn't significant, just have an empty stepx
  #build a table of the results for each ordistep result; 
  #should have run included and multiple rows if multiple variables
  if(is.data.frame(stepx$anova)){
    ax<-data.frame(run=k,
                   var=rownames(stepx$anova),
                   AIC=stepx$anova[,2],
                   Fval=stepx$anova[,3],
                   orp=stepx$anova[,4])} else {
                     
                     ax<-data.frame(run=k,
                                    var='none',
                                    AIC=NA,
                                    Fval=NA,
                                    orp=NA)}
  ANO<-rbind(ANO, ax) # join run table with total table
} 
head(sumdataT) #see the CCA summary stats
# determine % inertia explained descriptive statistics (mean, sd, median)
(Keep<-sumdataT %>% as_tibble() %>%
    dplyr::select(-run, -X.Average.of.STDM..g.m.2..,-Benthic_CHLA_MG.M2,
           -Discharge.cms,-Dvar,-Discharge.cms.1,-Dvar.1,
           -X.Average.of.STDM..g.m.2...1,-Benthic_CHLA_MG.M2.1) %>%
    mutate(perconstraint=ss.constr.chi/ss.tot.chi*100,
           perunconstrain=ss.unconst.chi/ss.tot.chi*100,
           perpartial=ss.partial.chi/ss.tot.chi*100) %>%
    dplyr::select(-ss.constr.chi, -ss.unconst.chi, -ss.partial.chi) %>%
  summarise_if(is.numeric, .funs=c(mean,sd, median)))
write_csv(sumdataT, "bootstrappedCCAres.csv") #save this so I don't have to run it

sumdataT %>% 
  dplyr::select(X.Average.of.STDM..g.m.2..:Dvar.1) %>%
  summarize_all(mean)

sumdataT %>% 
  dplyr::select(X.Average.of.STDM..g.m.2..:Dvar.1) %>%
  summarize_all(sd)

head(ANO) #see the results from ordisteps
# look at which variables are most common and the average p
ANO %>% group_by(var) %>% 
  summarize(n=n(),
            meanp=mean(orp,na.rm=T)) %>%
  arrange(desc(n)) 
# look at which variables were considered most important (chosen first)
ANO %>% group_by(run) %>% slice(1) %>% 
  ungroup() %>% group_by(var) %>% summarize(countvar=n()) %>% 
  arrange(desc(countvar)) %>%
  mutate(per_firstsig=countvar/sum(countvar)*100)
write_csv(ANO, "bootstrappedmodelselect.csv")

# which taxa had the lowest and highest loadings on CCA1 respectively
cc1tlow<-spDat1 %>% column_to_rownames() %>% rowMeans(na.rm=T) %>% 
  sort() %>% head(3) %>% names()
cc1thigh<-spDat1 %>% column_to_rownames() %>% rowMeans(na.rm=T) %>% 
  sort() %>% tail(3) %>% names()

# which taxa had the lowest and highest loadings on CCA2 respectively
cc2tlow<-spDat2 %>% column_to_rownames() %>% rowMeans(na.rm=T) %>% 
  sort() %>% head(3) %>% names()
cc2thigh<-spDat2 %>% column_to_rownames() %>% rowMeans(na.rm=T) %>% 
  sort() %>% tail(3) %>% names()

### graphing -----------------------------

# environmental loadings on each axis
# dataframe for graphing
BootEnvGraphD<-sumdataT %>% 
  dplyr::select(run, X.Average.of.STDM..g.m.2..,Benthic_CHLA_MG.M2,Discharge.cms,Dvar,
         X.Average.of.STDM..g.m.2...1, Benthic_CHLA_MG.M2.1, Discharge.cms.1, Dvar.1) %>%
  gather(Variable, value, -run) %>%
  mutate(axis=case_when(Variable %in% c("X.Average.of.STDM..g.m.2..",
                                        "Benthic_CHLA_MG.M2",
                                        "Discharge.cms","Dvar")~"CCA1",
                        Variable %in% c("X.Average.of.STDM..g.m.2...1",
                                        "Benthic_CHLA_MG.M2.1",
                                        "Discharge.cms.1","Dvar.1")~"CCA2"),
         Label=case_when(substr(Variable,1,4)=="X.Av"~"Mussel BM",
                         substr(Variable,1,4)=="Bent"~"Chl.a",
                         substr(Variable,1,4)=="Disc"~"Discharge",
                         substr(Variable,1,4)=="Dvar"~"Substrate")) %>%
  group_by(run, Label) %>% 
  dplyr::select(-Variable) %>% spread(axis, value)
  
# getting the means to add to the graph  
BootEnvGraphSum<- BootEnvGraphD %>% group_by(Label) %>%
  summarize(meanCC1=mean(CCA1),
            sdCC1=sd(CCA1),
            meanCC2=mean(CCA2),
            sdCC2=sd(CCA2))
			
# plot for environmental variable
ggplot() +
  scale_x_continuous(name="CCA 1")+
  scale_y_continuous(name="CCA 2") +
  #high transparency for each bootstrap run
  geom_segment(data=BootEnvGraphD, 
               aes(x = 0, y = 0, xend = CCA1, yend = CCA2,
                   color=Label), alpha=.05,
               arrow = arrow(length = unit(0.0, "cm")))+
  # adding the means on top to see how they compare
  geom_segment(data = BootEnvGraphSum,
             aes(x = 0, y = 0, xend = meanCC1, yend = meanCC2),
             color='black', lwd=1.05,
               arrow = arrow(length = unit(0.0, "cm")))+
  scale_color_viridis_d(guide=F) + 
  facet_grid(~Label)+ #avoid overlapping environmental variables by facetting
  theme_bw()
ggsave("Figures/S3bootenv.tiff", width=6.5, height=2.6)

# graphing the species data
colnames(spDat1) <- c("Taxa", seq(1:1000)) #renaming these nightmare of a dataframe
#building a long data for CCA1
BootSpGraph1<-spDat1 %>% group_by(Taxa) %>% gather(run, value,-Taxa) %>% 
  mutate(axis="CCA1") 
colnames(spDat2) <- c("Taxa", seq(1:1000))#renaming these nightmare of a dataframe
#building a long data for CCA2
BootSpGraph2<-spDat2 %>% group_by(Taxa) %>% gather(run, value,-Taxa) %>% 
  mutate(axis="CCA2")
 #joining both CCA1 and CCA2 to create the long dataframe
BootSpGraph<-rbind(BootSpGraph1, BootSpGraph2) %>% group_by(run, Taxa) %>%
  spread(axis, value) %>%
  filter(Taxa %in% c(cc1thigh,cc1tlow,cc2thigh, cc2tlow)) %>% 
  left_join(FTaxaTable) #adding this for good graph labels
# finding the mean value to plot on top
BootSpGraphSum<-BootSpGraph %>% group_by(Taxa, Family, Tlabel, Order) %>%
  summarize(CCA1=mean(CCA1, na.rm=T),CCA2=mean(CCA2, na.rm=T))

# plot of species loadings on CCA1 and CCA2
ggplot() +
  #low transparency for bootstrapped values
  geom_point(data=BootSpGraph, aes(x=CCA1, y=CCA2, color=Order), 
             alpha=0.05)+
  #adding the mean values on top
  geom_point(data=BootSpGraphSum, aes(x=CCA1, y=CCA2), 
             fill='darkgrey', color='black',shape=21, size=2)+
  scale_x_continuous(name="CCA 1")+
  scale_y_continuous(name="CCA 2") +
  #adding the environmental variables to look at correlationsh
  geom_segment(data = BootEnvGraphSum,
               aes(x = 0, y = 0, xend = meanCC1, yend = meanCC2),
               color='black', lwd=.5,
               arrow = arrow(length = unit(0.0, "cm")))+
  scale_color_viridis_d()+
  facet_wrap(~Tlabel)+ theme_bw()+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  theme(legend.direction = 'horizontal',legend.position = 'bottom')
ggsave("Figures/S4bootspp.tiff", width=6, height=6.3)
