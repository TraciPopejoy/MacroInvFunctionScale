install.packages("boot",dep=TRUE)
library("boot")

df<-samp2016count %>% spread(Taxa, Count) %>% 
  filter(SamplingSeason == "Fall2016")%>%
  dplyr::select(-SamplingSeason, 
                -Hem.Gerridae,
                -Dip.Chaoboridae,#only in 1
                -Dip.Orthorrhaphous, -Dip.Thaumaleidae, -Lep.Pyralidae, #only in 1
                #-Neu.Sisyridae,-Col.Hydrophilidae,-Dip.Tupulidae, #only in 2
                #-Tri.Odontoceridae, -`FLAT WORMS`, #only in 2 samples
                -Biv.Unionidae) %>% #don't care about these
  replace(is.na(.),0) %>%
  mutate(Reach=substr(Sample, 1,5),
         Site=substr(Sample, 1,2),
         Treat=substr(Sample, 4,5))%>%
  select(Reach, Site, Treat, Sample, everything()) %>% ungroup() %>%
  left_join(field.env)

un <- function(data){
  cca.samp<-cca(data[,-c(1:4,55:69)]~`Average of STDM (g.m-2)`+
                  Benthic_CHLA_MG.M2+Discharge.cms+Dvar+
                  Condition(HUC12num), data, scaling=2)
  ss<-summary(cca.samp)
  
  #cca.samp0<-cca(data[,-c(1:4,55:69)]~1+Condition(HUC12num), data, scaling=2)
  an<-anova(cca.samp) #if the global model is significant, can do forward model selection
  #if(an$`Pr(>F)`[1] < 0.05 ) {
  #  stepx<-ordistep(cca.samp0,scope=formula(cca.samp), direction="forward", 
  #                  trace=F)
  #  stepx$anova$padj <- p.adjust(stepx$anova$`Pr(>F)`, method = 'holm', n = 4)
  #}else{stepx<-NULL}
  #if(is.data.frame(stepx$anova)){
  #  ax<-data.frame(var=rownames(stepx$anova),
  #                 AIC=stepx$anova[,2],
  #                 Fval=stepx$anova[,3],
  #                 orp=stepx$anova[,4],
  #                 padj=stepx$anova[,5])}
  #else {
  #  ax<-data.frame(var='none',
  #                 AIC=NA,
  #                 Fval=NA,
  #                 orp=NA,
  #                 padj=NA)}
  sumdata<-c(t(ss$biplot[,1]),#t1:4 Mussel, Ben, Dis, Dvar
             t(ss$biplot[,2]),#t5:8 Mussle, Ben, Dis, Dvar
             ss$constr.chi/ss$tot.chi*100,#t9
             ss$unconst.chi/ss$tot.chi*100, #t10
             ss$partial.chi/ss$tot.chi*100, #t11
             ss$tot.chi, #t12
             an$`Pr(>F)`[1], #t13
             ss$concont$importance[2,], #t14:t17
             ss$species[rownames(ss$species)=="Tri.Odontoceridae",1][1],#t18
             ss$species[rownames(ss$species)=="Tri.Odontoceridae",2][1], 
             ss$species[rownames(ss$species)=="Dip.Similidae",1][1], #t20
             ss$species[rownames(ss$species)=="Dip.Similidae",2][1],
             ss$species[rownames(ss$species)=="Od.Damselflies",1][1], #t22
             ss$species[rownames(ss$species)=="Od.Damselflies",2][1],
             ss$species[rownames(ss$species)=="Tri.Helicophyche",1][1],#t24
             ss$species[rownames(ss$species)=="Tri.Helicophyche",2][1],
             ss$partial.chi, #t26
             ss$constr.chi #t27
             ) 
  #giveme<-list(ccaResults=sumdata, 
  #             anovaResults=ax)
  return(sumdata)
  }

#use the function stratified to get the correct sampling
library(devtools)
source_gist("https://gist.github.com/mrdwab/6424112")

rgen <- function(df,stratified){
  stratified(df, group = "Reach", size = 1)  
 }
  
test.boot <- boot(df, un, 1000, sim = "parametric", ran.gen = rgen)
print(test.boot)
boot.ci(test.boot, type="norm",index=11, conf=.5) # % space 
boot.ci(test.boot, type="norm",index=9) # % constrained condition 
boot.ci(test.boot, type="norm",index=10) # % unconstrained condition 
boot.ci(test.boot, type="norm",index=14) # % CCA1
boot.ci(test.boot, type="basic",index=15) # % CCA2
boot.ci(test.boot, type="basic",index=16) # % CCA3 
boot.ci(test.boot, type="basic",index=17) # % CCA4
boot.graph<-data.frame(
  C1MusMean=boot.ci(test.boot, type="basic",index=1)$t0,
  C1Mus05=boot.ci(test.boot, type="basic",index=1)$basic[4],
  C1Mus95=boot.ci(test.boot, type="basic",index=1)$basic[5],
  C1ChlMean=boot.ci(test.boot, type="basic",index=2)$t0,
  C1Chl05=boot.ci(test.boot, type="basic",index=2)$basic[4],
  C1Chl95=boot.ci(test.boot, type="basic",index=2)$basic[5],  
  C1DisMean=boot.ci(test.boot, type="basic",index=3)$t0,
  C1Dis05=boot.ci(test.boot, type="basic",index=3)$basic[4],
  C1Dis95=boot.ci(test.boot, type="basic",index=3)$basic[5],  
  C1SubMean=boot.ci(test.boot, type="basic",index=4)$t0,
  C1Sub05=boot.ci(test.boot, type="basic",index=4)$basic[4],
  C1Sub95=boot.ci(test.boot, type="basic",index=4)$basic[5],  
  C2MusMean=boot.ci(test.boot, type="basic",index=5)$t0,
  C2Mus05=boot.ci(test.boot, type="basic",index=5)$basic[4],
  C2Mus95=boot.ci(test.boot, type="basic",index=5)$basic[5],
  C2ChlMean=boot.ci(test.boot, type="basic",index=6)$t0,
  C2Chl05=boot.ci(test.boot, type="basic",index=6)$basic[4],
  C2Chl95=boot.ci(test.boot, type="basic",index=6)$basic[5],
  C2DisMean=boot.ci(test.boot, type="basic",index=7)$t0,
  C2Dis05=boot.ci(test.boot, type="basic",index=7)$basic[4],
  C2Dis95=boot.ci(test.boot, type="basic",index=7)$basic[5],
  C2SubMean=boot.ci(test.boot, type="basic",index=8)$t0,
  C2Sub05=boot.ci(test.boot, type="basic",index=8)$basic[4],
  C2Sub95=boot.ci(test.boot, type="basic",index=8)$basic[5])


boot.graph1<-as.data.frame(t(boot.graph)) %>%
  mutate(Variable=substring(rownames(.), 3,5),
         CCaxis=substring(rownames(.), 1,2),
         stattype=substring(rownames(.), 6,7)) %>%
  mutate(CCaxis.stat=paste(CCaxis,stattype, sep=".")) %>%
  select(-CCaxis, -stattype) %>%
  spread(CCaxis.stat, V1) %>%
  mutate(Label=c("Benthic Chl a", "Discharge",
                  "Mussel Dry Mass", "Substrate Variability"))

#from CCA.R

FieldCCA

ggplot() +
  #geom_point(data=tidy.ccaFSp, aes(x=CCA1, y=CCA2), color="grey")+
  #geom_point(data=tidy.ccaFspespe, aes(x=CCA1,y=CCA2)) +
  scale_x_continuous(name="CCA 1")+
  scale_y_continuous(name="CCA 2") +
  #geom_segment(data = tidy.ccaFcon[tidy.ccaFcon$Score=="biplot",],
  #             aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
  #             arrow = arrow(length = unit(0.0, "cm")), 
  #             color="red")+
  geom_segment(data = boot.graph1,
               aes(x = 0, y = 0, xend = C1.Me, yend = C2.Me),
               arrow = arrow(length = unit(0.0, "cm")))+
  geom_errorbarh(data=boot.graph1, aes(xmin=C1.05, xmax=C1.95, y=C2.Me, color=Variable), height=0)+
  geom_linerange(data=boot.graph1, aes(ymin=C2.05, ymax=C2.95, x=C1.Me, color=Variable))+
  geom_text_repel(data=boot.graph1, aes(x=C1.Me, y=C2.Me, label=Variable))+
  #geom_text_repel(data=frep[frep$Score=="biplot",], aes(x=CCA1, y=CCA2, label=Tlabel), size=2)+
  scale_color_viridis_d() +theme_bw()

bootres<-as.data.frame(test.boot$t)
names(bootres)<-c("C1Mus","C1Chl","C1Dis","C1Sub",
                  "C2Mus","C2Chl","C2Dis","C2Sub")
bootres1<-bootres[,1:8] %>% mutate(nob=1:n()) %>% group_by(nob) %>%
  gather(Variable, value, -nob) %>%
  mutate(axis=substr(Variable, 1,2),
         var=substr(Variable, 3,5)) %>% 
  spread(axis, value) %>% group_by(nob, var) %>%
  summarise_if(is.numeric,list(~.[which(!is.na(.))])) %>%
  mutate(Lab=case_when(var=="Chl"~"Benthic Chl a", 
                       var=="Dis"~"Discharge",
                       var=="Mus"~"Mussel Dry Mass",
                       var=="Sub"~"Substrate Variability"))
bootres1
ggplot() +
  scale_x_continuous(name="CCA 1")+
  scale_y_continuous(name="CCA 2") +
  geom_segment(data = bootres1, 
               aes(x = 0, y = 0, xend = C1, yend = C2,
                   color=var), alpha=.07,
               arrow = arrow(length = unit(0.0, "cm")))+
   #geom_segment(data = boot.graph1,
   #           aes(x = 0, y = 0, xend = C1.Me, yend = C2.Me,
   #                color=Variable),linetype='dashed', lwd=1.1,
   #             arrow = arrow(length = unit(0.0, "cm")))+
  scale_color_viridis_d(guide=F) + facet_wrap(~Lab)+ theme_bw()
ggsave("Figures/SupBootCCA.tiff")
