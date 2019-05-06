source("InvDataframes.R")

F.ComDens16<-F.ComDens %>% filter(SamplingSeason=="Fall2016")
E.ComDens12<-E.ComDens %>% filter(Week=="w12")
# Gamma & Alpha Diversity --------

# Field Data
#head(F.ComDens)
sum(colSums(F.ComDens16[,-c(1:7,65)])>1) #gamma diversity?

F.Talpha<-data.frame(F.ComDens[,c(1:7)],
           richness=specnumber(F.ComDens[,-c(1:7, 65)]),
           SimpsonsI=diversity(F.ComDens[,-c(1:7,65)], "simpson")) %>%
  group_by(Reach, SamplingSeason, Treatment, Season, Mussel.g.m2) %>%
  summarize(meanID.npm2=mean(InvDensity.npm2R),
            meanR=mean(richness),
            meanSimp=mean(SimpsonsI)) %>%
  mutate(TreatFAC=factor(Treatment, levels=c("NM","MR")))
mean(F.Talpha$meanR);sd(F.Talpha$meanR)
ggplot(F.Talpha, aes(x=Treatment, y=meanSimp))+
  geom_boxplot(fill="grey")+facet_wrap(~Season)+
  ylab("Simpson Diversity Index")
ggplot(F.Talpha, aes(x=Treatment, y=meanR))+
  geom_boxplot(fill="grey")+facet_wrap(~Season)+
  ylab("Richness")
FIn<-ggplot(F.Talpha, aes(x=TreatFAC, y=meanID.npm2))+
  geom_boxplot(fill="grey")+
  scale_y_log10(name=expression("Invertebrate # m "^-2))+
  scale_x_discrete(name="Treatment",
                   labels=c("Control\nReach","Mussel\nReach"))

#Enclosure Data
#head(E.ComDens)
sum(colSums(E.ComDens12[,-c(1:8)])>1) #gamma diversity?

E.Talpha<-data.frame(E.ComDens12[,c(1:8)],
                     InvDensity.npm2=rowSums(E.ComDens12[,-c(1:8)]),
                     richness=specnumber(E.ComDens12[,-c(1:8)]),
                     SimpsonsI=diversity(E.ComDens12[,-c(1:8)], "simpson"),
                     TreatFAC=factor(E.ComDens12$TreatA, 
                                    levels=c("CTRL","ACTL","ACTS","AMBL","AMBS"),
                                    labels=c("Control",
                                             'Actinonaias\nLive',
                                             'Actinonaias\nSham',
                                             'Amblema\nLive',
                                             'Amblema\nSham')))
ggplot(E.Talpha, aes(x=TreatA, y=SimpsonsI))+
  geom_boxplot(fill="grey")+facet_wrap(~Week)+theme_krementz()
ggplot(E.Talpha, aes(x=TreatA, y=richness))+
  geom_boxplot(fill="grey")+facet_wrap(~Week)+theme_krementz()
EIn<-ggplot(E.Talpha, aes(x=TreatFAC, y=InvDensity.npm2))+
  geom_boxplot(fill="grey")+
  scale_y_log10(name=expression("Invertebrate # m "^-2))+
  scale_x_discrete(name="Treatment")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))


#Shell Data
#head(E.ComDens)
sum(colSums(S.ComDens[,-c(1:7)])>0) #gamma diversity?

S.Talpha<-data.frame(S.ComDens[,c(1:7)],
                     InvDensity.npcm2=rowSums(S.ComDens[,-c(1:7)]),
                     richness=specnumber(S.ComDens[,-c(1:7)]),
                     SimpsonsI=diversity(S.ComDens[,-c(1:7)], "simpson")) %>%
  left_join(ShellChl)
ggplot(S.Talpha, aes(x=TreatA, y=SimpsonsI))+
  geom_boxplot(fill="grey")+facet_wrap(~SheType, scales="free_x")+theme_krementz()
ggplot(S.Talpha, aes(x=TreatA, y=richness))+
  geom_boxplot(fill="grey")+facet_wrap(~SheType, scales="free_x")+theme_krementz()
SIn<-ggplot(S.Talpha, aes(x=SheType, y=InvDensity.npcm2*10000))+
  geom_boxplot(fill="grey")+
  scale_y_log10(name=expression("Invertebrate # m "^-2),)+
  scale_x_discrete(name="Treatment", labels=c('Actinonaias\nLive','Actinonaias\nSham',
                                              'Amblema\nLive','Amblema\nSham'))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
library(cowplot)
plot_grid(FIn,EIn,SIn, nrow=1,labels="AUTO")
ggsave("./Figures/AbundTreat.tiff",width=10, height=4)

FIbm<-ggplot(F.Talpha, aes(x=Mussel.g.m2, y=meanID.npm2))+
  geom_point()+
  geom_smooth(method="lm",color="black", linetype="dashed")+
  scale_y_log10(name=expression("Invertebrates # m "^-2))+
  scale_x_continuous(trans="log1p", breaks=c(0,5,10,50,200),
                    name=expression("Mussel biomass g m "^-2))
EIbm<-ggplot(E.Talpha, aes(x=MusselBiomass.g.m2, y=InvDensity.npm2))+
  geom_point()+
  geom_smooth(method="lm",color="black", linetype="dashed")+
  scale_y_log10(name=expression("Invertebrates # m "^-2))+
  scale_x_continuous(trans="log1p", breaks=c(0,5,10,50,150,300),
                     name=expression("Mussel biomass g m "^-2))
SIbm<-ggplot(S.Talpha, aes(x=TShellSurArea.cm2, y=InvDensity.npcm2*10000))+
  geom_point()+
  geom_smooth(method="lm",color="black", linetype="dashed")+
  scale_y_log10(name=expression("Invertebrates # m "^-2))+
  scale_x_continuous(trans="log1p", breaks=c(500,1000,1500,2500),
                     name=expression("Shell Area sampled cm "^2))
plot_grid(FIbm, EIbm, SIbm, nrow = 1, labels="AUTO")
ggsave("./Figures/MussAbund.tiff",width=10, height=4)

#Table1
F.TalphaS<-F.Talpha %>% group_by(Treatment,SamplingSeason) %>%
  summarize(meanRt=mean(meanR),
            sdR=sd(meanR),
            meanIDt.npm2=mean(meanID.npm2),
            sdID=sd(meanID.npm2, na.rm=T))
E.TalphaS<-E.Talpha %>% group_by(TreatA,Week) %>%
  summarize(meanRt=mean(richness),
            sdR=sd(richness),
            meanIDt.npm2=mean(InvDensity.npm2),
            sdID=sd(InvDensity.npm2, na.rm=T),
            meanMBM=mean(MusselBiomass.g.m2))
S.TalphaS<-S.Talpha %>% group_by(TreatA,ShellSpecies) %>%
  summarize(meanRt=mean(richness),
            sdR=sd(richness),
            meanIDt.npm2=mean(InvDensity.npcm2)*10000,
            sdID=sd(InvDensity.npcm2, na.rm=T)*10000)
tab1<-rbind(F.TalphaS, E.TalphaS, S.TalphaS)
write.csv(tab1, "Table1.InvertBasics.csv")

# Beta Diversity ---------

#Field Data
F.Tbeta<-betadiver(F.ComDens[,-c(1:7,65)], method="sor")
(F.Beta.mod <- betadisper(F.Tbeta, paste(F.ComDens$Treatment,F.ComDens$SamplingSeason)))
## Perform test
anova(F.Beta.mod)
## Tukey's Honest Significant Differences
plot(TukeyHSD(F.Beta.mod))
boxplot(F.Beta.mod)

#Enclosure Data
E.Tbeta<-betadiver(E.ComDens[,-c(1:8)], method="sor")
(E.Beta.mod <- betadisper(E.Tbeta, E.ComDens$TreatA))
## Perform test
anova(E.Beta.mod)
## Tukey's Honest Significant Differences
plot(TukeyHSD(E.Beta.mod)) 
boxplot(E.Beta.mod)

#Shell Data
S.Tbeta<-betadiver(S.ComDens[,-c(1:7)], method="sor")
(S.Beta.mod <- betadisper(S.Tbeta, S.ComDens$TreatA))
## Perform test
anova(S.Beta.mod)
## Tukey's Honest Significant Differences
plot(TukeyHSD(S.Beta.mod)) 
boxplot(S.Beta.mod)

# Most Abundant Families -----------
library(devtools)
source_gist("524eade46135f6348140")

# Field Data
F.ComDens %>% dplyr::select(-InvDensity.npm2,-Mussel.g.m2,-Other.other) %>%
  gather(Taxa, Abundance, -FSamID,-SamplingSeason,-Season,-Reach,-Treatment) %>% 
  group_by(Taxa) %>%
  summarize(mean.A=mean(Abundance), sd.A=sd(Abundance),sum.A=sum(Abundance)) %>%
  arrange(desc(sum.A))
#results: ChironomidaeL, Hydropsychidae, Caenidae, ElmL, Oligo
# Enclosure Data
E.ComDens %>% dplyr::select(-Enc2,-Week,-TreatA,-Type,-Spp,-MusselBiomass.g.m2) %>%
  gather(Taxa, Abundance,-TEid,-Enc) %>% 
  group_by(Taxa) %>%
  summarize(mean.A=mean(Abundance), sd.A=sd(Abundance),sum.A=sum(Abundance)) %>%
  arrange(desc(sum.A))
#results: Heptageniidae, ChrionomidaeL, Polycentrop, Damselflies, Caenidae, ElmL
# Shell Data
S.ComDens %>% dplyr::select(-SamID,-Enc2,-ShellSpecies,-TreatA,-Type,-Spp, -SheType) %>%
  gather(Taxa, Abundance) %>% 
  group_by(Taxa) %>%
  summarize(mean.A=mean(Abundance), sd.A=sd(Abundance),sum.A=sum(Abundance)) %>%
  arrange(desc(sum.A))
#results: ChironomidaeL, Polycentropidae, Heptageniidae, ElmL

### So I am interested in:
###  Chironomidae, Heptageniidae, Caenidae, Polycentropidae, Elmidae
### function to plot regression on figure ####
# source: https://gist.github.com/kdauria/524eade46135f6348140
stat_smooth_func_with_pval <- function(mapping = NULL, data = NULL,
                                       geom = "smooth", position = "identity",
                                       ...,
                                       method = "auto",
                                       formula = y ~ x,
                                       se = TRUE,
                                       n = 80,
                                       span = 0.75,
                                       fullrange = FALSE,
                                       level = 0.95,
                                       method.args = list(),
                                       na.rm = FALSE,
                                       show.legend = NA,
                                       inherit.aes = TRUE,
                                       xpos = NULL,
                                       ypos = NULL,
                                       xpos2 = NULL,
                                       ypos2 = NULL) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatSmoothFunc,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      method = method,
      formula = formula,
      se = se,
      n = n,
      fullrange = fullrange,
      level = level,
      na.rm = na.rm,
      method.args = method.args,
      span = span,
      xpos = xpos,
      ypos = ypos,
      xpos2 = xpos2,
      ypos2 = ypos2,
      ...
    )
  )
}

StatSmoothFunc <- ggproto("StatSmooth", Stat,
                          
                          setup_params = function(data, params) {
                            # Figure out what type of smoothing to do: loess for small datasets,
                            # gam with a cubic regression basis for large data
                            # This is based on the size of the _largest_ group.
                            if (identical(params$method, "auto")) {
                              max_group <- max(table(data$group))
                              
                              if (max_group < 1000) {
                                params$method <- "loess"
                              } else {
                                params$method <- "gam"
                                params$formula <- y ~ s(x, bs = "cs")
                              }
                            }
                            if (identical(params$method, "gam")) {
                              params$method <- mgcv::gam
                            }
                            
                            params
                          },
                          
                          compute_group = function(data, scales, method = "auto", formula = y~x,
                                                   se = TRUE, n = 80, span = 0.75, fullrange = FALSE,
                                                   xseq = NULL, level = 0.95, method.args = list(),
                                                   na.rm = FALSE, xpos=NULL, ypos=NULL, 
                                                   xpos2=NULL, ypos2=NULL) {
                            if (length(unique(data$x)) < 2) {
                              # Not enough data to perform fit
                              return(data.frame())
                            }
                            
                            if (is.null(data$weight)) data$weight <- 1
                            
                            if (is.null(xseq)) {
                              if (is.integer(data$x)) {
                                if (fullrange) {
                                  xseq <- scales$x$dimension()
                                } else {
                                  xseq <- sort(unique(data$x))
                                }
                              } else {
                                if (fullrange) {
                                  range <- scales$x$dimension()
                                } else {
                                  range <- range(data$x, na.rm = TRUE)
                                }
                                xseq <- seq(range[1], range[2], length.out = n)
                              }
                            }
                            # Special case span because it's the most commonly used model argument
                            if (identical(method, "loess")) {
                              method.args$span <- span
                            }
                            
                            if (is.character(method)) method <- match.fun(method)
                            
                            base.args <- list(quote(formula), data = quote(data), weights = quote(weight))
                            model <- do.call(method, c(base.args, method.args))
                            
                            m = model
                            eq1 <- substitute(italic(y) == a + b %.% italic(x), 
                                              list(a = format(coef(m)[1], digits = 3), 
                                                   b = format(coef(m)[2], digits = 3)))
                            
                            eq2 <- substitute(italic(r)^2~"="~r2*","~~italic(p)~"="~pval, 
                                              list(r2 = format(summary(m)$r.squared, digits = 3),
                                                   pval = format(summary(m)$coef[2,4], digits = 3)))
                            
                            func_string1 = as.character(as.expression(eq1))
                            func_string2 = as.character(as.expression(eq2))
                            
                            if(is.null(xpos)) xpos = min(data$x)*0.9
                            if(is.null(ypos)) ypos = max(data$y)*0.9
                            if(is.null(xpos2)) xpos2 = xpos
                            if(is.null(ypos2)) ypos2 = max(data$y)*0.6
                            
                            data.frame(x = rbind(xpos, xpos2), 
                                       y = rbind(ypos, ypos2), 
                                       label = rbind(func_string1, func_string2))
                            
                          },
                         required_aes = c("x", "y")
)
### PLOTTING families in samples ####
### because they are abundant in most of my samples
library(ggsci)
#ChironomidaeL
F.ComDens2<-F.ComDens %>% mutate(logMB.g.m2=log(Mussel.g.m2))
F.ComDens2$logMB.g.m2 <- replace(F.ComDens2$logMB.g.m2,F.ComDens2$logMB.g.m2==-Inf,0)
ggplot(F.ComDens2, aes(x=logMB.g.m2, y=Dip.ChironomidaeL))+
  geom_smooth(method="lm", color="black", se=F)+
  geom_point(size=3, aes(color=SamplingSeason))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  scale_color_futurama(name="Season")+theme_krementz()+
  theme(legend.position="bottom", legend.direction="horizontal")+
  labs(y="Chironomid n/m2", x="log2 Mussel Biomass g/m2")
E.ComDens2<-E.ComDens %>% mutate(logMB.g.m2=log(MusselBiomass.g.m2))
E.ComDens2$logMB.g.m2 <- replace(E.ComDens2$logMB.g.m2,E.ComDens2$logMB.g.m2==-Inf,0)
ggplot(E.ComDens2, aes(x=logMB.g.m2, y=Dip.ChironomidaeL))+
  geom_smooth(method="lm", color="black")+
  geom_point(size=3, aes(color=TreatA))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  scale_color_futurama()+theme_krementz()+
  theme(legend.position="bottom", legend.direction="horizontal")+
  labs(y="Chironomid n/m2", x="log2 Mussel Biomass g/m2")
ggplot(S.ComDens, aes(x=SheType,y=Dip.ChironomidaeL))+
  geom_boxplot()+
  geom_point(size=3, aes(color=TreatA))+
  scale_color_futurama()+theme_krementz()+
  theme(legend.position="bottom", legend.direction="horizontal")
S.ComDens2<-S.ComDens %>% left_join(MDat)
ggplot(S.ComDens2, aes(x=TShellSurArea.cm2,y=Dip.ChironomidaeL))+
  geom_smooth(method="lm", color="black", se=F)+
  geom_point(size=3, aes(color=TreatA, shape=SheType))+
  #stat_smooth_func(geom="text",method="lm",hjust=-.7,parse=TRUE) +
  scale_shape_discrete(name="Shell Type")+
  scale_color_futurama(name="Treat")+theme_krementz()+
  labs(y="Chironomid n/cm2", x="Shell Area Sampled cm2")

#Heptageniidae
ggplot(F.ComDens2, aes(x=logMB.g.m2, y=Eph.Heptageniidae))+
  geom_smooth(method="lm", color="black", se=F)+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point(size=3, aes(color=SamplingSeason))+
  scale_color_futurama()+theme_krementz()+
  theme(legend.position="bottom", legend.direction="horizontal")+
  labs(y="Heptageniidae n/m2",x="log2 Mussel Biomass g/m2")
ggplot(E.ComDens2, aes(x=logMB.g.m2, y=Eph.Heptageniidae))+
  geom_smooth(method="lm", color="black")+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point(size=3, aes(color=TreatA))+
  scale_color_futurama()+theme_krementz()+
  labs(y="Heptageniidae n/m2",x="log2 Mussel Biomass g/m2")
ggplot(S.ComDens, aes(x=SheType,y=Eph.Heptageniidae))+
  geom_boxplot()+
  geom_point(size=2, aes(color=TreatA))+
  scale_color_futurama()
ggplot(S.ComDens2, aes(x=TShellSurArea.cm2,y=Eph.Heptageniidae))+
  geom_smooth(method="lm", color="black", se=F)+
  geom_point(size=3, aes(color=TreatA, shape=SheType))+
  #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  scale_shape_discrete(name="Shell Type")+
  scale_color_futurama(name="Treat")+theme_krementz()+
  labs(y="Heptageniidae n/cm2",x="Shell Area Sampled cm2")

#Caenidae
ggplot(F.ComDens2, aes(x=logMB.g.m2, y=Eph.Trcorythidae))+
  geom_smooth(method="lm", color="black", se=F)+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point(size=3, aes(color=SamplingSeason))+
  scale_color_futurama(guide=F)+theme_krementz()+
  labs(y="Caenidae n/m2",x="log2 Mussel Biomass g/m2")
ggplot(E.ComDens2, aes(x=logMB.g.m2, y=Eph.Trcorythidae))+
  geom_smooth(method="lm", color="black", se=F)+
  geom_point(size=3, aes(color=TreatA))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  scale_color_futurama(guide=F)+theme_krementz()+
  labs(y="Caenidae n/m2",x="log2 Mussel Biomass g/m2")
ggplot(S.ComDens, aes(x=SheType,y=Eph.Trcorythidae))+
  geom_boxplot()+
  geom_point(size=2, aes(color=TreatA))+
  scale_color_futurama()
ggplot(S.ComDens2, aes(x=TShellSurArea.cm2,y=Eph.Trcorythidae))+
  geom_smooth(method="lm", color="black", se=F)+
  geom_point(size=3, aes(color=TreatA, shape=SheType))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  scale_shape_discrete(name="Shell Type", guide=F)+
  scale_color_futurama(name="Treat", guide=F)+theme_krementz()+
  labs(y="Caenidae n/cm2",x="Shell Area Sampled cm2")

#Polycentropidae
ggplot(F.ComDens2, aes(x=logMB.g.m2, y=Tri.Polycentropidae))+
  geom_smooth(method="lm", color="black", se=F)+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point(size=3, aes(color=SamplingSeason))+
  scale_color_futurama(guide=F)+theme_krementz()+
  labs(y="Polycentropidae n/m2",x="log2 Mussel Biomass g/m2")
ggplot(E.ComDens2, aes(x=logMB.g.m2, y=Tri.Polycentropidae))+
  geom_smooth(method="lm", color="black", se=F)+
  geom_point(size=3, aes(color=TreatA))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  scale_color_futurama(guide=F)+theme_krementz()+
  labs(y="Polycentropidae n/m2",x="log2 Mussel Biomass g/m2")
ggplot(S.ComDens, aes(x=SheType,y=Tri.Polycentropidae))+
  geom_boxplot()+
  geom_point(size=2, aes(color=TreatA))+
  scale_color_futurama()
ggplot(S.ComDens2, aes(x=TShellSurArea.cm2,y=Tri.Polycentropidae))+
  geom_smooth(method="lm", color="black", se=F)+
  geom_point(size=3, aes(color=TreatA, shape=SheType))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  scale_shape_discrete(name="Shell Type", guide=F)+
  scale_color_futurama(name="Treat", guide=F)+theme_krementz()+
  labs(y="Polycentropidae n/cm2",x="Shell Area Sampled cm2")

#Elmidae
ggplot(F.ComDens2, aes(x=logMB.g.m2, y=Col.ElmL))+
  geom_smooth(method="lm", color="black", se=F)+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point(size=3, aes(color=SamplingSeason))+
  scale_color_futurama(guide=F)+theme_krementz()+
  labs(y="Elmidae n/m2",x="log2 Mussel Biomass g/m2")
ggplot(E.ComDens2, aes(x=logMB.g.m2, y=Col.ElmL))+
  geom_smooth(method="lm", color="black", se=F)+
  geom_point(size=3, aes(color=TreatA))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  scale_color_futurama(guide=F)+theme_krementz()+
  labs(y="Elmidae n/m2",x="log2 Mussel Biomass g/m2")
ggplot(S.ComDens, aes(x=SheType,y=Col.ElmL))+
  geom_boxplot()+
  geom_point(size=2, aes(color=TreatA))+
  scale_color_futurama()
ggplot(S.ComDens2, aes(x=TShellSurArea.cm2,y=Col.ElmL))+
  geom_smooth(method="lm", color="black", se=F)+
  geom_point(size=3, aes(color=TreatA, shape=SheType))+
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  scale_shape_discrete(name="Shell Type", guide=F)+
  scale_color_futurama(name="Treat", guide=F)+theme_krementz()+
  labs(y="Elmidae n/cm2",x="Shell Area Sampled cm2")

# Fuzzy ordination --------------------
source("EnvDataframes.R")

#CCA
CCAField.data<-F.ComDens %>% left_join(FieldDischarge) %>% 
  left_join(FieldChlA) %>% 
  mutate(Site=substr(Reach, 1,2))
CCAField.data$Mussel.g.m2<-replace_na(CCAField.data$Mussel.g.m2,0)
CCAField.data[CCAField.data$FSamID=="GL-MR.Fall2015","Discharge.cms"]<-0.0294 #usgs
CCAField.data[CCAField.data$FSamID=="GL-NM.Fall2015","Discharge.cms"]<-0.0274 #usgs
CCAField.data[CCAField.data$FSamID=="KT-MR.Fall2015","Discharge.cms"]<-0.109 
CCAField.data[CCAField.data$FSamID=="L3-MR.Summer2015","Discharge.cms"]<-0.819 
CCAField.data[CCAField.data$FSamID=="L3-NM.Summer2015","Discharge.cms"]<-0.556 

names(CCAField.data[,-c(1:7,70:73)])
cca.F<-cca(CCAField.data[,-c(1:7,70:73)]~
             CCAField.data$Mussel.g.m2+CCAField.data$Discharge.cms+
             CCAField.data$Season)
cca.F
ccaF.plot<-plot(cca.F)
cca.F.sum<-summary(cca.F)
#str(cca.F.sum)
sort(cca.F.sum$biplot[,1])
sort(cca.F.sum$species[,1], decreasing=T)[1:6]
sort(cca.F.sum$species[,1], decreasing=F)[1:6]
sort(cca.F.sum$biplot[,2])
sort(cca.F.sum$species[,2], decreasing=T)[1:6]
sort(cca.F.sum$species[,2], decreasing=F)[1:6]
cca.F.sum$concont
cca.F.sum$constr.chi/cca.F.sum$tot.chi*100

# Mussel Biomass is already in all of the data sets.
library(fso)
head(FieldDischarge)
head(EncDischarge)
head(EncChlAraw)
head(ShellChl)

names(F.ComPer[,-c(1:6)])#those columns just identifiers, and others
FieldFuzzyData<-F.ComPer %>% left_join(FieldDischarge) %>%
  left_join(FieldChlA)
head(FieldFuzzyData[,-c(1:6)])
F.Fdist<-dsvdis(FieldFuzzyData[,-c(1:6,61:63)], 'bray')
FieldFuzzyData$Mussel.g.m2<-replace_na(FieldFuzzyData$Mussel.g.m2,0)
FieldFuzzyData[FieldFuzzyData$FSamID=="GL-MR.Fall2015","Discharge.cms"]<-0.0294 #usgs
FieldFuzzyData[FieldFuzzyData$FSamID=="GL-NM.Fall2015","Discharge.cms"]<-0.0274 #usgs
FieldFuzzyData[FieldFuzzyData$FSamID=="KT-MR.Fall2015","Discharge.cms"]<-0.109 
FieldFuzzyData[FieldFuzzyData$FSamID=="L3-MR.Summer2015","Discharge.cms"]<-0.819 
FieldFuzzyData[FieldFuzzyData$FSamID=="L3-NM.Summer2015","Discharge.cms"]<-0.556 

F.fsoM<-fso(FieldFuzzyData$Mussel.g.m2, F.Fdist, permute = 10000)
summary(F.fsoM)
F.fso<-fso(FieldFuzzyData$Discharge.cms, F.Fdist, permute = 10000)
summary(F.fso) #p = 0.09
F.fso.plot<-FieldFuzzyData %>% mutate(mu=F.fso$mu) %>%
  dplyr::select(-WC_CHLA_MG.L,-Benthic_CHLA_MG.M2) %>%
  gather(Taxa, Dens, -FSamID, -Reach,-SamplingSeason,
         -Treatment,-Season,-Mussel.g.m2,-Discharge.cms,-mu) %>%
  left_join(FTaxaTable)%>%
  arrange(desc(mu))
ggplot(F.fso.plot, aes(x=Discharge.cms, y=mu))+
  geom_point(aes(color=Season), size=3)+geom_smooth(method="lm")

greyC<-gray.colors(length(unique(F.fso.plot$Order)),
                   start = 0.02, end = 0.97, gamma = 2.2, alpha = NULL)
greyC2<-c(greyC[1:5],"red",greyC[7],"purple",greyC[9:16],"orange",greyC[18],"blue",
          greyC[1],"orange",greyC[21:22])
Fmu_table <- table(F.fso.plot$mu)
Fmu_levels <- names(Fmu_table)[order(Fmu_table)]
F.fso.plot$OrdMu<-factor(F.fso.plot$mu,levels=Fmu_levels)
ggplot(data = F.fso.plot, aes(x =OrdMu, y = Dens, fill=Order)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values=greyC2, name="Species")+
  coord_flip()+
  scale_x_discrete(name="Ordinated Sites", labels=unique(F.fso.plot$FSamID)) +
  theme(axis.text.x = element_text(size=9,color="black"),
        axis.title.y=element_text(size=20),
        plot.background = element_blank(),
        panel.border=element_blank(),
        panel.grid.major= element_line(colour=NA), 
        panel.grid.minor=element_line(colour=NA),
        title=element_text(size=20),
        panel.background = element_rect(fill = "white"),
        axis.line.x=element_line(colour="black"),
        axis.line.y=element_line(colour="black"))
# Dip, eph, tri, col, anne, odo, ple
F.mfso<-mfso(~FieldFuzzyData$Mussel.g.m2+FieldFuzzyData$Discharge.cms+
               FieldFuzzyData$WC_CHLA_MG.L,
           F.Fdist, permute = 10000)
summary(F.mfso)
plot(F.mfso)
F.mfso.plot<-data.frame(Reach=FieldFuzzyData[!is.na(FieldFuzzyData$WC_CHLA_MG.L)&
                                             !is.na(FieldFuzzyData$Discharge.cms),"FSamID"],
                        MB.mu=F.mfso$mu[,1],
                        Dis.mu=F.mfso$mu[,2],
                        Ch.mu=F.mfso$mu[,3]) %>% left_join(FieldFuzzyData[,c(1:6)])
ggplot(F.mfso.plot, aes(x=MB.mu, y=Ch.mu, color=Treatment, shape=Season))+
  geom_point()

### Enclosure
names(E.ComPer[,-c(1:7)])#those columns just identifiers, and others
EncFuzzyData<-E.ComPer %>% left_join(EncDischarge) %>%
  left_join(EncChlAraw)
head(EncFuzzyData[,-c(1:7,53,54)])
E.Fdist<-dsvdis(EncFuzzyData[,-c(1:7,53,54)], 'bray')
E.fso<-fso(EncFuzzyData$MusselBiomass.g.m2, E.Fdist, permute = 10000)
summary(E.fso)
E.fso.plot<-EncFuzzyData %>% mutate(mu=E.fso$mu) %>%
  dplyr::select(-Discharge.cms, -ChlAdensity) %>%
  gather(Taxa, Dens, -TEid, -Enc2,-Week,-TreatA,-Type,-Spp,-MusselBiomass.g.m2,-mu) %>%
  left_join(FTaxaTable)%>%
  arrange(desc(mu))
ggplot(E.fso.plot, aes(x=MusselBiomass.g.m2, y=mu))+
  geom_point(size=3)+geom_smooth(method="lm")

greyC<-gray.colors(length(unique(E.fso.plot$Order)),
                   start = 0.02, end = 0.97, gamma = 2.2, alpha = NULL)
greyCE<-c(greyC[1:4], "red", greyC[5], "purple",greyC[6:15],"blue", greyC[17:18],
          "orange", greyC[19:22])
Emu_table <- table(E.fso.plot$mu)
Emu_levels <- names(Emu_table)[order(Emu_table)]
E.fso.plot$OrdMu<-factor(E.fso.plot$mu,levels=Emu_levels)
ggplot(data = E.fso.plot, aes(x =OrdMu, y = Dens, fill=Order)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values=greyCE, name="Species")+
  coord_flip()+
  scale_x_discrete(name="Ordinated Sites", labels=unique(E.fso.plot$TEid)) +
  theme(axis.text.x = element_text(size=9,color="black"),
        axis.title.y=element_text(size=20),
        plot.background = element_blank(),
        panel.border=element_blank(),
        panel.grid.major= element_line(colour=NA), 
        panel.grid.minor=element_line(colour=NA),
        title=element_text(size=20),
        panel.background = element_rect(fill = "white"),
        axis.line.x=element_line(colour="black"),
        axis.line.y=element_line(colour="black"))

E.mfso<-mfso(~EncFuzzyData$MusselBiomass.g.m2+EncFuzzyData$Discharge.cms+
               EncFuzzyData$ChlAdensity, E.Fdist)
summary(E.mfso)
plot(E.mfso)
E.mfso.plot<-data.frame(Reach=EncFuzzyData[!is.na(EncFuzzyData$ChlAdensity)&
                                               !is.na(EncFuzzyData$Discharge.cms),"TEid"],
                        MB.mu=E.mfso$mu[,1],
                        Dis.mu=E.mfso$mu[,2],
                        Ch.mu=E.mfso$mu[,3]) %>% left_join(EncFuzzyData[,c(1:7)])
ggplot(E.mfso.plot, aes(x=MB.mu, y=Ch.mu, color=TreatA, shape=Type))+
  geom_point()

### Shell
names(S.ComPer[,-c(1:6)])#those columns just identifiers, and others
ShellFuzzyData<-S.ComPer %>% left_join(EncDischarge[EncDischarge$Week=="w12",]) %>%
  left_join(ShellChl)
#have a -0.5 in chlorophyll, fixing
ShellFuzzyData[ShellFuzzyData$SamID=="C09.AMB","ChlAdensity"]<-0
head(ShellFuzzyData[,-c(1:6,38:41)])
S.Fdist<-dsvdis(ShellFuzzyData[,-c(1:6,38:41)], 'bray')
S.fso<-fso(~ShellFuzzyData$ChlAdensity+ShellFuzzyData$Discharge.cms+
             ShellFuzzyData$Type+ShellFuzzyData$ShellSpecies+
             ShellFuzzyData$TShellSurArea.cm2, S.Fdist,
           permute = 10000)
S.fso<-fso(~ShellFuzzyData$Type+ShellFuzzyData$ChlAdensity, S.Fdist, permute=10000)
summary(S.fso)
S.fso.plot<-ShellFuzzyData %>% mutate(muType=S.fso$mu[,1],
                                      muChlA=S.fso$mu[,2]) %>%
  dplyr::select(-Discharge.cms) %>%
  gather(Taxa, Dens, -SamID, -Enc2,-Week,-TreatA,-Type,-Spp,-ShellSpecies,
         -TShellSurArea.cm2,-ChlAdensity,-muType, -muChlA) %>%
  left_join(FTaxaTable) %>% arrange(desc(muType))
ggplot(S.fso.plot, aes(x=ChlAdensity, y=muChlA))+
  geom_point(size=3)+geom_smooth(method="lm")

greyC<-gray.colors(length(unique(S.fso.plot$Taxa)),
                   start = 0.02, end = 0.97, gamma = 2.2, alpha = NULL)
greyCS<-c(greyC[1:3], "red", greyC[7],"purple",greyC[9:12], "blue", greyC[19:22])
Smu_table <- table(S.fso.plot$muChlA)
Smu_levels <- names(Smu_table)[order(Smu_table)]
S.fso.plot$OrdMu<-factor(S.fso.plot$muChlA,levels=Smu_levels)
shell.labels<-S.fso.plot %>% group_by(SamID) %>%
  slice(1) %>% dplyr::select(SamID, Enc2, ShellSpecies, 
                             TreatA, Type, Spp, muChlA, ChlAdensity) %>%
  arrange(desc(muChlA))
ggplot(data = S.fso.plot, aes(x =OrdMu, y = Dens, fill=Order)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values=greyCS, name="Species")+
  coord_flip()+
  scale_x_discrete(name="Ordinated Sites", 
                   labels=paste(shell.labels$TreatA, round(shell.labels$ChlAdensity,2)))+
                   #labels=paste(shell.labels$SamID, round(shell.labels$muType,2))) +
  theme(axis.text.x = element_text(size=9,color="black"),
        axis.title.y=element_text(size=20),
        plot.background = element_blank(),
        panel.border=element_blank(),
        panel.grid.major= element_line(colour=NA), 
        panel.grid.minor=element_line(colour=NA),
        title=element_text(size=20),
        panel.background = element_rect(fill = "white"),
        axis.line.x=element_line(colour="black"),
        axis.line.y=element_line(colour="black"))
S.mfso<-mfso(~ShellFuzzyData$Type+ShellFuzzyData$Discharge.cms+
               ShellFuzzyData$ChlAdensity,
             S.Fdist, permute = 1000)
summary(S.mfso)
plot(S.mfso)
S.mfso.plot<-data.frame(Reach=ShellFuzzyData[!is.na(ShellFuzzyData$ChlAdensity)&
                                             !is.na(ShellFuzzyData$Discharge.cms),
                                             "SamID"],
                        MB.mu=S.mfso$mu[,1],
                        Dis.mu=S.mfso$mu[,2],
                        Ch.mu=S.mfso$mu[,3]) %>% left_join(ShellFuzzyData[,c(1:6)])
ggplot(S.mfso.plot, aes(x=MB.mu, y=Ch.mu, color=Type, shape=ShellSpecies))+
  geom_point(size=2)


# ACTUAL fuzzy analysis ------------
names(F.BioDens[,-c(1:6)])#those columns just identifiers, and others
FieldFuzzyDataBM<-F.BioDens %>% left_join(FieldDischarge) %>%
  left_join(FieldChlA)
FBM.Fdist<-dsvdis(FieldFuzzyDataBM[,-c(1:6,54:56)], 'bray')
FieldFuzzyDataBM$Mussel.g.m2<-replace_na(FieldFuzzyDataBM$Mussel.g.m2,0)
FieldFuzzyDataBM[FieldFuzzyDataBM$FSamID=="GL-MR.Fall2015","Discharge.cms"]<-0.0294 #usgs
FieldFuzzyDataBM[FieldFuzzyDataBM$FSamID=="GL-NM.Fall2015","Discharge.cms"]<-0.0274 #usgs
FieldFuzzyDataBM[FieldFuzzyDataBM$FSamID=="KT-MR.Fall2015","Discharge.cms"]<-0.109 
FieldFuzzyDataBM[FieldFuzzyDataBM$FSamID=="L3-MR.Summer2015","Discharge.cms"]<-0.819 
FieldFuzzyDataBM[FieldFuzzyDataBM$FSamID=="L3-NM.Summer2015","Discharge.cms"]<-0.556 

FBM.fso<-fso(FieldFuzzyDataBM$Mussel.g.m2, FBM.Fdist, permute = 10000)
summary(FBM.fso)
plot(FBM.fso)
FBM.fso.plot<-FieldFuzzyDataBM %>% mutate(mu=FBM.fso$mu) %>%
  dplyr::select(-WC_CHLA_MG.L,-Benthic_CHLA_MG.M2) %>%
  gather(Taxa, Dens, -FSamID, -Reach,-SamplingSeason,
         -Treatment,-Season,-Mussel.g.m2,-Discharge.cms,-mu) %>%
  left_join(FTaxaTable)%>%
  arrange(desc(mu))
ggplot(FBM.fso.plot, aes(x=Mussel.g.m2, y=mu))+
  geom_point(aes(color=Season), size=3)+geom_smooth(method="lm")+
  scale_x_log10()
FBM.mfso<-mfso(~FieldFuzzyDataBM$Mussel.g.m2+FieldFuzzyDataBM$Discharge.cms+
               FieldFuzzyDataBM$WC_CHLA_MG.L,
             FBM.Fdist, permute = 10000)
summary(FBM.mfso)

### Enclosure
names(E.ComPer[,-c(1:7)])#those columns just identifiers, and others
EncFuzzyData<-E.ComPer %>% left_join(EncDischarge) %>%
  left_join(EncChlAraw)
head(EncFuzzyData[,-c(1:7,53,54)])
E.Fdist<-dsvdis(EncFuzzyData[,-c(1:7,53,54)], 'bray')
E.fso<-fso(EncFuzzyData$MusselBiomass.g.m2, E.Fdist, permute = 10000)
summary(E.fso)
E.fso.plot<-EncFuzzyData %>% mutate(mu=E.fso$mu) %>%
  dplyr::select(-Discharge.cms, -ChlAdensity) %>%
  gather(Taxa, Dens, -TEid, -Enc2,-Week,-TreatA,-Type,-Spp,-MusselBiomass.g.m2,-mu) %>%
  left_join(FTaxaTable)%>%
  arrange(desc(mu))
ggplot(E.fso.plot, aes(x=MusselBiomass.g.m2, y=mu))+
  geom_point(size=3)+geom_smooth(method="lm")

greyC<-gray.colors(length(unique(E.fso.plot$Order)),
                   start = 0.02, end = 0.97, gamma = 2.2, alpha = NULL)
greyCE<-c(greyC[1:4], "red", greyC[5], "purple",greyC[6:15],"blue", greyC[17:18],
          "orange", greyC[19:22])
Emu_table <- table(E.fso.plot$mu)
Emu_levels <- names(Emu_table)[order(Emu_table)]
E.fso.plot$OrdMu<-factor(E.fso.plot$mu,levels=Emu_levels)
ggplot(data = E.fso.plot, aes(x =OrdMu, y = Dens, fill=Order)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values=greyCE, name="Species")+
  coord_flip()+
  scale_x_discrete(name="Ordinated Sites", labels=unique(E.fso.plot$TEid)) +
  theme(axis.text.x = element_text(size=9,color="black"),
        axis.title.y=element_text(size=20),
        plot.background = element_blank(),
        panel.border=element_blank(),
        panel.grid.major= element_line(colour=NA), 
        panel.grid.minor=element_line(colour=NA),
        title=element_text(size=20),
        panel.background = element_rect(fill = "white"),
        axis.line.x=element_line(colour="black"),
        axis.line.y=element_line(colour="black"))

E.mfso<-mfso(~EncFuzzyData$MusselBiomass.g.m2+EncFuzzyData$Discharge.cms+
               EncFuzzyData$ChlAdensity, E.Fdist)
summary(E.mfso)
plot(E.mfso)
E.mfso.plot<-data.frame(Reach=EncFuzzyData[!is.na(EncFuzzyData$ChlAdensity)&
                                             !is.na(EncFuzzyData$Discharge.cms),"TEid"],
                        MB.mu=E.mfso$mu[,1],
                        Dis.mu=E.mfso$mu[,2],
                        Ch.mu=E.mfso$mu[,3]) %>% left_join(EncFuzzyData[,c(1:7)])
ggplot(E.mfso.plot, aes(x=MB.mu, y=Ch.mu, color=TreatA, shape=Type))+
  geom_point()

### Shell
names(S.ComPer[,-c(1:6)])#those columns just identifiers, and others
ShellFuzzyData<-S.ComPer %>% left_join(EncDischarge[EncDischarge$Week=="w12",]) %>%
  left_join(ShellChl)
#have a -0.5 in chlorophyll, fixing
ShellFuzzyData[ShellFuzzyData$SamID=="C09.AMB","ChlAdensity"]<-0
head(ShellFuzzyData[,-c(1:6,38:41)])
S.Fdist<-dsvdis(ShellFuzzyData[,-c(1:6,38:41)], 'bray')
S.fso<-fso(~ShellFuzzyData$ChlAdensity+ShellFuzzyData$Discharge.cms+
             ShellFuzzyData$Type+ShellFuzzyData$ShellSpecies+
             ShellFuzzyData$TShellSurArea.cm2, S.Fdist,
           permute = 10000)
S.fso<-fso(~ShellFuzzyData$Type+ShellFuzzyData$ChlAdensity, S.Fdist, permute=10000)
summary(S.fso)
S.fso.plot<-ShellFuzzyData %>% mutate(muType=S.fso$mu[,1],
                                      muChlA=S.fso$mu[,2]) %>%
  dplyr::select(-Discharge.cms) %>%
  gather(Taxa, Dens, -SamID, -Enc2,-Week,-TreatA,-Type,-Spp,-ShellSpecies,
         -TShellSurArea.cm2,-ChlAdensity,-muType, -muChlA) %>%
  left_join(FTaxaTable) %>% arrange(desc(muType))
ggplot(S.fso.plot, aes(x=ChlAdensity, y=muChlA))+
  geom_point(size=3)+geom_smooth(method="lm")

greyC<-gray.colors(length(unique(S.fso.plot$Taxa)),
                   start = 0.02, end = 0.97, gamma = 2.2, alpha = NULL)
greyCS<-c(greyC[1:3], "red", greyC[7],"purple",greyC[9:12], "blue", greyC[19:22])
Smu_table <- table(S.fso.plot$muChlA)
Smu_levels <- names(Smu_table)[order(Smu_table)]
S.fso.plot$OrdMu<-factor(S.fso.plot$muChlA,levels=Smu_levels)
shell.labels<-S.fso.plot %>% group_by(SamID) %>%
  slice(1) %>% dplyr::select(SamID, Enc2, ShellSpecies, 
                             TreatA, Type, Spp, muChlA, ChlAdensity) %>%
  arrange(desc(muChlA))
ggplot(data = S.fso.plot, aes(x =OrdMu, y = Dens, fill=Order)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values=greyCS, name="Species")+
  coord_flip()+
  scale_x_discrete(name="Ordinated Sites", 
                   labels=paste(shell.labels$TreatA, round(shell.labels$ChlAdensity,2)))+
  #labels=paste(shell.labels$SamID, round(shell.labels$muType,2))) +
  theme(axis.text.x = element_text(size=9,color="black"),
        axis.title.y=element_text(size=20),
        plot.background = element_blank(),
        panel.border=element_blank(),
        panel.grid.major= element_line(colour=NA), 
        panel.grid.minor=element_line(colour=NA),
        title=element_text(size=20),
        panel.background = element_rect(fill = "white"),
        axis.line.x=element_line(colour="black"),
        axis.line.y=element_line(colour="black"))
S.mfso<-mfso(~ShellFuzzyData$Type+ShellFuzzyData$Discharge.cms+
               ShellFuzzyData$ChlAdensity,
             S.Fdist, permute = 1000)
summary(S.mfso)
plot(S.mfso)
S.mfso.plot<-data.frame(Reach=ShellFuzzyData[!is.na(ShellFuzzyData$ChlAdensity)&
                                               !is.na(ShellFuzzyData$Discharge.cms),
                                             "SamID"],
                        MB.mu=S.mfso$mu[,1],
                        Dis.mu=S.mfso$mu[,2],
                        Ch.mu=S.mfso$mu[,3]) %>% left_join(ShellFuzzyData[,c(1:6)])
ggplot(S.mfso.plot, aes(x=MB.mu, y=Ch.mu, color=Type, shape=ShellSpecies))+
  geom_point(size=2)


### exploring fuck my life ###
F.ComPer
test<-F.ComPer %>% gather(Taxa, Per, -FSamID, -Reach, - SamplingSeason, -Treatment, 
                    -Season,-Mussel.g.m2) %>%
  left_join(FTaxaTable[,1:4]) %>% group_by(FSamID, Reach, SamplingSeason, Mussel.g.m2, Order) %>%
  summarize(RA=sum(Per)) %>% spread(Order, RA)

testFData<-test %>% left_join(FieldDischarge) %>%
  left_join(FieldChlA)
testFdist<-dsvdis(testFData[,-c(1:4,27:29)], 'bray')

testFData$Mussel.g.m2<-replace_na(testFData$Mussel.g.m2,0)
test.fso<-fso(testFData$Mussel.g.m2, testFdist)
summary(test.fso)
testF.fso.plot<-testFData %>% ungroup() %>% mutate(mu=test.fso$mu) %>%
  dplyr::select(-Discharge.cms, -WC_CHLA_MG.L,-Benthic_CHLA_MG.M2) %>%
  gather(Taxa, Dens, -FSamID, -Reach,-SamplingSeason,-Mussel.g.m2,-mu) %>%
  arrange(desc(mu))
ggplot(testF.fso.plot, aes(x=Mussel.g.m2, y=mu))+
  geom_point(size=3)+geom_smooth(method="lm")

greyC<-gray.colors(length(unique(testF.fso.plot$Taxa)),
                   start = 0.02, end = 0.97, gamma = 2.2, alpha = NULL)
Fmu_table <- table(testF.fso.plot$mu)
Fmu_levels <- names(Fmu_table)[order(Fmu_table)]
testF.fso.plot$OrdMu<-factor(testF.fso.plot$mu,levels=Fmu_levels)
ggplot(data = testF.fso.plot, aes(x =OrdMu, y = Dens, fill=Taxa)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values=greyC, name="Species", guide=F)+
  coord_flip()+
  scale_x_discrete(name="Ordinated Sites", labels=unique(F.fso.plot$FSamID)) +
  theme(axis.text.x = element_text(size=9,color="black"),
        axis.title.y=element_text(size=20),
        plot.background = element_blank(),
        panel.border=element_blank(),
        panel.grid.major= element_line(colour=NA), 
        panel.grid.minor=element_line(colour=NA),
        title=element_text(size=20),
        panel.background = element_rect(fill = "white"),
        axis.line.x=element_line(colour="black"),
        axis.line.y=element_line(colour="black"))

F.mfso<-mfso(~FieldFuzzyData$Mussel.g.m2+FieldFuzzyData$Discharge.cms+
               FieldFuzzyData$WC_CHLA_MG.L,
             F.Fdist)
summary(F.mfso)
plot(F.mfso)
F.mfso.plot<-data.frame(Reach=FieldFuzzyData[!is.na(FieldFuzzyData$WC_CHLA_MG.L)&
                                               !is.na(FieldFuzzyData$Discharge.cms),"FSamID"],
                        MB.mu=F.mfso$mu[,1],
                        Dis.mu=F.mfso$mu[,2],
                        Ch.mu=F.mfso$mu[,3]) %>% left_join(FieldFuzzyData[,c(1:6)])
ggplot(F.mfso.plot, aes(x=MB.mu, y=Ch.mu, color=Treatment, shape=Season))+
  geom_point()
#### Biomass Amount #####
head(F.BioDens)
F.BioDens.plot<-F.BioDens
F.BioDens.plot$TotalBM<-rowSums(F.BioDens[,7:53])
ggplot(F.BioDens.plot, aes(x=Treatment, y=TotalBM))+
  geom_boxplot(fill="grey")+theme_krementz()+
  scale_y_log10()
ggplot(F.BioDens.plot, aes(x=Mussel.g.m2, y=TotalBM))+
  geom_point(aes(color=Season), size=3)+
  geom_smooth(method="lm", color="black", se=F)+
  theme_krementz()+NULL
  scale_y_log10()+scale_x_log10()

View(F.BioDens.plot %>% group_by(Treatment, SamplingSeason) %>%
  summarize(sumBM=round(mean(TotalBM),2), sdBM=round(sd(TotalBM),2)))

E.BioDens.plot<-E.BioDens
E.BioDens.plot$TotalBM<-rowSums(E.BioDens[,8:48], na.rm=T)
View(E.BioDens.plot %>% group_by(TreatA, Week) %>%
       summarize(sumBM=round(mean(TotalBM, na.rm=T),2), 
                 sdBM=round(sd(TotalBM, na.rm=T),2)))
S.BioDens.plot<-S.BioDens %>% dplyr::select(-Cladocera.miscD)
S.BioDens.plot$TotalBM<-rowSums(S.BioDens[,8:32], na.rm=T)*10000
View(S.BioDens.plot %>% group_by(TreatA, ShellSpecies) %>%
       summarize(sumBM=round(mean(TotalBM, na.rm=T),2), 
                 sdBM=round(sd(TotalBM, na.rm=T),2)))
