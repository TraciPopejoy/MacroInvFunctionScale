# Mantel Test for Enclosure Data ==============================
# Build a enclosure raster dataframe ------------------
library(sp)
#enclosure treatment information from CCA script
enc.tot<-CCAEnc.data[,c(1:8,41:48)] 
enc.tot$invtot<-rowSums(CCAEnc.data[,-c(1:8,41:48)]) # get total invert density
#create a roster dataframe; this has the spacing and enclosure information
EnclosureRaster <- data.frame(enc=
                                c("A01",NA,"A02",NA,"A03",NA,"A04",NA,"A05",NA,
                                  "A06",NA,"A07",NA,"A08",NA,"A09",NA,"A10",NA,
                                  rep(NA,20),
                                  NA,"B01",NA,"B02",NA,"B03",NA,"B04",NA,"B05",NA,
                                  "B06",NA,"B07",NA,"B08",NA,"B09",NA,"B10",
                                  rep(NA,20),
                                  "C01",NA,"C02",NA,"C03",NA,"C04",NA,"C05",NA,
                                  "C06",NA,"C07",NA,"C08",NA,"C09",NA,"C10",NA,
                                  rep(NA,20),
                                  NA,"D01",NA,"D02",NA,"D03",NA,"D04",NA,"D05",NA,
                                  "D06",NA,"D07",NA,"D08",NA,"D09",NA,"D10",
                                  rep(NA,20),
                                  "E01",NA,"E02",NA,"E03",NA,"E04",NA,"E05",NA,
                                  "E06",NA,"E07",NA,"E08",NA,"E09",NA,"E10",NA),
                              z = seq(180),
                              xc = c(rep(1,20),rep(2,20),rep(3,20),rep(4,20),rep(5,20),
                                     rep(6,20),rep(7,20),rep(8,20),rep(9,20)),
                              yc = c(rep(seq(1:20),9)))%>%
  left_join(enc.tot, by=c('enc'='Enc2'))
coordinates(EnclosureRaster) = ~xc+yc #promotes to SpatialPointsDataFrame
gridded(EnclosureRaster) = TRUE #promotes to SpatialPixelsDataFrame
df <- as(EnclosureRaster, "SpatialGridDataFrame") 
image(df["invtot"]) # bad plot of spatial distribution of macroinvertebrate data
cc = coordinates(df)
a=df[["TreatA"]]
zc=as.character(a)
zc[is.na(zc)]="." #keeps NA by putting . for text call below
text(cc[,1],cc[,2],zc) #adds the enclosure names to image above

# Mantel test to look at spatial autocorrelation --------------------
# remove NAs (blanks to get proper spatial orientation)
ERred<-EnclosureRaster[!is.na(EnclosureRaster$enc),] 
spat.dists <- dist(coordinates(ERred)) #calculate euclidian distance of enclosure space
inv.dists <- dist(ERred$invtot) #calculate euclidian distance of invertebrate samples
library(ade4)
mantel.rtest(spat.dists, inv.dists, nrepet = 9999)
#marginally not significant
#null hypothesis rejected = little spatial correlation

# Fancy graph for Supplement ----------------------
# follows this tutorial:  https://timogrossenbacher.ch/2016/12/beautiful-thematic-maps-with-ggplot2-only/
#bind dataframe without NA and then use fortify to make it long
ERredF<-fortify(cbind(ERred@data, ERred@coords)) 
labelsX <- c()

#build quantiles to maximize color variation within graph
quantiles <- quantile(ERredF$invtot, 
                      probs = seq(0, 1, length.out = 6 + 1))

# here I define custom labels (the default ones would be ugly)
for(idx in 1:length(quantiles)){
  labelsX <- c(labelsX, paste0(round(log10(quantiles[idx]), 2), 
                               " – ", 
                               round(log10(quantiles[idx + 1]), 2)))
}
# I need to remove the last label 
# because that would be something like "66.62 - NA"
labelsX <- labelsX[1:length(labelsX)-1]

# here I actually create a new variable on the dataset with the quantiles
ERredF$invquant <- cut(ERredF$invtot, 
                       breaks = quantiles, 
                       labels = labelsX, 
                       include.lowest = T)

# create the good map - Figure S1
ggplot(ERredF, aes(x=xc, y=yc, fill=invquant)) + 
  geom_raster()+
  geom_text(aes(label=TreatA), color="white", fontface="bold")+
  scale_y_reverse()+ scale_x_reverse()+
  scale_fill_viridis_d(name = expression(paste("Log"[10]," ind. m"^-2)), 
                       direction = -1,
                       guide = guide_legend(keyheight = unit(5, units = "mm"),
                                            title.position = 'top', reverse = T))+ 
  theme_minimal()+
  ylab("Downstream               Upstream") + 
  xlab("North Bank             South Bank")
ggsave("Figures/S1Mantelinvertbiomassenc.jpg", width=5, height=6)

# CCA for Supplement --------------
enc.envR<-CCAEnc.data[,c(1:8, 41:48)] %>% ungroup()%>%
  left_join(ERredF[,c(1,19,20)], by=c("Enc2"="enc"))%>% 
  mutate_if(is.numeric,scale) %>% 
  #changing type into binary to reduce variable redundancy
  mutate(Live=case_when(Type=="Control"~"No",
                        Type=="Sham"~"No",
                        Type=="Live"~"Yes")) 

# run the cca
cca.E.cond<-cca(enc.com~ACT+AMB+Live+
             Discharge.cms+ChlAdensity+Dvar+Condition(xc+yc), enc.envR, scaling=2)

vif.cca(cca.E.cond) #check for redundant environmental variables
cca.E.cond
ccaE.plot<-plot(cca.E)
cca.E.sum<-summary(cca.E)
