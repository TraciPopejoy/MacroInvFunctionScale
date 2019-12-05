library(ggplot2)
library(hydroMap);library(dataRetrieval);library(raster)
library(maptools); library(maps)
# starting from line 28 in EnvDataframes.R
# good data frames already: gagesites, FieldSpData
# shape files for the waterbasins: SEokHuc8, SEokHuc12
plot(SEOKhu8[SEOKhu8$HUC8 %in% FieldSpData$HUC8,])
goodWB<-SEOKhu8[SEOKhu8$HUC8 %in% FieldSpData$HUC8,]

#grab all streams 4 order and above within the bounding box
Stream4 <- getFlowLines(c(-96, -94.6, 33.6,35),
                            streamorder = 4)
#pull the specific rivers 
riv<-c("Kiamichi River","Glover River","Little River")
ImpRiv<-Stream4[Stream4@data$gnis_name %in% riv,]

#get reservoirs
klakes<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/KiamichiShape", 
                 layer="NHDWaterbody")
llakes<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/LittleShape", 
                layer="NHDWaterbody")
mlakes<-readOGR(dsn="C:/Users/Owner/Documents/GISfile/MountainForkShape", 
        layer="NHDWaterbody")
names(mlakes)<-c("objectID","PERMANENT_", "FDATE",
                 "RESOLUTION", "GNIS_ID","GNIS_NAME",
                 "AREASQKM", "ELEVATION", "REACHCODE",
                 "FTYPE", "FCODE","visibility",
                 "SHAPE_LENG", "SHAPE_AREA")
kiares<-spRbind(klakes[grep('Hugo', klakes$GNIS_NAME),],
klakes[grep('Sardis',klakes$GNIS_NAME),])
lires<-spRbind(mlakes[640,-c(1,12)],
               llakes[1904,])
oklakes<-spRbind(kiares, lires)

#check I have what I need
plot(goodWB, col='darkgrey')
plot(Stream4, col='lightblue', add=T)
plot(ImpRiv, add=T)
plot(oklakes, add=T, col='lightblue')

library(ggmap)
#these calls 'tidy' the data so I can plot it with ggplot
#usa <- map_data("usa")
#states <- map_data("state")
wbtidy<-broom::tidy(goodWB) #lots of warnings about character strings
OKflow<-broom::tidy(Stream4) #lots of warnings about character strings
oklakestidy<-broom::tidy(oklakes)
imp<-broom::tidy(ImpRiv)
MB<-as_tibble(FieldSpData@data) %>% cbind(FieldSpData@coords)%>%
  mutate(Lab=case_when(Reach=="ENC-NA"~"ENC",
                        T~paste(Reach)))
tidygage<-as_tibble(gagesites@data) %>% cbind(gagesites@coords) %>%
  rename("SiteID"="site_no","Longitude"="dec_long_va","Latitude"="dec_lat_va") %>%
  mutate(Treatment="gage", Lab=NA)
points<-rbind(MB[,c(1,2,10:12)], tidygage[,c(2,11:14)])
points[is.na(points$Treatment),2]<-"Enc."
points$treFac<-factor(points$Treatment, levels=c("MR","NM","Enc.","gage"))
#devtools::install_github("3wen/legendMap")
library(legendMap); library(ggrepel)

okmap<-ggplot()+
  geom_path(data=OKflow, aes(x=long, y=lat, group=group), 
            color="lightblue")+
  geom_path(data=imp, aes(x=long, y=lat, group=group),
            color="darkblue")+
  geom_path(data = wbtidy, aes(x = long, y = lat, group=group),
            color="darkgrey")+
  geom_polygon(data=oklakestidy, aes(x=long,y=lat,group=group), 
               fill="lightblue")+
  #geom_polygon(data=states, aes(x=long, y=lat, group=group), 
  #             fill=NA, color="darkgrey")+
  geom_point(data=points,aes(x=Longitude, y=Latitude, 
                 shape=treFac),
             size=2)+
  xlab("Longitude")+ylab("Latitude")+
  coord_map(xlim=c(-95.86, -94.42),ylim=c(33.80,34.85))+theme_bw()+
  scale_shape_manual(name="", guide = 'legend',
                       values=c(24,25,3,15),
                       labels = c('Mussel','Control','Enc.','USGS Gage')) +
  scale_bar(lon = -95.80, lat = 33.85, 
            distance_lon = 10, distance_lat = 3, 
          distance_legend = 6, dist_unit = "km", legend_size=2.2,
          arrow_length = 10, arrow_distance=9, arrow_north_size = 4)+
  theme(panel.border = element_rect(colour = "black", fill=NA),
        panel.background = element_rect(fill = NA, colour = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(), axis.title = element_blank(),
        rect = element_blank(),
        plot.margin = unit(0 * c(-1.5, -1.5, -1.5, -1.5), "lines"),
        legend.position="bottom")
okmap
ggsave("Map.tiff", okmap,dpi=300, width=3.5, height=3.5)

okmap+geom_text_repel(data=points[!is.na(points$Lab),],
                      aes(x=Longitude, y=Latitude, 
                          label=Lab),  size=2.2)
ggsave("Maplabel.tiff",dpi=300, width=3.5, height=3.5)

leftmap<-plot_grid(okmap, bigmap, rel_heights =c(.75,.25), ncol=1)
totalMap<-plot_grid(leftmap,gamap, nrow=1, labels="AUTO")
ggsave("testingfuck.tiff", totalMap, height=3.5)
streams<-plot_grid(okmap,gamap, nrow=1, labels="AUTO")
ggsave("streams.tiff",streams,height=3.5)

library(gridExtra)
p1<-arrangeGrob(bigmap, okmap, gamap, layout_matrix = rbind(c(2,3),c(1,3)),
                respect=TRUE, heights=c())
ggsave("testmap.tiff",p1)

library(LifeTables)
data(MLTobs) 
test.mx.m <- mlt.mx[,1]
# build the life table 
lt.mx(nmx=test.mx.m, sex="male")
