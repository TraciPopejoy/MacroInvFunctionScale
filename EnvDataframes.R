library(tidyverse); library(readxl)
### covariables - environment and spatial
# Spatial Information =================
FieldSpData<-read_excel("./data/HotSpotLocations.xlsx") %>% 
  mutate(Reach=paste(SiteID, Treatment, sep="-")) %>% 
  dplyr::select(-`Drainage - HUC 8`)

# HUC12 (or HUC8) for spatial autocorrelation 
library(sf)
SEOKhu12<- read_sf(dsn="C:/Users/Owner/Documents/GISfile/SEOkNHD/Shape", layer="WBDHU12")
SEOKhu8<- read_sf(dsn="C:/Users/Owner/Documents/GISfile/SEOkNHD/Shape", layer="WBDHU8")

FieldSpData<-st_as_sf(FieldSpData, 
                      coords = c('Longitude', 'Latitude'), 
                      crs=st_crs(SEOKhu12))
FieldSpData<-st_join(FieldSpData, SEOKhu8 %>% dplyr::select(HUC8, NAME)) %>%
  st_join(SEOKhu12 %>% dplyr::select(HUC12)) %>%
  mutate(HUC12num=as.numeric(paste(HUC12)),
         HUC8num=as.numeric(paste(HUC8)))

library(dataRetrieval)
gagesites <- readNWISdata(site=c("07337900","07338500","07335790","07335700"), 
                          service = "site")
coordinates(gagesites)<-c("dec_long_va", "dec_lat_va")

#figure out distances between points
FieldSpData[,1:3]
library(raster) #too lazy to code it up properly, brute force
distancesbtwnR<-c('K7'=pointDistance(FieldSpData[1,], FieldSpData[2,],longlat=T),
'KS'=pointDistance(FieldSpData[3,], FieldSpData[4,],longlat=T),
'KT'=pointDistance(FieldSpData[5,], FieldSpData[6,],longlat=T),
'L3'=pointDistance(FieldSpData[7,], FieldSpData[8,],longlat=T),
'LY'=pointDistance(FieldSpData[9,], FieldSpData[10,],longlat=T),
'K2'=pointDistance(FieldSpData[11,], FieldSpData[12,],longlat=T),
'GL'=pointDistance(FieldSpData[13,], FieldSpData[14,],longlat=T))
distancesbtwnR
mean(distancesbtwnR);sd(distancesbtwnR)
min(distancesbtwnR);max(distancesbtwnR)

#Discharge ===========================
FieldDischarge<-read_excel("./data/Field Discharge.xlsx") %>% 
  group_by(Reach, SamplingSeason) %>%
  summarize(Discharge.cms=mean(`Discharge (M3/s)`, na.rm=T))
####suplementing values find a better way! ###
#based on F16 data from respective reach (or S16 for L3)
FieldDischarge %>% left_join(FieldSpData %>% as_tibble()) %>%
  filter(SamplingSeason=="Fall2016") %>%
  ungroup() %>% group_by(River,`USGS gage Identifier`)%>%
  summarize(median(Discharge.cms), IQR(Discharge.cms))

EncDVw09<-read.csv("../FEn17/FEn17_data/EncPhysDisFEn17OK.csv") %>% 
  rename(Enclosure=`ï..Enclosure`) %>% left_join(TreatENC) %>%
  mutate(Discharge.cms=0.25*V.mps*Depth.m,
         Week="w09") %>%
  dplyr::select(Enc2,Week,Discharge.cms)
EncDVw12<-read_excel("../FEn17/FEn17_data/videowithabiotics.xlsx") %>% 
  left_join(TreatENC, by=c("Unit"="Enclosure")) %>%
  mutate(Discharge.cms=0.25*Velocity*Depth,
         Week="w12") %>% filter(Month=="Oct") %>%
  dplyr::select(Enc2,Week,Discharge.cms)
EncDischarge<-rbind(EncDVw09,EncDVw12)
EncDtran<-read_excel("../FEn17/FEn17_data/FieldEncDataSum2017V2.xlsx", 
                     sheet="DischargeT") %>%
  filter(Date >"2017-09-17 00:00:00") 
EncDtran %>% summarize(meanSeptD=mean(totalQ.cms, na.rm=T))

# Chlorophyll Abundance ===============
FieldChlA<-read_excel("./data/Field Chlorophyll.xlsx") %>% slice(-14)
#K2 MR Sum2015 was duplicated, don't know why
EncChlAraw<-read_csv("../FEn17/FEn17_data/ChlDataFEn17OK.csv") %>%
  dplyr::rename(WeekBad = Week, Enclosure = `Enclosure#`) %>%
  mutate(Week=case_when(WeekBad==9~"w09",
                        WeekBad==4~"w04",
                        WeekBad==12~"w12")) %>%
  filter(Type=="Tiles") %>% 
  dplyr::select(-Type) %>%left_join(TreatENC, by=("Enclosure")) %>%
  mutate(Area=5.939574, # centimeter squared
         ChlAdensity=26.7*((`664nm`-fir750nm)-(`665nm`-sec750nm))*(Vacetone/Area)*1,
         TEid=paste(Week,Enc2, sep="")) %>%
  dplyr::select(TEid,Week,Enc2,TreatA,ChlAdensity)
head(MDat)
SlurryData<-read_csv("../FEn17/FEn17_data/SubSamFEn17OK.csv") %>%
  dplyr::rename(WeekBad = Week, Filter=`Filter#`) %>%
  mutate(Week=case_when(WeekBad==8~"w09",
                        WeekBad==4~"w04",
                        WeekBad==12~"w12"))
SlurryData$BucketVol<-as.numeric(paste(SlurryData$BucketVol)) #tell volume to stop being a factor/text
ShellChl<-read_csv("../FEn17/FEn17_data/ChlDataFEn17OK.csv") %>%
  dplyr::rename(WeekBad = Week, Enclosure = `Enclosure#`) %>%
  mutate(Week=case_when(WeekBad==9~"w09",
                        WeekBad==4~"w04",
                        WeekBad==12~"w12"))%>%
  filter(Type=="shellsampling") %>%  dplyr::select(-Type,-WeekBad,-Week)%>%
  left_join(SlurryData, by=c("Enclosure"="Filter"))%>%
  dplyr::select(-`Basket#`)%>%dplyr::rename(Filter=Enclosure)%>%
  left_join(TreatENC, by=c("Enclosure.y"="Enclosure")) %>%
  mutate(SamID=paste(Enc2, substr(Type.x,3,6),sep="."))
ShellChl$SamID<-recode(ShellChl$SamID, C10.0ACT="C10.ACT",C10.0AMB="C10.AMB",
                       E10.0ACT="E10.ACT",E10.0AMB="E10.AMB")
ShellChl<-ShellChl %>% left_join(MDat[,3:4], by="SamID") %>%
  mutate(BucketVol2=as.numeric(paste(BucketVol)),
         ChlAdensity=26.7*((`664nm`-fir750nm)-(`665nm`-sec750nm))*(Vacetone/VolFilrwe)*(BucketVol2/TShellSurArea.cm2)*1)%>%
  dplyr::select(SamID,Enc2,ChlAdensity, TShellSurArea.cm2)

# Pebble Counts ==============
peb.raw<-read_excel("./data/Pebble Counts.xlsx", sheet = "Pebble Counts Reprocessed")
#install_github("bceaton/GSDtools")
library(GSDtools)
peb.wolfman<-peb.raw[,1:3]  %>% rename(Treatment = Reach) %>% 
  mutate(Reach=paste(SiteID, Treatment, sep="-"))
for(k in 1:nrow(peb.raw)){
  store<-MakeCFD(as.matrix(peb.raw[k,4:104]))
  wolf<-WolmanCI(store, n=100, P=c(10,50,90,60))
  sub<-data.frame(D10=wolf[1,2],
           D10low=wolf[1,3],
           D10high=wolf[1,4],
           D50=wolf[2,2],
           D50low=wolf[2,3],
           D50high=wolf[2,4],
           D90=wolf[3,2],
           D90low=wolf[3,3],
           D90high=wolf[3,4],
           D60=wolf[4,2])
  peb.wolfman[k,5:14]<-sub
}

ggplot(store, aes(x=size, y=probs))+geom_point()+geom_line()+
  scale_x_log10(breaks=c(0.1,0.3,1,3,10,30,100))

peb.enc<-read_excel("./data/Enc_pebbles.xlsx", sheet = "fixed")
peb.enc.sum<-peb.enc[,1:2]  
for(k in 1:nrow(peb.enc)){
  store<-MakeCFD(as.matrix(peb.enc[k,3:27]))
  wolf<-WolmanCI(store, n=25, P=c(10,50,90,60))
  sub<-data.frame(D10=wolf[1,2],
                  D10low=wolf[1,3],
                  D10high=wolf[1,4],
                  D50=wolf[2,2],
                  D50low=wolf[2,3],
                  D50high=wolf[2,4],
                  D90=wolf[3,2],
                  D90low=wolf[3,3],
                  D90high=wolf[3,4],
                  D60=wolf[4,2]
                  )
  peb.enc.sum[k,3:12]<-sub
}

peb.ENC.sum<-peb.enc.sum %>% mutate(Dvar=D60/D10) %>%left_join(TreatENC)
ggplot(peb.ENC.sum, aes(x=TreatA, y=D50))+geom_boxplot()
ggplot(peb.ENC.sum, aes(x=TreatA, y=Dvar))+geom_boxplot()

# Summary =========
#joined data frames
Fenv.data<-FieldChlA %>% full_join(FieldDischarge) %>% 
  left_join(FieldSpData %>% as_tibble()) %>%
  left_join(peb.wolfman) %>% mutate(Dvar=D60/D10) %>%
  dplyr::select(Reach,SamplingSeason,WC_CHLA_MG.L,Benthic_CHLA_MG.M2,
         Discharge.cms,HUC8num,HUC12num,Dvar,D50,D90)

Eenv.data<-EncChlAraw %>% left_join(EncDischarge) %>% 
  left_join(TreatENC) %>% left_join(peb.ENC.sum) %>%
  dplyr::select(TEid,Week,Enc2,TreatA,ChlAdensity,Discharge.cms,Type,Spp,
         Dvar, D50, D90)
