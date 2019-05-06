### covariables - environment and spatial
# Spatial Information =================
FieldSpData<-read_excel("./data/HotSpotLocations.xlsx") %>% 
  mutate(Reach=paste(SiteID, Treatment, sep="-"))
library(sp)
coordinates(FieldSpData)<-c("Longitude","Latitude")
plot(FieldSpData)

# Discharge ===========================
FieldDischarge<-read_excel("./data/Field Discharge.xlsx") %>% 
  group_by(Reach, SamplingSeason) %>%
  summarize(Discharge.cms=mean(`Discharge (M3/s)`, na.rm=T))
####suplementing values find a better way! ###
#based on F16 data from respective reach (or S16 for L3)

EncDVw09<-read.csv("../FEn17//FEn17_data/EncPhysDisFEn17OK.csv") %>% 
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

