# Code to figure out USGS averages at different seasons =======
library(dataRetrieval);library(lubridate);library(readxl)
# Hydrology Analysis ----------------------
# parameter 00060 is cfs daily mean
discharge <- readNWISdata(site=c("07337900","07338500","07335790","07335700"), 
                    parameterCd = c("00060"),
                    startDate="1960-01-01", endDate="2018-12-30") %>%
  mutate(Gage=case_when(site_no=="07337900"~"Glover",
                        site_no=="07338500"~"Little@Lukfata",
                        site_no=="07335790"~"Kimichi@Clayton",
                        site_no=="07335700"~"Kiamichi@BigCedar"),
         Date=ymd(dateTime),
         Month=month(Date),
         Year=year(Date),
         Day=day(Date),
         MonDay=paste(Month,Day, sep="."),
         Season=case_when(Month==11|Month==12|Month==1|Month==2~"Winter",
                          Month==3|Month==4|Month==5~"Spring",
                          Month==6|Month==7|Month==8~"Summer",
                          Month==9|Month==10~"Fall"))
# just look at the fall data 
View(discharge %>% group_by(Gage,Season) %>% 
  summarise(meanD=mean(X_00060_00003)*0.0283168466, #convert to cms
            sdD=sd(X_00060_00003)*0.0283168466,
            medD=median(X_00060_00003)*0.0283168466,
            iqrD=IQR(X_00060_00003)*0.0283168466,
            minD=min(Date),
            maxD=max(Date)) %>%
  filter(Season == "Fall")) 
min(discharge$dateTime);max(discharge$dateTime)

# goal: figure out % above or below during my experiment
# downloaded the % amounts for each cms value at each gage from the USGS website
perStreamCond<-read_excel("./data/USGS dvstat.xlsx", sheet=2) %>%
  mutate(MonDay=paste(month_nu,day_nu, sep="."),
         siteNO=as.character(paste("0",site_no, sep="")))
# join our discharge data to the perStreamCondition dataframe 
monster<-left_join(discharge, perStreamCond, by=c("site_no"="siteNO","MonDay"))

#make a function to predict where the recorded discharges falls on the % rate 
testRate<-function(disCurve, discharge){
  testrate<-data.frame(per=c(5,10,20,25,50,75,80,90,95),dis=as.numeric(disCurve))
  line<-lm(log1p(dis)~per,data=testrate)
  x<-(log1p(discharge)-as.numeric(coefficients(line)[1]))/as.numeric(coefficients(line)[2])
}

# apply the function to each discharge measurement to figure % value that day
monsterDis<-monster %>% rowwise() %>% 
  mutate(StreamCon=testRate(c(p05_va,p10_va,p20_va,p25_va,
                              p50_va,p75_va,p80_va,p90_va,p95_va),X_00060_00003))

# summarize Stream Condition for each year
# so this reports the % of time flow is above a certain flow
DisSumYear<-monsterDis %>% group_by(Gage, Year) %>%
  summarize(meanYearDis.cfs=mean(X_00060_00003))
DisSumSeason<-monsterDis %>% group_by(Gage, Year, Season) %>%
  summarize(meanDis.cfs=mean(X_00060_00003),
            meanStreamCon.per=mean(StreamCon))

# make a big table of both raw discharge values and stream condition
covTableDis<-left_join(DisSumYear,DisSumSeason) %>% select(everything(), meanYearDis.cfs) %>%
  filter(Season!="Winter" & Season!="Spring")
write.csv(covTableDis, file="DischargeSummary.csv") #save that tabel

# Pulling out just the USGS discharge during enclosure experiment
FEn17Dis<-monsterDis %>% 
  filter(Gage=="Kiamichi@BigCedar", Date < "2017-10-08", Date >"2017-05-01")%>%
  mutate(timeperiod=case_when(Date <="2017-07-02"& Date >"2017-05-01"~"start",
                                Date <="2017-8-10"& Date >"2017-07-02"~"w04",
                                Date <="2017-9-15"& Date >"2017-08-10"~"w09",
                                Date <="2017-10-08"& Date >"2017-09-15"~"w12"))%>%
  group_by(Gage, Year, timeperiod)%>% 
  summarize(meanYearDis.cfs=mean(X_00060_00003),
            meanStreamCon.per=mean(StreamCon))


# Temperature Analysis ---------------------------------------------
# 00010 is water temperature, C
# use readNWISdata to pull this data into R from USGS website
# all times available
temp12<-readNWISdata(site=c("07337900","07338500","07335790","07335700"), 
                          parameterCd = c("00010"),
                          startDate="1984-01-01", endDate="2018-12-30",
                          service="dv") %>%
  mutate(Gage=case_when(site_no=="07337900"~"Glover",
                        site_no=="07338500"~"Little@Lukfata",
                        site_no=="07335790"~"Kimichi@Clayton",
                        site_no=="07335700"~"Kiamichi@BigCedar"),
         Date=ymd(dateTime),
         Month=month(Date),
         Year=year(Date),
         Season=case_when(Month==12|Month==11|Month==1|Month==2~"Winter",
                          Month==3|Month==4|Month==5~"Spring",
                          Month==6|Month==7|Month==8~"Summer",
                          Month==9|Month==10~"Fall"))
						  
# limiting the search to around midday to avoid extreme temperatures
tempear<-readNWISqw(site=c("07337900","07338500","07335790","07335700"), 
             parameterCd = c("00010"),
             startDate="1950-01-01", endDate="2018-12-30") %>%
  select("agency_cd","sample_dt","sample_tm","site_no","parm_cd","result_va") %>%
  filter(sample_tm > "11:00" | is.na(sample_tm)) %>%
  mutate(Gage=case_when(site_no=="07337900"~"Glover",
                        site_no=="07338500"~"Little@Lukfata",
                        site_no=="07335790"~"Kimichi@Clayton",
                        site_no=="07335700"~"Kiamichi@BigCedar"),
         Date=ymd(sample_dt),
         Month=month(Date),
         Year=year(Date),
         Season=case_when(Month==12|Month==1|Month==2~"Winter",
                          Month==3|Month==4|Month==5~"Spring",
                          Month==6|Month==7|Month==8~"Summer",
                          Month==9|Month==10|Month==11~"Fall")) %>%
  rename("X_00010_00003"="result_va")
  
# joining the two into one table  
tempUSGS<-full_join(temp12, tempear)
tempUSGS %>% group_by(Gage,Season) %>% 
  summarise(meanT=mean(X_00010_00003),
            medT=median(X_00010_00003),
            iqrT=IQR(X_00010_00003),
            minD=min(Date),
            maxD=max(Date)) %>%
  filter(Season == "Fall")
min(temperature$dateTime);max(temperature$dateTime)

# building the table for the manuscript
TempSumSeason<-temperature %>% group_by(Gage, Year, Season) %>%
  summarize(meanTemp.C=mean(X_00010_00003)) %>%
  spread(Season,meanTemp.C)
#covTableDis<-left_join(DisSumYear,DisSumSeason)
#write.csv(covTableDis, file="DischargeSummary.csv")

#pulling in the enclosure data to compare it
EncPC<-read_excel("../FEn17/FEn17_data/FieldEncDataSum2017V2.xlsx", 
                     sheet="Physical-Chem") 
EncPC %>% filter(Date >"2017-09-17 00:00:00", Variable=="Temperature")
