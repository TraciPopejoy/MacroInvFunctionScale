
#flowLines <- getFlowLines(mapRange, streamorder, filePath)

##### Hydrology Analysis #####
library(dataRetrieval);library(lubridate);library(readxl)
# 00060 is cfs daily mean
discharge <- readNWISdata(site=c("07337900","07338500","07335790","07335700"), 
                    parameterCd = c("00060"),
                    startDate="2015-01-01", endDate="2017-12-30") %>%
  mutate(Gage=case_when(site_no=="07337900"~"Glover",
                        site_no=="07338500"~"Little@Lukfata",
                        site_no=="07335790"~"Kimichi@Clayton",
                        site_no=="07335700"~"Kiamichi@BigCedar"),
         Date=ymd(dateTime),
         Month=month(Date),
         Year=year(Date),
         Day=day(Date),
         MonDay=paste(Month,Day, sep="."),
         Season=case_when(Month==12|Month==1|Month==2~"Winter",
                          Month==3|Month==4|Month==5~"Spring",
                          Month==6|Month==7|Month==8~"Summer",
                          Month==9|Month==10|Month==11~"Fall"))
perStreamCond<-read_excel("dvstat.xlsx", sheet=2) %>%
  mutate(MonDay=paste(month_nu,day_nu, sep="."),
         siteNO=as.character(paste("0",site_no, sep="")))

monster<-left_join(discharge, perStreamCond, by=c("site_no"="siteNO","MonDay"))

testRate<-function(disCurve, discharge){
  testrate<-data.frame(per=c(5,10,20,25,50,75,80,90,95),dis=as.numeric(disCurve))
  line<-lm(log1p(dis)~per,data=testrate)
  x<-(log1p(discharge)-as.numeric(coefficients(line)[1]))/as.numeric(coefficients(line)[2])
}

monsterDis<-monster %>% rowwise() %>% 
  mutate(StreamCon=testRate(c(p05_va,p10_va,p20_va,p25_va,
                              p50_va,p75_va,p80_va,p90_va,p95_va),X_00060_00003))

DisSumYear<-monsterDis %>% group_by(Gage, Year) %>%
  summarize(meanYearDis.cfs=mean(X_00060_00003))
DisSumSeason<-monsterDis %>% group_by(Gage, Year, Season) %>%
  summarize(meanDis.cfs=mean(X_00060_00003),
            meanStreamCon.per=mean(StreamCon))


covTableDis<-left_join(DisSumYear,DisSumSeason) %>% select(everything(), meanYearDis.cfs) %>%
  filter(Season!="Winter" & Season!="Spring")
write.csv(covTableDis, file="DischargeSummary.csv")

FEn17Dis<-monsterDis %>% 
  filter(Gage=="Kiamichi@BigCedar", Date < "2017-10-08", Date >"2017-05-01")%>%
  mutate(timeperiod=case_when(Date <="2017-07-02"& Date >"2017-05-01"~"start",
                                Date <="2017-8-10"& Date >"2017-07-02"~"w04",
                                Date <="2017-9-15"& Date >"2017-08-10"~"w09",
                                Date <="2017-10-08"& Date >"2017-09-15"~"w12"))%>%
  group_by(Gage, Year, timeperiod)%>% 
  summarize(meanYearDis.cfs=mean(X_00060_00003),
            meanStreamCon.per=mean(StreamCon))

# 00010 is cfs daily mean
temperature<-readNWISdata(site=c("07337900","07338500","07335790","07335700"), 
                          parameterCd = c("00010"),
                          startDate="2015-01-01", endDate="2017-12-30",
                          service="dv") %>%
  mutate(Gage=case_when(site_no=="07337900"~"Glover",
                        site_no=="07338500"~"Little@Lukfata",
                        site_no=="07335790"~"Kimichi@Clayton",
                        site_no=="07335700"~"Kiamichi@BigCedar"),
         Date=ymd(dateTime),
         Month=month(Date),
         Year=year(Date),
         Season=case_when(Month==12|Month==1|Month==2~"Winter",
                          Month==3|Month==4|Month==5~"Spring",
                          Month==6|Month==7|Month==8~"Summer",
                          Month==9|Month==10|Month==11~"Fall"))
TempSumSeason<-temperature %>% group_by(Gage, Year, Season) %>%
  summarize(meanTemp.C=mean(X_00010_00003)) %>%
  spread(Season,meanTemp.C)
covTableDis<-left_join(DisSumYear,DisSumSeason)
write.csv(covTableDis, file="DischargeSummary.csv")
