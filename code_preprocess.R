library(dplyr)
library(reshape)
library(stringr)
setwd("~/OneDrive - University of Pittsburgh/XiaLab/ADBiomarker/manuscript/20240903/")
data=read.csv("data/rawdata.csv")
DATA=data
#### Basic information(nrow=3200,n_id=3200) ####
data=data[,1:136]
dim(data)
n_miss=apply(data[,3:136], MARGIN=1, FUN=function(x) {sum(is.na(x))+sum(x=="",na.rm = T)})
table(n_miss)
data=data[n_miss!=134,]
dim(data); length(unique(data$id_participant)) ####3885 ppts
data$dob=as.Date(data$dob)
data$date_firstsx=as.Date(data$date_firstsx)
data$date_msdx=as.Date(data$date_msdx)
write.csv(data, "data/BasicInformation.csv",row.names = F)

#### PDDS ####
data=DATA[,c(1,2,1168:1220)]
data=data %>% filter(!is.na(pdds_number))
table(data$redcap_event_name,useNA = "ifany")
name1=colnames(data)[grepl("date",colnames(data))==T][1:25]
name2=colnames(data)[grepl("value",colnames(data))==T ]
data1=data %>% dplyr::select(c(id_participant,redcap_event_name,pdds_number,any_of(name1)))
long1=melt(data1,id.vars = c("id_participant","redcap_event_name","pdds_number"))
long1=long1[,-4]
colnames(long1)[4]="pdds_date"
data2=data %>% dplyr::select(c(id_participant,redcap_event_name,pdds_number,any_of(name2)))
long2=melt(data2,id.vars = c("id_participant","redcap_event_name","pdds_number"))
long2=long2[,-4]
colnames(long2)[4]="PDDS"
long=long1
long$PDDS=long2$PDDS
long=long%>% arrange(id_participant)
long=long %>%  group_by(id_participant)%>% mutate(PDDS_n=1:length(id_participant)) 
long=long[,c(1,2,3,6,4,5)]
long=long %>% filter(PDDS_n<=pdds_number)
dim(long); length(unique(long$id_participant)) #13,553 PDDS; 2063 ppts
write.csv(long, "data/pdds.csv",row.names = F)
#### EDSS ####
data=DATA[,c(1,2,911,916)]
data=data %>% filter(is.na(edss_score)==F)
write.csv(data, "data/edss.csv",row.names = F)
### Functional test ###
data=DATA[,c(1,2,834:874)]
data<- data %>% filter( sdmt_date!="" |pasat_date!=""| hpt_date!=""| footwalk_date!="")
dim(data); length(unique(data$id_participant)) ###2199 rows, 1358 ppts
write.csv(data, "data/FunctionalTest.csv",row.names = F)
#### OCT ####
names=colnames(DATA)[grepl("oct",colnames(DATA))==T]
feature=names[6:19]
data=DATA[,c("id_participant","report1_date",feature)]
data=data %>% mutate(report1_date=as.Date(report1_date),
                     oct_n=1) %>%
  filter(is.na(report1_date)==F)
summary(data$report1_date)
data$report1_date[which(data$report1_date=="1015-08-24")]="2015-08-24"
data$report1_date[which(data$report1_date=="1026-09-22")]="2016-09-22"
dim(data)
oct=data
data=DATA[,c("id_participant",paste0("report1_date_",2),paste0(feature,"_2"))]
colnames(data)=c("id_participant","report1_date",feature)
data=data %>% mutate(report1_date=as.Date(report1_date),
                     oct_n=2) %>%
  filter(is.na(report1_date)==F)
summary(data$report1_date)
dim(data)
oct=rbind(oct, data)
for (i in 3:24){
  data=DATA[,c("id_participant",paste0("report1_date_",i),paste0(feature,"_",i))]
  colnames(data)=c("id_participant","report1_date",feature)
  data=data %>% mutate(report1_date=as.Date(report1_date),
                       oct_n=i) %>%
    filter(is.na(report1_date)==F)
  print( summary(data$report1_date))
  print( dim(data))
  oct=rbind(oct, data)
}
table(oct$oct_n)
oct=oct %>% arrange(id_participant,report1_date)
oct=oct[,c(1,2,17,3:16)]
oct=oct %>% mutate(rnfl_left_oct=ifelse(rnfl_left_oct==9999, NA, rnfl_left_oct),
                   rnfl_right_oct=ifelse(rnfl_right_oct==9999, NA, rnfl_right_oct),
                   macular_left_oct=ifelse(macular_left_oct==9999, NA, macular_left_oct),
                   macular_right_oct=ifelse(macular_right_oct==9999, NA, macular_right_oct),
                   mac_thickness_left_oct=ifelse(mac_thickness_left_oct==9999, NA, mac_thickness_left_oct),
                   mac_thickness_right_oct=ifelse(mac_thickness_right_oct==9999, NA, mac_thickness_right_oct),
                   min_gcl_left_oct=ifelse(min_gcl_left_oct==9999, NA, min_gcl_left_oct),
                   min_gcl_right_oct=ifelse(min_gcl_right_oct==9999, NA, min_gcl_right_oct),
                   average_gcl_left_oct=ifelse(average_gcl_left_oct==9999, NA, average_gcl_left_oct),
                   average_gcl_right_oct=ifelse(average_gcl_right_oct==9999, NA, average_gcl_right_oct),
                   maximum_gcl_left_oct=ifelse(maximum_gcl_left_oct==9999, NA, maximum_gcl_left_oct),
                   maximum_gcl_right_oct=ifelse(maximum_gcl_right_oct==9999, NA, maximum_gcl_right_oct) )
write.csv(oct, "data/oct.csv",row.names = F)

#### Neuroquant ####
data=read.csv("data/MRINeuroQuant.csv")
data=data %>% mutate(TBV= as.numeric(brain_vol),
                     TBV_per=brain_per,
                     ICV=total_icv,
                     GM=left_gray_vol+right_gray_vol,
                     GM_per= left_gray_per+right_gray_per,
                     TM=total_thal_vol,
                     TM_per=total_thal_per ) %>%
  select(c(mri_date, id_promote , ICV, TBV, TBV_per, GM, GM_per, TM, TM_per  )) %>%
  filter(is.na(ICV)==F)
write.csv(data, "data/neuroquant.csv",row.names = F)




