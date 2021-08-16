library(tidyverse)
library(lubridate)
library(pROC)
library(slider)
library(fields)

setwd("yourpath")


###Bangladesh
#Load Bangladesh App data and AFE data
AFE <- read_csv("yourpath/dlbg with AFEs.csv")
AppDat <- read_csv("yourpath/DiaPRData2020.01.01.csv")



#Plot to show vomiting issue
qplot(AdmitDate,jitter(as.numeric(`Form 2::Vomit1`)),data=AppDat %>% select(AdmitDate,`Form 2::Vomit1`) %>% mutate(`Form 2::Vomit1`=case_when(`Form 2::Vomit1`=="Yes"~1,TRUE~0))) +
  ylab("Jittered Vomiting") + xlab("Date")


#Define etiology using AFe threshold 0.5
AFE_0.5=AFE %>% select(ID,ends_with("_afe")) %>% mutate_at(vars(ends_with("_afe")),~.>.5) 


Type=c("viral","bacterial","viral","bacterial","parasite","parasite","parasite","parasite","viral",
       "viral","bacterial","viral","bacterial","bacterial","bacterial","bacterial","bacterial","bacterial")

#Check that columns are assigned properly
cbind(colnames(AFE_0.5)[-1],Type)

#Assign etiology to IDs
AFE_0.5$viral=apply(AFE_0.5[,-1],1, function(x) any(Type[which(x==T)]=="viral"))
AFE_0.5$bacterial=apply(AFE_0.5[,-1],1, function(x) any(Type[which(x==T)]=="bacterial"))
AFE_0.5$parasite=apply(AFE_0.5[,-1],1, function(x) any(Type[which(x==T)]=="parasite"))

#Define viral only and any bacteria (just in case) 
AFE_0.5 = AFE_0.5 %>%  replace(is.na(.), 0) %>% mutate(viral_only=case_when(viral==T & bacterial==F & parasite==F ~ 1, TRUE ~0),any_bac=case_when(bacterial==T~1,TRUE~0))


# select columns of interest and combine partial and exclusive breastfeeding
dat=AppDat %>% select(ID=`Study ID`,AdmitDate,AgeMonths,contains("bloodstoolreported",ignore.case=F),contains('vomit',ignore.case=F),
                      contains('breastfeeding',ignore.case = F),contains("MUAC"),Gender,
                      DiarrheaDays,
                      DiarrheaHours,
                      DiarrheaEpisodes,
                      `Form 2::fever1.1`,
                      `Form 4::Medications`,
                      `Form 4::MotherEducation`,
                      `Form 4::FatherEducation`,
                      `Form 4::PeopleAtHome`) %>%
  mutate_at(vars(contains('breastfeed')),~case_when(. %in% c(1,2)~1,TRUE~0))

# combine distinct clinical data and AFe data
dat=AFE_0.5 %>% distinct %>% inner_join(dat %>% distinct,by=c("ID")) 

#get pathogen frequences and combinations
freq_paths=t(dat %>% select(ends_with("afe")) %>% summarize_all(~sum(.)))
path_names=colnames(dat %>% select(ends_with("afe")))
path_names=substr(path_names,1,nchar(path_names)-4)
path_combs=table(unlist(apply(dat %>% select(ends_with("afe")),1,function(x) paste(path_names[which(x)],collapse="+"))))
df_paths=data.frame(Pathogen=as.character(names(path_combs)),Frequency=as.numeric(path_combs))
df_paths$Pathogen=as.character(df_paths$Pathogen)
df_paths[1,1]="None"
bang_paths=df_paths %>% arrange(desc(Frequency))
#write.csv(bang_paths,file="bang_paths.csv")

# Clean data so it looks like GEMS prediction model 
dat$f4b_muac=apply(dat %>% select(contains("MUAC")),1,mean) # average MUAC measurements 
dat=dat %>% mutate(f4a_drh_blood=case_when(`Form 2::bloodstoolreported1.1`==1&`Form 3::bloodstoolreported2.1`==1~1,TRUE~0)) %>% 
  select(ID,viral,bacterial,parasite,viral_only,AdmitDate,base_age=AgeMonths,f4a_drh_blood,f4a_drh_vomit=`Form 2::vomit1.1`,f4b_muac,any_breast_fed=`Form 2::breastfeeding1.1`,
         Gender,
         DiarrheaDays,
         DiarrheaHours,
         DiarrheaEpisodes,
         `Form 2::fever1.1`,
         `Form 4::Medications`,
         `Form 4::MotherEducation`,
         `Form 4::FatherEducation`,
         `Form 4::PeopleAtHome`) %>%
  mutate_at(c("f4a_drh_blood","f4a_drh_vomit","any_breast_fed"),as.factor)


# Put date on numerical scale using GEMS origin date
dat$dt=as.numeric(mdy(dat$AdmitDate))-13848

#Use age as in GEMS
dat=dat %>% mutate(base_age = case_when(base_age <= 24 ~ base_age,
                                        base_age>24 & base_age < 36 ~ as.double(24),
                                        base_age>=36 & base_age < 48 ~ as.double(36),
                                        base_age>=48 & base_age < 60 ~ as.double(48)),
                   Seasonal_sine=sin(2*pi*.$dt/365.25),
                   Seasonal_cosine=cos(2*pi*.$dt/365.25))

#Gets rid of unknown etiology cases
dat = dat %>% filter(!(viral==0 & bacterial==0 & parasite==0)) #Do we want to do this? 


#write.csv(dat_avg,"dhaka_weather.csv")
# load weather over study period and combine with main data 
dat_avg=read.csv("dhaka_weather.csv")
#dat=dat %>% left_join(dat_avg,by="dt")
clim_band=14
temp.aggr=dat_avg %>% group_by(dt) %>% summarize(temp=mean(temp)) %>% tidyr::complete(dt=seq(min(dt),max(dt))) %>% mutate(temp.aggr=slide_index_dbl(.x=temp,.i=dt,.f=~mean(.x,na.rm=T),.before=clim_band,.after=-1))
rain.aggr=dat_avg %>% group_by(dt) %>% summarize(rain=mean(rain)) %>% tidyr::complete(dt=seq(min(dt),max(dt))) %>% mutate(rain.aggr=slide_index_dbl(.x=rain,.i=dt,.f=~mean(.x,na.rm=T),.before=clim_band,.after=-1))
dat=dat %>% left_join(temp.aggr,by="dt") %>% left_join(rain.aggr,by="dt")

#If we want to exclude when vomiting was incorrectly collected
dat_bangladesh=dat %>% filter(dt>=4385)

### Mali
#Load Mali clinical and AFe data
AFE <- read_csv("yourpath/Mali_Afe.csv")
AppDat <- read_csv("yourpath/base_smartphone.csv")


#Mali weather
weather=read_csv("Mali_weather.csv")
#Rename from French to English and take average of MUAC
dat=AppDat %>% select(ID=Code_Patient,Date_Inclusion,base_age=`Age-mois`,f4a_drh_vomit=Vomissement1,MUAC1,MUAC2,any_breast_fed=Allaitement1,f4a_drh_blood=`Sang-Selle1`,
                      Gender=Sexe,
                      DiarrheaDays=Jours_Diarrhee,
                      DiarrheaHours=Heure_Diarrhee,
                      DiarrheaEpisodes=Nombre_Selle,
                      `Form 2::fever1.1`=Fievre1,
                      `Form 4::Medications`=Medicaments,
                      `Form 4::MotherEducation`=Educ_Mere,
                      `Form 4::FatherEducation`=Educ_Pere,
                      `Form 4::PeopleAtHome`=Pers_menage,
                      `Medic-tradi`)
dat$f4b_muac=apply(dat%>% select(contains("MUAC")),1,mean)
#Set AFe threshold of 0.5
AFE_0.5=AFE %>% select(ID,ends_with("_afe")) %>% mutate_at(vars(ends_with("_afe")),~.>.5) 


Type=c("viral","bacterial","viral","bacterial","parasite","parasite","parasite","parasite","viral",
       "viral","bacterial","viral","bacterial","bacterial","bacterial","bacterial","bacterial","bacterial")

cbind(colnames(AFE_0.5)[-1],Type)



AFE_0.5$viral=apply(AFE_0.5[,-1],1, function(x) any(Type[which(x==T)]=="viral"))
AFE_0.5$bacterial=apply(AFE_0.5[,-1],1, function(x) any(Type[which(x==T)]=="bacterial"))
AFE_0.5$parasite=apply(AFE_0.5[,-1],1, function(x) any(Type[which(x==T)]=="parasite"))
# 

#Define viral only etiology
AFE_0.5 = AFE_0.5 %>%  replace(is.na(.), 0) %>% mutate(viral_only=case_when(viral==T & bacterial==F & parasite==F ~ 1, TRUE ~0),any_bac=case_when(bacterial==T~1,TRUE~0),
                                                       No_known=case_when(viral==F & bacterial==F & parasite==F ~ 1, TRUE ~0))


dat=AFE_0.5 %>% inner_join(dat,by=c("ID")) 

#Pathogen combination and Frequency
freq_paths=t(dat %>% select(ends_with("afe")) %>% summarize_all(~sum(.)))
path_names=colnames(dat %>% select(ends_with("afe")))
path_names=substr(path_names,1,nchar(path_names)-4)
path_combs=table(unlist(apply(dat %>% select(ends_with("afe")),1,function(x) paste(path_names[which(x)],collapse="+"))))
df_paths=data.frame(Pathogen=as.character(names(path_combs)),Frequency=as.numeric(path_combs))
df_paths$Pathogen=as.character(df_paths$Pathogen)
df_paths[1,1]="None"
mali_paths=df_paths %>% arrange(desc(Frequency))

#write.csv(mali_paths,"mali_paths.csv")

#dat=dat %>% select(-ends_with("afe"))

dat$Date_Inclusion[109]="2/18/2020" # fixing this date

dat = dat %>% mutate(date=case_when(substr(Date_Inclusion,2,2)=="/"~mdy(Date_Inclusion),TRUE~dmy(Date_Inclusion)))
dat = dat %>% arrange(date)


#Put dates on GEMS scale
dat$dt=as.numeric(dat$date)-13848


dat=dat %>% mutate(base_age = case_when(base_age <= 24 ~ base_age,
                                        base_age>24 & base_age < 36 ~ as.double(24),
                                        base_age>=36 & base_age < 48 ~ as.double(36),
                                        base_age>=48 & base_age < 60 ~ as.double(48)),
                   Seasonal_sine=sin(2*pi*.$dt/365.25),
                   Seasonal_cosine=cos(2*pi*.$dt/365.25))

dat=dat %>% mutate_at(vars(c("f4a_drh_vomit","any_breast_fed","f4a_drh_blood")),~as.factor(ifelse(.=="Oui",1,0)))

#Combine Mali weather with clinical data
clim_band=14
weather$dt=as.numeric(weather$date)-13848


temp.aggr=weather %>% group_by(dt) %>% summarize(temp=mean(temp)) %>% tidyr::complete(dt=seq(min(dt),max(dt))) %>% mutate(temp.aggr=slide_index_dbl(.x=temp,.i=dt,.f=~mean(.x,na.rm=T),.before=clim_band,.after=-1))
rain.aggr=weather %>% group_by(dt) %>% summarize(rain=mean(rain)) %>% tidyr::complete(dt=seq(min(dt),max(dt))) %>% mutate(rain.aggr=slide_index_dbl(.x=rain,.i=dt,.f=~mean(.x,na.rm=T),.before=clim_band,.after=-1))
dat=dat %>% left_join(temp.aggr,by="dt") %>% left_join(rain.aggr,by="dt")

dat_mali = dat %>% filter(!(viral==0 & bacterial==0 & parasite==0)) 
###

dat_bangladesh=dat_bangladesh %>% rename(date=AdmitDate) %>% mutate(Site=6,date=mdy(date))
dat_mali=dat_mali %>% select(one_of(colnames(dat_bangladesh))) %>% mutate(Site=2)


#Combine Mali and bangladesh 
dat_combo=rbind(dat_bangladesh,dat_mali)



# Analysis using the prediction models from the App

out=readRDS("yourpath//CurPatMod")

#Current patient (CP) Only
dat_combo$pred=predict(out,newdata=dat_combo,type="response")

xs=readRDS("yourpath//xs")
OR=readRDS("yourpath//CurPatOR")
closest=predict(out,newdata=dat_combo,type="response") %>% purrr::map(function(y) vapply(y,function(x) which.min(abs(xs-x)),c(1)))
LR=OR[unlist(closest)]    
dat_combo$LR=LR
CurPatAUC=with(dat_combo,roc(viral_only,LR,plot=T,ci=T));paste0(round(CurPatAUC$auc,3)," (",round(CurPatAUC$ci[1],3),"-",round(CurPatAUC$ci[3],3),")")
dat_combo$LRpred=with(dat_combo,plogis(LR))

dat_combo$LRobs_prop=as.numeric(unlist(dat_combo$LRpred %>% purrr::map(function(z)  
  dat_combo %>% filter(LRpred>z-.05,LRpred<z+.05) %>% summarize(mean(viral_only)))))


CurPatCalib=summary(lm(dat_combo$LRobs_prop~dat_combo$LRpred))$coef[1:2,1]
confint(lm(dat_combo$LRobs_prop~dat_combo$LRpred))


#Seasonal + CP 
mods2=readRDS("yourpath//SeasMod")
OR_Seas=readRDS("yourpath//OR_Seas")
closest_seas=dat_combo %>% split(.$Site) %>% purrr::map(~predict(mods2[[.$Site[1]]],newdata=.,type="response")) %>% purrr::map(function(y) vapply(y,function(x) which.min(abs(xs-x)),c(1)))
num=dat_combo %>% split(.$Site) %>% purrr::map(~.$Site[1])
LR_Seas=map2(closest_seas,num, function(x,y) OR_Seas[[y]][as.vector(unlist(x))])
dat_combo$LR_Seas=c(LR_Seas[['6']],LR_Seas[['2']])
SeasAUC=with(dat_combo,roc(viral_only,LR+LR_Seas,plot=T,ci=T))

#
dt=as.numeric(ymd(today()))-13848

#
dat_combo$Spred=with(dat_combo,plogis(LR+LR_Seas))

dat_combo$Sobs_prop=as.numeric(unlist(dat_combo$Spred %>% purrr::map(function(z)  
  dat_combo %>% filter(Spred>z-.05,Spred<z+.05) %>% summarize(mean(viral_only)))))

SeasCalib=summary(lm(unlist(obs_prop)~dat_combo$Spred))$coef[1:2,1]

confint(lm(unlist(obs_prop)~dat_combo$Spred))



#Weather + CP
mods3=readRDS("C://Users//u6020766//Google Drive//UU//VIDA//PTOmods//ClimMods")
OR_Clim=readRDS("C://Users//u6020766//Google Drive//UU//VIDA//PTOmods//OR_Clim")
closest_clim=dat_combo %>% split(.$Site) %>% purrr::map(~predict(mods3[[.$Site[1]]],newdata=.,type="response")) %>% purrr::map(function(y) vapply(y,function(x) which.min(abs(xs-x)),c(1)))
LR_clim=map2(closest_clim,num, function(x,y) OR_Clim[[y]][as.vector(unlist(x))])
dat_combo$LR_clim=c(LR_clim[['6']],LR_clim[['2']])
ClimAUC=with(dat_combo,roc(viral_only,LR+LR_clim,plot=T,ci=T))

with(dat_combo %>% filter(Site==6),roc(viral_only,LR+LR_clim,plot=T,ci=T))

dat_combo$Cpred=with(dat_combo,plogis(LR+LR_clim))

dat_combo$Climobs_prop=as.numeric(unlist(dat_combo$Cpred %>% purrr::map(function(z)  
  dat_combo %>% filter(Cpred>z-.05,Cpred<z+.05) %>% summarize(mean(viral_only)))))

ClimCalib=summary(lm(unlist(obs_prop)~dat_combo$Cpred))$coef[1:2,1]

confint(lm(unlist(obs_prop)~dat_combo$Cpred))

#Historical Patient + CP
priordf=read_csv("C:/Users/u6020766/Google Drive/UU/GEMS/gems_prior_v2.csv")
priordf

dat_combo=dat_combo %>% mutate(daymo=(as.numeric(date) %% 365)) %>% arrange(daymo)
dat_combo=dat_combo %>% left_join(priordf,by=c("daymo","Site"="site"))
PrAUC=with(dat_combo,roc(viral_only,LR+preOdds,plot=T,ci=T))
dat_combo$PrPred=plogis(dat_combo$LR+dat_combo$preOdds)

dat_combo$Probs_prop=as.numeric(unlist(dat_combo$PrPred %>% purrr::map(function(z)  
  dat_combo %>% filter(PrPred>z-.05,PrPred<z+.05) %>% summarize(mean(viral_only)))))

PrCalib=summary(lm(unlist(obs_prop)~dat_combo$PrPred))$coef[1:2,1]
confint(lm(unlist(obs_prop)~dat_combo$PrPred))


#Recent Patient + CP 

weighted_avg_fun=function(x,y) {
  #wts=rep(1,clin_band)
  #wts=1/(1:clin_band)
  wts=Wendland(1:clin_band,theta=clin_band,k=1,dimension=1)
  if(length(y)==1) {NA} else {
    today=y[length(y)]
    inds=today-y
    x=x[-length(x)]
    weighted.mean(x,w=wts[inds],na.rm=T)}
}

clin_band=28
test_dat=dat_combo %>% group_by(dt,Site) %>% summarize(mean_pred=mean(LRpred,na.rm=T))
test_dat
test_split_wAvg=test_dat %>% split(.$Site) %>% purrr::map(. %>% mutate(lens=slide_index(.x=dt,.i=dt,~.x,.before=clin_band,.after=-1,.names_to="tst")) %>% unnest(lens) %>% 
  mutate(dts=dt-lens) %>% dplyr::select(dt,dts) %>% group_by(dt) %>% nest() %>% right_join(test_dat,by="dt") %>% ungroup() %>% 
  mutate(roll_mean=slide_index2(.x=mean_pred,.y=dt,.i=dt,~weighted_avg_fun(.x,.y),.before=clin_band,.after=0)) %>% unnest(roll_mean) %>% replace(is.na(.), .5))

test_preOdds=bind_rows(test_split_wAvg) %>% mutate(preOdds_v2=log(roll_mean/(1-roll_mean)))

dat_combo=dat_combo  %>% left_join(test_preOdds,by=c("dt","Site")) 

RecPatAUC=with(dat_combo,roc(viral_only,LR+preOdds_v2,plot=T,ci=T))
dat_combo$RecPatPred=plogis(dat_combo$LR+dat_combo$preOdds_v2)


dat_combo$RecPobs_prop=as.numeric(unlist(dat_combo$RecPatPred %>% purrr::map(function(z)  
  dat_combo %>% filter(RecPatPred>z-.05,RecPatPred<z+.05) %>% summarize(mean(viral_only)))))

qplot(dat_combo$RecPatPred,unlist(obs_prop)) + geom_line(aes(x=c(0,1),y=c(0,1))) + scale_x_continuous(breaks=seq(0,1,by=.05),labels=seq(0,1,by=.05))
RecPatCalib=summary(lm(unlist(obs_prop)~dat_combo$RecPatPred))$coef[1:2,1]
confint(lm(unlist(obs_prop)~dat_combo$RecPatPred))



result=data.frame(row.names = c("Current Patient Only","Seasonal","Climate","Historical Patient","Recent Patient"),
           A=c(CurPatCalib[1],SeasCalib[1],ClimCalib[1],PrCalib[1],RecPatCalib[1]),
           B=c(CurPatCalib[2],SeasCalib[2],ClimCalib[2],PrCalib[2],RecPatCalib[2]),
           AUC=c(paste0(round(CurPatAUC$auc,3)," (",round(CurPatAUC$ci[1],3),"-",round(CurPatAUC$ci[3],3),")"),
                 paste0(round(SeasAUC$auc,3)," (",round(SeasAUC$ci[1],3),"-",round(SeasAUC$ci[3],3),")"),
                 paste0(round(ClimAUC$auc,3)," (",round(ClimAUC$ci[1],3),"-",round(ClimAUC$ci[3],3),")"),
                 paste0(round(PrAUC$auc,3)," (",round(PrAUC$ci[1],3),"-",round(PrAUC$ci[3],3),")"),
                 paste0(round(RecPatAUC$auc,3)," (",round(RecPatAUC$ci[1],3),"-",round(RecPatAUC$ci[3],3),")"))) %>% 
  select(AUC,A,B)

colnames(dat_combo)

a=ggplot(dat_combo,aes(LRpred,LRobs_prop)) + geom_point() + geom_abline(intercept=0,slope=1) + 
  scale_x_continuous(limit=c(0,1),breaks=seq(0,1,by=.25),labels=seq(0,1,by=.25)) + scale_y_continuous(limit=c(0,1)) + theme_bw() +
  ylab("Observed Proportion") + xlab("Current Patient Only Prediction")

b=ggplot(dat_combo,aes(Spred,Sobs_prop)) + geom_point() + geom_abline(intercept=0,slope=1) + 
  scale_x_continuous(limit=c(0,1),breaks=seq(0,1,by=.25),labels=seq(0,1,by=.25)) + scale_y_continuous(limit=c(0,1)) + theme_bw() +
  ylab("Observed Proportion") + xlab("Seasonal Prediction")

c=ggplot(dat_combo,aes(Cpred,Climobs_prop)) + geom_point() + geom_abline(intercept=0,slope=1) + 
  scale_x_continuous(limit=c(0,1),breaks=seq(0,1,by=.25),labels=seq(0,1,by=.25)) + scale_y_continuous(limit=c(0,1)) + theme_bw() +
  ylab("Observed Proportion") + xlab("Climate Prediction")

d=ggplot(dat_combo,aes(PrPred,Probs_prop)) + geom_point() + geom_abline(intercept=0,slope=1) + 
  scale_x_continuous(limit=c(0,1),breaks=seq(0,1,by=.25),labels=seq(0,1,by=.25)) + scale_y_continuous(limit=c(0,1)) + theme_bw() +
  ylab("Observed Proportion") + xlab("Historical Patient Prediction")

e=ggplot(dat_combo,aes(RecPatPred,RecPobs_prop)) + geom_point() + geom_abline(intercept=0,slope=1) + 
  scale_x_continuous(limit=c(0,1),breaks=seq(0,1,by=.25),labels=seq(0,1,by=.25)) + scale_y_continuous(limit=c(0,1)) + theme_bw() +
  ylab("Observed Proportion") + xlab("Recent Patient Prediction")


a=ggplot(dat_combo,aes(LRpred,viral_only)) + geom_point() + geom_smooth(method="loess") + geom_abline(intercept=0,slope=1) + 
  scale_x_continuous(limit=c(0,1),breaks=seq(0,1,by=.25),labels=seq(0,1,by=.25)) + scale_y_continuous(limit=c(0,1)) + theme_bw() +
  ylab("Observed Proportion") + xlab("Current Patient Only Prediction")

b=ggplot(dat_combo,aes(Spred,viral_only)) + geom_point() + geom_smooth(method="loess") + geom_abline(intercept=0,slope=1) + 
  scale_x_continuous(limit=c(0,1),breaks=seq(0,1,by=.25),labels=seq(0,1,by=.25)) + scale_y_continuous(limit=c(0,1)) + theme_bw() +
  ylab("Observed Proportion") + xlab("Seasonal Prediction")

c=ggplot(dat_combo,aes(Cpred,viral_only))+ geom_point() + geom_smooth(method="loess") + geom_abline(intercept=0,slope=1) + 
  scale_x_continuous(limit=c(0,1),breaks=seq(0,1,by=.25),labels=seq(0,1,by=.25)) + scale_y_continuous(limit=c(0,1)) + theme_bw() +
  ylab("Observed Proportion") + xlab("Climate Prediction")

d=ggplot(dat_combo,aes(PrPred,viral_only)) + geom_point() + geom_smooth(method="loess") + geom_abline(intercept=0,slope=1) + 
  scale_x_continuous(limit=c(0,1),breaks=seq(0,1,by=.25),labels=seq(0,1,by=.25)) + scale_y_continuous(limit=c(0,1)) + theme_bw() +
  ylab("Observed Proportion") + xlab("Historical Patient Prediction")

e=ggplot(dat_combo,aes(RecPatPred,viral_only)) + geom_point() + geom_smooth(method="loess") + geom_abline(intercept=0,slope=1) + 
  scale_x_continuous(limit=c(0,1),breaks=seq(0,1,by=.25),labels=seq(0,1,by=.25)) + scale_y_continuous(limit=c(0,1)) + theme_bw() +
  ylab("Observed Proportion") + xlab("Recent Patient Prediction")

grid.arrange(a,b,c,d,e,ncol=3)

write.csv(result,file="overall_res.csv")


#Results split by site 
site_auc=dat_combo %>% split(.$Site) %>% purrr::map(function(x) x %>% summarize(
                                                           roc(viral_only,LR,ci=T)$auc,
                                                           roc(viral_only,LR+LR_Seas,ci=T)$auc,
                                                           roc(viral_only,LR+LR_clim,ci=T)$auc,
                                                           roc(viral_only,LR+preOdds,ci=T)$auc,
                                                           roc(viral_only,LR+preOdds_v2,ci=T)$auc))
site_auc_lo=dat_combo %>% split(.$Site) %>% purrr::map(function(x) x %>% summarize(
                                                          roc(viral_only,LR,ci=T)$ci[1],
                                                           roc(viral_only,LR+LR_Seas,ci=T)$ci[1],
                                                           roc(viral_only,LR+LR_clim,ci=T)$ci[1],
                                                           roc(viral_only,LR+preOdds,ci=T)$ci[1],
                                                           roc(viral_only,LR+preOdds_v2,ci=T)$ci[1]))
site_auc_hi=dat_combo %>% split(.$Site) %>% purrr::map(function(x) x %>% summarize(
                                                           roc(viral_only,LR,ci=T)$ci[3],
                                                           roc(viral_only,LR+LR_Seas,ci=T)$ci[3],
                                                           roc(viral_only,LR+LR_clim,ci=T)$ci[3],
                                                           roc(viral_only,LR+preOdds,ci=T)$ci[3],
                                                           roc(viral_only,LR+preOdds_v2,ci=T)$ci[3]))

site_AUCs=data.frame(row.names = c("Current Patient Only","Seasonal","Climate","Historical Patient","Recent Patient"),
           matrix(paste0(round(unlist(site_auc),3)," (",
             round(unlist(site_auc_lo),3)," - ",
                   round(unlist(site_auc_hi),3),")"),ncol=2)) %>% rename("Mali"=X1,"Bangladesh"=X2)
write.csv(site_AUCs,file="site_AUCs.csv")



data.frame(Prediction=c(plogis(LR),plogis(LR2)),Etiology=c(rep("Known Etiology",126),rep("Unknown Etiology",89))) %>% ggplot(aes(x=Prediction)) + xlab("Current Patient Prediction") + 
  geom_density() + facet_wrap(~Etiology) + theme_bw()

grid.arrange(qplot(LR),qplot(LR2))

##Descriptive summary using dat_combo
dat_combo$`Form 4::MotherEducation`
colnames(dat_combo)
DS=t(dat_combo %>% summarize(n=n(),
                        Age=paste0(median(base_age)," (",IQR(base_age),")"),
                        Blank=NA,
                        M=paste0(sum(Gender=="Male"|Gender=="M")," (",100*round(sum(Gender=="Male"|Gender=="M")/n,3),"%)"),
                        F=paste0(sum(Gender=="Female"|Gender=="F")," (",100*round(sum(Gender=="Female"|Gender=="F")/n,3),"%)"),
                        DD=paste0(median(DiarrheaDays+DiarrheaHours/24,na.rm=T)," (",round(IQR(DiarrheaDays+DiarrheaHours/24,na.rm=T),3),")"),
                        Epi=paste0(median(DiarrheaEpisodes)," (",IQR(DiarrheaEpisodes),")"),
                        Blank2=NA,
                        B1=paste0(sum(f4a_drh_blood==1)," (",100*round(sum(f4a_drh_blood==1)/n,3),"%)"),
                        B0=paste0(sum(f4a_drh_blood==0)," (",100*round(sum(f4a_drh_blood==0)/n,3),"%)"),
                        Blank3=NA,
                        F1=paste0(sum(`Form 2::fever1.1`=="1"|`Form 2::fever1.1`=="Oui")," (",100*round(sum(`Form 2::fever1.1`=="1"|`Form 2::fever1.1`=="Oui")/n,3),"%)"),
                        F0=paste0(sum(`Form 2::fever1.1`=="0"|`Form 2::fever1.1`=="Non")," (",100*round(sum(`Form 2::fever1.1`=="0"|`Form 2::fever1.1`=="Non")/n,3),"%)"),
                        Blank4=NA,
                        V1=paste0(sum(f4a_drh_vomit==1)," (",100*round(sum(f4a_drh_vomit==1)/n,3),"%)"),
                        V1b=paste0(sum(f4a_drh_vomit==1)," (",100*round(sum(f4a_drh_vomit==1)/n,3),"%)"),
                        V0=paste0(sum(f4a_drh_vomit==0)," (",100*round(sum(f4a_drh_vomit==0)/n,3),"%)"),
                        Blank5=NA,
                        BF1=paste0(sum(any_breast_fed==1)," (",100*round(sum(any_breast_fed==1)/n,3),"%)"),
                        BF0=paste0(sum(any_breast_fed==0)," (",100*round(sum(any_breast_fed==0)/n,3),"%)"),
                        MUAC=paste0(median(f4b_muac)," (",IQR(f4b_muac),")"),
                        Blank6=NA,
                        PM1=paste0(sum(`Form 4::Medications`=="Yes"|`Form 4::Medications`=="Oui")," (",100*round(sum(`Form 4::Medications`=="Yes"|`Form 4::Medications`=="Oui")/n,3),"%)"),
                        PM2=paste0(sum(`Form 4::Medications`=="Non"|`Form 4::Medications`=="No")," (",100*round(sum(`Form 4::Medications`=="No"|`Form 4::Medications`=="Non")/n,3),"%)"),
                        Blank7=NA,
                        ME=paste0(median(`Form 4::MotherEducation`)," (",IQR(`Form 4::MotherEducation`),")"),
                        FE=paste0(median(`Form 4::FatherEducation`)," (",IQR(`Form 4::MotherEducation`),")"),
                        PLH=paste0(median(`Form 4::PeopleAtHome`)," (",IQR(`Form 4::PeopleAtHome`),")")
                        
))

DS_B=t(dat_combo %>% filter(Site==6) %>% summarize(n=n(),
                                Age=paste0(median(base_age)," (",IQR(base_age),")"),
                                Blank=NA,
                                M=paste0(sum(Gender=="Male"|Gender=="M")," (",100*round(sum(Gender=="Male"|Gender=="M")/n,3),"%)"),
                                F=paste0(sum(Gender=="Female"|Gender=="F")," (",100*round(sum(Gender=="Female"|Gender=="F")/n,3),"%)"),
                                DD=paste0(median(DiarrheaDays+DiarrheaHours/24,na.rm=T)," (",round(IQR(DiarrheaDays+DiarrheaHours/24,na.rm=T),3),")"),
                                Epi=paste0(median(DiarrheaEpisodes)," (",IQR(DiarrheaEpisodes),")"),
                                Blank2=NA,
                                B1=paste0(sum(f4a_drh_blood==1)," (",100*round(sum(f4a_drh_blood==1)/n,3),"%)"),
                                B0=paste0(sum(f4a_drh_blood==0)," (",100*round(sum(f4a_drh_blood==0)/n,3),"%)"),
                                Blank3=NA,
                                F1=paste0(sum(`Form 2::fever1.1`=="1"|`Form 2::fever1.1`=="Oui")," (",100*round(sum(`Form 2::fever1.1`=="1"|`Form 2::fever1.1`=="Oui")/n,3),"%)"),
                                F0=paste0(sum(`Form 2::fever1.1`=="0"|`Form 2::fever1.1`=="Non")," (",100*round(sum(`Form 2::fever1.1`=="0"|`Form 2::fever1.1`=="Non")/n,3),"%)"),
                                Blank4=NA,
                                V1=paste0(sum(f4a_drh_vomit==1)," (",100*round(sum(f4a_drh_vomit==1)/n,3),"%)"),
                                V1b=paste0(sum(f4a_drh_vomit==1)," (",100*round(sum(f4a_drh_vomit==1)/n,3),"%)"),
                                V0=paste0(sum(f4a_drh_vomit==0)," (",100*round(sum(f4a_drh_vomit==0)/n,3),"%)"),
                                Blank5=NA,
                                BF1=paste0(sum(any_breast_fed==1)," (",100*round(sum(any_breast_fed==1)/n,3),"%)"),
                                BF0=paste0(sum(any_breast_fed==0)," (",100*round(sum(any_breast_fed==0)/n,3),"%)"),
                                MUAC=paste0(median(f4b_muac)," (",IQR(f4b_muac),")"),
                                Blank6=NA,
                                PM1=paste0(sum(`Form 4::Medications`=="Yes"|`Form 4::Medications`=="Oui")," (",100*round(sum(`Form 4::Medications`=="Yes"|`Form 4::Medications`=="Oui")/n,3),"%)"),
                                PM2=paste0(sum(`Form 4::Medications`=="Non"|`Form 4::Medications`=="No")," (",100*round(sum(`Form 4::Medications`=="No"|`Form 4::Medications`=="Non")/n,3),"%)"),
                                Blank7=NA,
                                ME=paste0(median(`Form 4::MotherEducation`)," (",IQR(`Form 4::MotherEducation`),")"),
                                FE=paste0(median(`Form 4::FatherEducation`)," (",IQR(`Form 4::MotherEducation`),")"),
                                PLH=paste0(median(`Form 4::PeopleAtHome`)," (",IQR(`Form 4::PeopleAtHome`),")")
                                
))
DS_M=t(dat_combo %>% filter(Site==2) %>% summarize(n=n(),
                                                   Age=paste0(median(base_age)," (",IQR(base_age),")"),
                                                   Blank=NA,
                                                   M=paste0(sum(Gender=="Male"|Gender=="M")," (",100*round(sum(Gender=="Male"|Gender=="M")/n,3),"%)"),
                                                   F=paste0(sum(Gender=="Female"|Gender=="F")," (",100*round(sum(Gender=="Female"|Gender=="F")/n,3),"%)"),
                                                   DD=paste0(median(DiarrheaDays+DiarrheaHours/24,na.rm=T)," (",round(IQR(DiarrheaDays+DiarrheaHours/24,na.rm=T),3),")"),
                                                   Epi=paste0(median(DiarrheaEpisodes)," (",IQR(DiarrheaEpisodes),")"),
                                                   Blank2=NA,
                                                   B1=paste0(sum(f4a_drh_blood==1)," (",100*round(sum(f4a_drh_blood==1)/n,3),"%)"),
                                                   B0=paste0(sum(f4a_drh_blood==0)," (",100*round(sum(f4a_drh_blood==0)/n,3),"%)"),
                                                   Blank3=NA,
                                                   F1=paste0(sum(`Form 2::fever1.1`=="1"|`Form 2::fever1.1`=="Oui")," (",100*round(sum(`Form 2::fever1.1`=="1"|`Form 2::fever1.1`=="Oui")/n,3),"%)"),
                                                   F0=paste0(sum(`Form 2::fever1.1`=="0"|`Form 2::fever1.1`=="Non")," (",100*round(sum(`Form 2::fever1.1`=="0"|`Form 2::fever1.1`=="Non")/n,3),"%)"),
                                                   Blank4=NA,
                                                   V1=paste0(sum(f4a_drh_vomit==1)," (",100*round(sum(f4a_drh_vomit==1)/n,3),"%)"),
                                                   V1b=paste0(sum(f4a_drh_vomit==1)," (",100*round(sum(f4a_drh_vomit==1)/n,3),"%)"),
                                                   V0=paste0(sum(f4a_drh_vomit==0)," (",100*round(sum(f4a_drh_vomit==0)/n,3),"%)"),
                                                   Blank5=NA,
                                                   BF1=paste0(sum(any_breast_fed==1)," (",100*round(sum(any_breast_fed==1)/n,3),"%)"),
                                                   BF0=paste0(sum(any_breast_fed==0)," (",100*round(sum(any_breast_fed==0)/n,3),"%)"),
                                                   MUAC=paste0(median(f4b_muac)," (",IQR(f4b_muac),")"),
                                                   Blank6=NA,
                                                   PM1=paste0(sum(`Form 4::Medications`=="Yes"|`Form 4::Medications`=="Oui")," (",100*round(sum(`Form 4::Medications`=="Yes"|`Form 4::Medications`=="Oui")/n,3),"%)"),
                                                   PM2=paste0(sum(`Form 4::Medications`=="Non"|`Form 4::Medications`=="No")," (",100*round(sum(`Form 4::Medications`=="No"|`Form 4::Medications`=="Non")/n,3),"%)"),
                                                   Blank7=NA,
                                                   ME=paste0(median(`Form 4::MotherEducation`)," (",IQR(`Form 4::MotherEducation`),")"),
                                                   FE=paste0(median(`Form 4::FatherEducation`)," (",IQR(`Form 4::MotherEducation`),")"),
                                                   PLH=paste0(median(`Form 4::PeopleAtHome`)," (",IQR(`Form 4::PeopleAtHome`),")")
                                                   
))


write.csv(cbind(DS,DS_B,DS_M),"DS.csv")

