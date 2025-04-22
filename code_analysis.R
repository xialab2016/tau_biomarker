library(ggplot2)
library(patchwork)
library(tableone)
library(GGally)
library(dplyr)
library(nnet)
library(ggplot2)
library(patchwork)
library(MatchIt)
library(pROC)
library(lme4)
library(gee)
library(DescTools)
library(caret)
library(dunn.test)
library(geepack)

pls=read.csv("~/OneDrive - University of Pittsburgh/XiaLab/ADBiomarker/Data_Archive/plasma_CleanData_230828.csv")
pls=pls %>% filter(subtype_blood !="Related disorder")
pls=pls %>% filter(is.na(pTau217)==F)
#pls=pls %>% filter(id_participant !="PRT170468")
#### check for outlier ####
summary(pls$GFAP)
hist(pls$GFAP)
pls %>% arrange(GFAP) %>%mutate(id_participant=factor(id_participant,levels = id_participant))%>%
  ggplot( aes(x=id_participant,y=GFAP))+
  geom_point(shape=21)+
  geom_segment(aes(xend=id_participant, yend=0),
                alpha=0.5)
pls %>% arrange(NFL) %>%mutate(id_participant=factor(id_participant,levels = id_participant))%>%
  ggplot( aes(x=id_participant,y=NFL))+
  geom_point(shape=21)+
  geom_segment(aes(xend=id_participant, yend=0),
               alpha=0.5)

#### Plasma marker in MS and controls (MS vs. control) ####
pls=pls %>% mutate(race_ethnicity=ifelse(race=="White" & ethnicity=="Not Hispanic or Latino","Non-Hispanic White","Other" ))
pls=pls %>% mutate(subtype_blood1=ifelse(subtype_blood %in% c("RIS","RRMS"),"RRMS",
                                         ifelse(subtype_blood=="SPMS","SPMS",
                                                ifelse(subtype_blood=="PPMS","PPMS","Control"))))
pls$subtype_blood1=factor(pls$subtype_blood1, levels = c("Control","RRMS","PPMS","SPMS"))
table(pls$subtype_blood1) 
#Control    RRMS    PPMS    SPMS 
#20     117      24      19 
pls$subtype_blood2=as.character(pls$subtype_blood1)
pls$subtype_blood2[pls$subtype_blood2 %in% c("PPMS","SPMS")]="PMS"
pls$subtype_blood2=factor(pls$subtype_blood2,levels = c("Control","RRMS","PMS"))
pls$subtype_blood3=as.character(pls$subtype_blood2)
pls$subtype_blood3[pls$subtype_blood2 %in% c("PMS","RRMS")]="MS"
pls$subtype_blood3=factor(pls$subtype_blood3,levels = c("Control","MS"))

tmp=pls %>% mutate(`Plasma p-tau181`=pTau181,
                   `Plasma p-tau217`=pTau217,
                   `Plasma t-tau`=Total.Tau,
                   `Plasma NfL`=NFL,
                   `Plasma GFAP`=GFAP,
                   Age=blood_age,
                   subtype=factor(subtype_blood3, levels = c("MS","Control")))
p_cor=ggpairs(data=tmp,columns = c(28,23:27),
              ggplot2::aes(colour = subtype,alpha=0.4),
              upper = list(continuous = wrap("cor", method = "spearman")))
p_cor=p_cor+theme_bw()
p_cor

tab=CreateTableOne(vars = c("Age","sex","race_ethnicity","Plasma p-tau181", "Plasma p-tau217","Plasma t-tau", "Plasma NfL","Plasma GFAP"),
                   strata = "subtype",
                   addOverall = T,
                   data=tmp,
                   factorVars = c("sex","race_ethnicity"),
                   includeNA = F,
                   argsApprox = list(correct = F), #chi-square test no correction
                   argsNormal = list(var.equal = T) # 2-sample t-test with equal var
)
tab=print(tab,
          exact=c("sex","race"),  #exact fisher test
          nonnormal = c("Age","Plasma p-tau181", "Plasma p-tau217","Plasma t-tau", "Plasma NfL","Plasma GFAP"),  
          catDigits = 2,
          contDigits = 2,
          pDigits = 3,
          formatOptions = list(big.mark = ","),
          showAllLevels = TRUE)

p1=ggplot(pls, aes(x=subtype_blood3, y=pTau181, fill=subtype_blood3)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs( x="",y = "Plasma p-tau181")+
  scale_fill_brewer(palette="Reds") + 
  theme_bw()+
  theme(legend.position = "none",
        axis.title.y = element_text(size=9, face = "bold"),
        axis.text.y = element_text(size=8, face = "bold"),
        axis.text.x = element_text(size=8, face = "bold") )
p2=ggplot(pls, aes(x=subtype_blood3, y=pTau217, fill=subtype_blood3)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs( x="",y = "Plasma p-tau217")+
  scale_fill_brewer(palette="Oranges") + 
  theme_bw()+
  theme(legend.position = "none",
        axis.title.y = element_text(size=9, face = "bold"),
        axis.text.y = element_text(size=8, face = "bold"),
        axis.text.x = element_text(size=8, face = "bold") )
p3=ggplot(pls, aes(x=subtype_blood3, y=Total.Tau, fill=subtype_blood3)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs( x="",y = "Plasma t-tau")+
  scale_fill_brewer(palette="Greens") + 
  theme_bw()+
  theme(legend.position = "none",
        axis.title.y = element_text(size=9, face = "bold"),
        axis.text.y = element_text(size=8, face = "bold"),
        axis.text.x = element_text(size=8, face = "bold") )
p4=ggplot(pls, aes(x=subtype_blood3, y=NFL, fill=subtype_blood3)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs( x="",y = "Plasma NfL")+
  scale_fill_brewer(palette="Blues") + 
  theme_bw()+
  theme(legend.position = "none",
        axis.title.y = element_text(size=9, face = "bold"),
        axis.text.y = element_text(size=8, face = "bold"),
        axis.text.x = element_text(size=8, face = "bold") )
p5=ggplot(pls, aes(x=subtype_blood3, y=GFAP, fill=subtype_blood3)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs( x="",y = "Plasma GFAP")+
  scale_fill_brewer(palette="Purples") + 
  theme_bw()+
  theme(legend.position = "none",
        axis.title.y = element_text(size=9, face = "bold"),
        axis.text.y = element_text(size=8, face = "bold"),
        axis.text.x = element_text(size=8, face = "bold") )

p1+p2+p3+p4+p5
wilcox.test(pls$pTau181 ~pls$subtype_blood3)
wilcox.test(pls$pTau217 ~pls$subtype_blood3)
wilcox.test(pls$Total.Tau ~pls$subtype_blood3)
wilcox.test(pls$NFL~pls$subtype_blood3)
wilcox.test(pls$GFAP ~pls$subtype_blood3)
pls$MS=ifelse(pls$subtype_blood3=="MS",1,0)
pls=pls %>% mutate(pTau181_scale=as.numeric(scale(pTau181)),
                   pTau217_scale=as.numeric(scale(pTau217)),
                   Total.Tau_scale=as.numeric(scale(Total.Tau)),
                   nfl_scale=as.numeric(scale(NFL)),
                   gfap_scale=as.numeric(scale(GFAP)),
                   lptau181=log(pTau181),
                   lptau217=log(pTau217),
                   lttau=log(Total.Tau),
                   lnfl=log(NFL),
                   lgfap=log(GFAP))
m1=glm(MS ~ pTau181_scale+blood_age+sex, data=pls, family = binomial()); summary(m1)
m2=glm(MS ~ pTau217_scale+blood_age+sex, data=pls, family = binomial()); summary(m2)
m3=glm(MS ~ Total.Tau_scale+blood_age+sex, data=pls, family = binomial()); summary(m3)
m4=glm(MS ~ nfl_scale+blood_age+sex, data=pls, family = binomial()); summary(m4)
m5=glm(MS ~ gfap_scale+blood_age+sex, data=pls, family = binomial()); summary(m5)

m1=lm(pTau181 ~ MS+blood_age+sex, data=pls); summary(m1)
m2=lm(pTau217 ~ MS+blood_age+sex, data=pls); summary(m2)
m3=lm(Total.Tau ~ MS+blood_age+sex, data=pls); summary(m3)
m4=lm(NFL ~ MS+blood_age+sex, data=pls); summary(m4)
m5=lm(GFAP ~ MS+blood_age+sex, data=pls); summary(m5)

m1=lm(lptau181 ~ MS+blood_age+sex, data=pls); summary(m1)
m2=lm(lptau217 ~ MS+blood_age+sex, data=pls); summary(m2)
m3=lm(lttau ~ MS+blood_age+sex, data=pls); summary(m3)
m4=lm(lnfl ~ MS+blood_age+sex, data=pls); summary(m4)
m5=lm(lgfap ~ MS+blood_age+sex, data=pls); summary(m5)



pls=pls %>% mutate( age_c=ifelse(blood_age<=45,"<=45 yrs",ifelse(blood_age<=60,"45-60 yrs",">60 yrs")),
                    age_c=factor(age_c,levels = c("<=45 yrs","45-60 yrs",">60 yrs")))

m1=lm(pTau181 ~ MS+age_c+sex, data=pls); summary(m1)
m2=lm(pTau217 ~ MS+age_c+sex, data=pls); summary(m2)
m3=lm(Total.Tau ~ MS+age_c+sex, data=pls); summary(m3)
m4=lm(NFL ~ MS+age_c+sex, data=pls); summary(m4)
m5=lm(GFAP ~ MS+age_c+sex, data=pls); summary(m5)




#### Plasma Tau and MS type ####
p1=ggplot(pls, aes(x=subtype_blood2, y=pTau181, fill=subtype_blood2)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs( x="",y = "Plasma p-tau181")+
  scale_fill_brewer(palette="Reds") + 
  theme_bw()+
  theme(legend.position = "none",
        axis.title.y = element_text(size=9, face = "bold"),
        axis.text.y = element_text(size=8, face = "bold"),
        axis.text.x = element_text(size=8, face = "bold") )
p2=ggplot(pls, aes(x=subtype_blood2, y=pTau217, fill=subtype_blood2)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs( x="",y = "Plasma p-tau217")+
  scale_fill_brewer(palette="Oranges") + 
  theme_bw()+
  theme(legend.position = "none",
        axis.title.y = element_text(size=9, face = "bold"),
        axis.text.y = element_text(size=8, face = "bold"),
        axis.text.x = element_text(size=8, face = "bold") )
p3=ggplot(pls, aes(x=subtype_blood2, y=Total.Tau, fill=subtype_blood2)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs( x="",y = "Plasma t-tau")+
  scale_fill_brewer(palette="Greens") + 
  theme_bw()+
  theme(legend.position = "none",
        axis.title.y = element_text(size=9, face = "bold"),
        axis.text.y = element_text(size=8, face = "bold"),
        axis.text.x = element_text(size=8, face = "bold") )
p4=ggplot(pls, aes(x=subtype_blood2, y=NFL, fill=subtype_blood2)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs( x="",y = "Plasma NfL")+
  scale_fill_brewer(palette="Blues") + 
  theme_bw()+
  theme(legend.position = "none",
        axis.title.y = element_text(size=9, face = "bold"),
        axis.text.y = element_text(size=8, face = "bold"),
        axis.text.x = element_text(size=8, face = "bold") )
p5=ggplot(pls, aes(x=subtype_blood2, y=GFAP, fill=subtype_blood2)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs( x="",y = "Plasma GFAP")+
  scale_fill_brewer(palette="Purples") + 
  theme_bw()+
  theme(legend.position = "none",
        axis.title.y = element_text(size=9, face = "bold"),
        axis.text.y = element_text(size=8, face = "bold"),
        axis.text.x = element_text(size=8, face = "bold") )

p1+p2+p3+p4+p5
kruskal.test(pls$pTau181 ~pls$subtype_blood2)
kruskal.test(pls$pTau217 ~pls$subtype_blood3)
kruskal.test(pls$Total.Tau ~pls$subtype_blood3)
kruskal.test(pls$NFL ~pls$subtype_blood3)
kruskal.test(pls$GFAP ~pls$subtype_blood3)

dunn.test(pls$pTau181,pls$subtype_blood2,kw=T,method = "bh")
dunn.test(pls$pTau217,pls$subtype_blood2,kw=T,method = "bh")
dunn.test(pls$Total.Tau,pls$subtype_blood2,kw=T,method = "bh")
dunn.test(pls$NFL,pls$subtype_blood2,kw=T,method = "bh")
dunn.test(pls$GFAP,pls$subtype_blood2,kw=T,method = "bh")

m <- multinom(subtype_blood2 ~ pTau181_scale+sex, data = pls)
summary(m)
exp(summary(m)$coefficients)
exp(summary(m)$coefficients-1.96*summary(m)$standard.errors)
exp(summary(m)$coefficients+1.96*summary(m)$standard.errors)
z <- summary(m)$coefficients/summary(m)$standard.errors
(1 - pnorm(abs(z), 0, 1)) * 2

m <- multinom(subtype_blood2 ~ pTau217_scale+sex, data = pls)
summary(m)
exp(summary(m)$coefficients)
exp(summary(m)$coefficients-1.96*summary(m)$standard.errors)
exp(summary(m)$coefficients+1.96*summary(m)$standard.errors)
z <- summary(m)$coefficients/summary(m)$standard.errors
(1 - pnorm(abs(z), 0, 1)) * 2

m <- multinom(subtype_blood2 ~ Total.Tau_scale+sex, data = pls)
summary(m)
exp(summary(m)$coefficients)
exp(summary(m)$coefficients-1.96*summary(m)$standard.errors)
exp(summary(m)$coefficients+1.96*summary(m)$standard.errors)
z <- summary(m)$coefficients/summary(m)$standard.errors
(1 - pnorm(abs(z), 0, 1)) * 2

m <- multinom(subtype_blood2 ~ nfl_scale+sex, data = pls)
summary(m)
exp(summary(m)$coefficients)
exp(summary(m)$coefficients-1.96*summary(m)$standard.errors)
exp(summary(m)$coefficients+1.96*summary(m)$standard.errors)
z <- summary(m)$coefficients/summary(m)$standard.errors
(1 - pnorm(abs(z), 0, 1)) * 2

m <- multinom(subtype_blood2 ~ gfap_scale+sex, data = pls)
summary(m)
exp(summary(m)$coefficients)
exp(summary(m)$coefficients-1.96*summary(m)$standard.errors)
exp(summary(m)$coefficients+1.96*summary(m)$standard.errors)
z <- summary(m)$coefficients/summary(m)$standard.errors
(1 - pnorm(abs(z), 0, 1)) * 2



m1=lm(pTau181 ~subtype_blood2+sex ,data=pls ); summary(m1)
m2=lm(NFL ~subtype_blood2+blood_age+sex ,data=pls ); summary(m2)

#### Plasma Tau and PMS vs. RRMS ####

pls_PMS=pls %>% filter(subtype_blood1 %in% c("SPMS","PPMS","RRMS")) %>% 
  mutate(outcome=ifelse(subtype_blood1 =="RRMS",0,1))
pls_PMS=pls_PMS %>% mutate(pTau181_scale2=as.numeric(scale(pTau181)),
                   pTau217_scale2=as.numeric(scale(pTau217)),
                   Total.Tau_scale2=as.numeric(scale(Total.Tau)),
                   nfl_scale2=as.numeric(scale(NFL)),
                   gfap_scale2=as.numeric(scale(GFAP)))
pls_PMS $duration=pls_PMS$blood_age -pls_PMS$age_ms_dx
pls_PMS $duration[which(pls_PMS$duration<0)]=0
summary(pls_PMS$duration)
pls_PMS=pls_PMS %>% mutate(duration_c=ifelse(duration <=5,"<=5 yrs", ifelse(duration<=15,"5-15 yrs",">15 yrs")),
                           age_c=ifelse(blood_age<=45,"<=45 yrs",ifelse(blood_age<=60,"45-60 yrs",">60 yrs")),
                           duration_c=factor(duration_c,levels = c("<=5 yrs","5-15 yrs",">15 yrs")),
                           age_c=factor(age_c,levels = c("<=45 yrs","45-60 yrs",">60 yrs")))



match_obj <- matchit(outcome ~ sex + age_c+race_ethnicity+duration_c,
                     data = pls_PMS, method="full",distance = "glm",
                     ratio =1,
                     replace = F)
plot(summary(match_obj), abs = FALSE)
datanew=match.data(match_obj)
table(datanew$outcome)

Quantile(datanew$blood_age[datanew$outcome==1])
Quantile(datanew$blood_age[datanew$outcome==0])
Quantile(datanew$blood_age[datanew$outcome==0], 
         datanew$weights[datanew$outcome==0])

Quantile(datanew$duration[datanew$outcome==1])
Quantile(datanew$duration[datanew$outcome==0])
Quantile(datanew$duration[datanew$outcome==0], 
         datanew$weights[datanew$outcome==0])

Quantile(datanew$GFAP[datanew$outcome==1])
Quantile(datanew$GFAP[datanew$outcome==0])
Quantile(datanew$GFAP[datanew$outcome==0], 
         datanew$weights[datanew$outcome==0])

table(datanew$age_c[datanew$outcome==1])
table(datanew$age_c[datanew$outcome==0])
tmp_prop=ifelse(datanew$outcome==0 & datanew$age_c=="<=45 yrs",1,0)
sum(tmp_prop*datanew$weights); sum(tmp_prop*datanew$weights)/117
tmp_prop=ifelse(datanew$outcome==0 & datanew$age_c=="45-60 yrs",1,0)
sum(tmp_prop*datanew$weights);sum(tmp_prop*datanew$weights)/117
tmp_prop=ifelse(datanew$outcome==0 & datanew$age_c==">60 yrs",1,0)
sum(tmp_prop*datanew$weights);sum(tmp_prop*datanew$weights)/117

table(datanew$race_ethnicity[datanew$outcome==1])
table(datanew$race_ethnicity[datanew$outcome==0])
tmp_prop=ifelse(datanew$outcome==0 & datanew$race_ethnicity=="Non-Hispanic White",1,0)
sum(tmp_prop*datanew$weights); sum(tmp_prop*datanew$weights)/117
tmp_prop=ifelse(datanew$outcome==0 & datanew$race_ethnicity=="Other",1,0)
sum(tmp_prop*datanew$weights);sum(tmp_prop*datanew$weights)/117

table(datanew$sex[datanew$outcome==1])
table(datanew$sex[datanew$outcome==0])
tmp_prop=ifelse(datanew$outcome==0 & datanew$sex=="M",1,0)
sum(tmp_prop*datanew$weights); sum(tmp_prop*datanew$weights)/117
tmp_prop=ifelse(datanew$outcome==0 & datanew$sex=="F",1,0)
sum(tmp_prop*datanew$weights);sum(tmp_prop*datanew$weights)/117

table(datanew$duration_c[datanew$outcome==1])
table(datanew$duration_c[datanew$outcome==0])
tmp_prop=ifelse(datanew$outcome==0 & datanew$duration_c=="<=5 yrs",1,0)
sum(tmp_prop*datanew$weights); sum(tmp_prop*datanew$weights)/117
tmp_prop=ifelse(datanew$outcome==0 & datanew$duration_c=="5-15 yrs",1,0)
sum(tmp_prop*datanew$weights);sum(tmp_prop*datanew$weights)/117
tmp_prop=ifelse(datanew$outcome==0 & datanew$duration_c==">15 yrs",1,0)
sum(tmp_prop*datanew$weights);sum(tmp_prop*datanew$weights)/117


m1=glm(outcome ~ pTau181_scale+duration+I(duration^2),data=datanew, weights = weights, family = binomial())
summary(m1);exp(coef(m1))
m2=glm(outcome ~ pTau217_scale+duration+I(duration^2),data=datanew, weights = weights, family = binomial())
summary(m2);exp(coef(m2)); exp(confint(m2, "pTau217_scale"))
m3=glm(outcome ~ Total.Tau_scale+duration+I(duration^2),data=datanew, weights = weights, family = binomial())
summary(m3);exp(coef(m3));exp(confint(m3, "Total.Tau_scale"))
m4=glm(outcome ~ nfl_scale+duration+I(duration^2),data=datanew, weights = weights, family = binomial())
summary(m4);exp(coef(m4));exp(confint(m4, "nfl_scale"))
m5=glm(outcome ~ gfap_scale+duration+I(duration^2),data=datanew, weights = weights, family = binomial())
summary(m5);exp(coef(m5));exp(confint(m5, "gfap_scale"))


m1=glm(outcome ~ pTau181_scale+sex+race_ethnicity+blood_age+duration,data=pls_PMS,  family = binomial()); summary(m1)
m2=glm(outcome ~ pTau217_scale+sex+race_ethnicity+blood_age+duration,data=pls_PMS,  family = binomial()); summary(m2)
m3=glm(outcome ~ Total.Tau_scale+sex+race_ethnicity+blood_age+duration,data=pls_PMS,  family = binomial()); summary(m3)
m4=glm(outcome ~ nfl_scale+sex+race_ethnicity+blood_age+duration,data=pls_PMS,  family = binomial()); summary(m4)
m5=glm(outcome ~ gfap_scale+sex+race_ethnicity+blood_age+duration,data=pls_PMS,  family = binomial()); summary(m5)


m_0=glm(outcome ~blood_age+sex+race_ethnicity+duration,data=pls_PMS,  family = binomial())
summary(m_0)
rocobj=roc(pls_PMS$outcome==1, predict(m_0,type="response"),ci = T)
auc(rocobj)

m_1=glm(outcome ~pTau181_scale+pTau217_scale+Total.Tau_scale,data=pls_PMS,  family = binomial())
summary(m_1)
rocobj=roc(pls_PMS$outcome==1, predict(m_1,type="response"),ci = T)
auc(rocobj)

m_2=glm(outcome ~nfl_scale+gfap_scale,data=pls_PMS,  family = binomial())
summary(m_2)
rocobj=roc(pls_PMS$outcome==1, predict(m_2,type="response"),ci = T)
auc(rocobj)

m_3=glm(outcome ~blood_age+sex+race_ethnicity+duration+
          pTau181_scale+pTau217_scale+Total.Tau_scale,data=pls_PMS,  family = binomial())
summary(m_3)
rocobj=roc(pls_PMS$outcome==1, predict(m_3,type="response"),ci = T)
auc(rocobj)
m_4=glm(outcome ~blood_age+sex+race_ethnicity+duration+
          nfl_scale+gfap_scale,data=pls_PMS,  family = binomial())
summary(m_4)
rocobj=roc(pls_PMS$outcome==1, predict(m_4,type="response"),ci = T)
auc(rocobj)

m_5=glm(outcome ~pTau181_scale+pTau217_scale+Total.Tau_scale+
          nfl_scale+gfap_scale,data=pls_PMS,  family = binomial())
summary(m_5)
rocobj=roc(pls_PMS$outcome==1, predict(m_5,type="response"),ci = T)
auc(rocobj)

m_6=glm(outcome ~blood_age+sex+race_ethnicity+duration+pTau181_scale+pTau217_scale+Total.Tau_scale+
          nfl_scale+gfap_scale,data=pls_PMS,  family = binomial())
summary(m_6)
rocobj=roc(pls_PMS$outcome==1, predict(m_6,type="response"),ci = T)
auc(rocobj)



m_full=glm(outcome ~age_c+sex+race_ethnicity+duration_c+ pTau181_scale+pTau217_scale,data=pls_PMS,  family = binomial())
summary(m_full)
rocobj=roc(pls_PMS$outcome==1, predict(m_full,type="response"),ci = T)
auc(rocobj)
ciobj <- ci.se(rocobj, specificities=seq(0, 1, l=25))




set.seed(1234)
set.seed(2023)
set.seed(2018)
index <- createDataPartition(pls_PMS$outcome, p=0.7, list=FALSE, times=1)
train_df <- pls_PMS[index,]; nrow (train_df)
train_df$outcome=factor(train_df$outcome,levels = c(0,1),labels = c("RRMS","PMS"))
test_df <- pls_PMS[-index,]; nrow (test_df)
test_df$outcome=factor(test_df$outcome,levels = c(0,1),labels = c("RRMS","PMS"))

ctrlspecs <- trainControl(method="repeatedcv", 
                          number=5, 
                          repeats = 20,
                          savePredictions="all",
                          classProbs=TRUE)

model0 <- train(outcome ~ blood_age+sex+race_ethnicity+duration, data=train_df,
                method="glm", 
                family=binomial, 
                trControl=ctrlspecs)
print(model0)
summary(model0)
varImp(model0)
predictions <- predict(model0, newdata=test_df)
confusionMatrix(data=predictions, test_df$outcome,positive = "PMS")


model1 <- train(outcome ~ pTau181_scale+pTau217_scale+Total.Tau_scale, data=train_df,
                method="glm", 
                family=binomial, 
                trControl=ctrlspecs)
print(model1)
summary(model1)
varImp(model1)
predictions <- predict(model1, newdata=test_df)
confusionMatrix(data=predictions, test_df$outcome,positive = "PMS")

model2 <- train(outcome ~ pTau181_scale+pTau217_scale+Total.Tau_scale, data=train_df,
                method="glm", 
                family=binomial, 
                trControl=ctrlspecs)
print(model2)
summary(model2)
varImp(model2)
predictions <- predict(model2, newdata=test_df)
confusionMatrix(data=predictions, test_df$outcome,positive = "PMS")

model3 <- train(outcome ~blood_age+sex+race_ethnicity+duration+ pTau181_scale+pTau217_scale+Total.Tau_scale, data=train_df,
                method="glm", 
                family=binomial, 
                trControl=ctrlspecs)
print(model3)
summary(model3)
varImp(model3)
predictions <- predict(model3, newdata=test_df)
confusionMatrix(data=predictions, test_df$outcome,positive = "PMS")

model4 <- train(outcome ~blood_age+sex+race_ethnicity+duration+nfl_scale+gfap_scale , data=train_df,
                method="glm", 
                family=binomial, 
                trControl=ctrlspecs)
print(model4)
summary(model4)
varImp(model4)
predictions <- predict(model4, newdata=test_df)
confusionMatrix(data=predictions, test_df$outcome,positive = "PMS")

model5 <- train(outcome ~nfl_scale+gfap_scale, data=train_df,
                method="glm", 
                family=binomial, 
                trControl=ctrlspecs)
print(model5)
summary(model5)
varImp(model5)
predictions <- predict(model5, newdata=test_df)
confusionMatrix(data=predictions, test_df$outcome,positive = "PMS")

model6 <- train(outcome ~blood_age+sex+race_ethnicity+duration+ pTau181_scale+pTau217_scale+Total.Tau_scale+nfl_scale+gfap_scale, data=train_df,
                method="glm", 
                family=binomial, 
                trControl=ctrlspecs)

predictions <- predict(model6, newdata=test_df)
confusionMatrix(data=predictions, test_df$outcome,positive = "PMS")

###### CSF tau in controls and MS #####
csf=read.csv("~/OneDrive - University of Pittsburgh/XiaLab/ADBiomarker/csf_CleanData_230828.csv")
csf=csf %>% filter(subtype_csf !="Related disorder")
csf=csf %>% mutate(race_ethnicity=ifelse(race=="White" & ethnicity=="Not Hispanic or Latino","Non-Hispanic White","Other" ))
csf=csf %>% mutate(subtype_csf1=ifelse(subtype_csf %in% c("CIS","RIS","RRMS"),"RRMS",
                                       ifelse(subtype_csf=="PPMS","PMS", "Control")))
csf$subtype_csf1=factor(csf$subtype_csf1, levels = c("Control","RRMS","PMS"),
                        labels = c("Control (n=3)","RRMS (n=19)","PMS (n=3)"))
p1=ggplot(csf, aes(x=subtype_csf1, y=pTau181, fill=subtype_csf1)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs( x="",y = "CSF p-tau181")+
  scale_fill_brewer(palette="Oranges") + 
  theme_classic()+
  theme(legend.position = "none",
        axis.title.y = element_text(size=9, face = "bold"),
        axis.text.y = element_text(size=8, face = "bold"),
        axis.text.x = element_text(size=8, face = "bold") )
p2=ggplot(csf, aes(x=subtype_csf1, y=pTau217, fill=subtype_csf1)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs( x="",y = "CSF p-tau217")+
  scale_fill_brewer(palette="Purples") + 
  theme_classic()+
  theme(legend.position = "none",
        axis.title.y = element_text(size=9, face = "bold"),
        axis.text.y = element_text(size=8, face = "bold"),
        axis.text.x = element_text(size=8, face = "bold") )
p3=ggplot(csf, aes(x=subtype_csf1, y=Total.Tau, fill=subtype_csf1)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs( x="",y = "CSF t-tau")+
  scale_fill_brewer(palette="Greys") + 
  theme_classic()+
  theme(legend.position = "none",
        axis.title.y = element_text(size=9, face = "bold"),
        axis.text.y = element_text(size=8, face = "bold"),
        axis.text.x = element_text(size=8, face = "bold") )

p1+p2+p3
kruskal.test(csf$pTau181 ~csf$subtype_csf1)
kruskal.test(csf$pTau217 ~csf$subtype_csf1)
kruskal.test(csf$Total.Tau ~csf$subtype_csf1)







################################################################################################
################################################################################################
################################################################################################

setwd("~/OneDrive - University of Pittsburgh/XiaLab/ADBiomarker/manuscript/20240903/")
pls=read.csv("../../Data_Archive/plasma_CleanData_MS_231206.csv")
#### longitudinal OCT ####
oct=read.csv("data/oct.csv")
oct=oct %>% filter(id_participant %in% pls$id_participant)
length(unique(oct$id_participant)) #84
data=left_join(pls,oct,by="id_participant")
data=data %>% mutate(diff=as.numeric((as.Date(report1_date)-as.Date(blood_date))/365.25))
summary(data$diff)
data=data %>% filter(!is.na(diff))
table(data$subtype,useNA="ifany")
data=data %>% mutate(type=factor(subtype, level=c("RRMS","PMS")))
data=data %>% filter(diff>=(-1)) 
length(unique(data$id_participant))
data=data %>% group_by(id_participant) %>% mutate(time=diff-min(diff))
data_left=data[,c("id_participant","report1_date",colnames(data)[grepl("left",colnames(data))])]
colnames(data_left)=c("id_participant","report1_date","rnfl","macular","mac_thickness","min_gcl","avg_gcl","max_gcl")
data_left$eye="L"
data_right=data[,c("id_participant","report1_date",colnames(data)[grepl("right",colnames(data))])]
colnames(data_right)=c("id_participant","report1_date","rnfl","macular","mac_thickness","min_gcl","avg_gcl","max_gcl")
data_right$eye="R"
data_eye=rbind(data_right,data_left)
tmp=data [!duplicated(data$id_participant),c("id_participant","pTau181","pTau217","Total.Tau","GFAP","NFL")]
tmp$pTau181_scale=as.numeric(scale(tmp$pTau181))
tmp$pTau217_scale=as.numeric(scale(tmp$pTau217))
tmp$Total.Tau_scale=as.numeric(scale(tmp$Total.Tau))
tmp$gfap_scale=as.numeric(scale(tmp$GFAP))
tmp$nfl_scale=as.numeric(scale(tmp$NFL))
tmp$ptau181_tert=ntile(tmp$pTau181,2)
tmp$ptau217_tert=ntile(tmp$pTau217,2)
tmp$ttau_tert=ntile(tmp$Total.Tau,2)
tmp$nfl_tert=ntile(tmp$NFL,2)
tmp$gfap_tert=ntile(tmp$GFAP,2)
                   
data2=left_join(data[,c("id_participant", "blood_age","sex","race_ethnicity","duration","blood_date","report1_date","time","diff","type")],
               tmp, by="id_participant")

data_long=left_join(data_eye, data2[,c("id_participant","blood_date","report1_date","time","diff","blood_age","sex","type","race_ethnicity","duration",
                                      "pTau181_scale","pTau217_scale","Total.Tau_scale","nfl_scale"  ,"gfap_scale",
                                      "ptau181_tert","ptau217_tert","ttau_tert","nfl_tert","gfap_tert")],
                    by=c("id_participant","report1_date"))

data_long=data_long %>% arrange(id_participant, eye,diff)
data_long$eye=factor(data_long$eye)
data_long$id_participant=factor(data_long$id_participant)

outcome=c("rnfl","macular","mac_thickness","min_gcl","avg_gcl")
marker=c("pTau181_scale","pTau217_scale","Total.Tau_scale","nfl_scale","gfap_scale")
covar="type+blood_age+sex+race_ethnicity+duration"
gee_result=data.frame(outcome=c(), marker=c(), beta=c(), se=c(), p=c())
i=0
for (o in outcome){
  for (m in marker){
    f=as.formula(paste(o,paste0(m,"+",covar), sep="~"))
    i=i+1
    gee_result[i,1]=o
    gee_result[i,2]=m
    gee<- gee(f,
              data = data_long,
              id = id_participant,
              family = gaussian(),
              corstr = "exchangeable")
    gee_result[i,3:5]=c(gee$coefficients[[2]], sqrt(gee$naive.variance[2,2]),2*pnorm(abs(gee$coefficients[[2]]/sqrt(gee$naive.variance[2,2])),lower.tail = F))
    
  }
}


#### functional tests ####
#### longitudianal SDMT ####
sdmt=read.csv("data/FunctionalTest.csv")
sdmt=sdmt %>% filter(id_participant %in% pls$id_participant)
sdmt$sdmt_correct=substr(sdmt$sdmt_score,1,2)
sdmt$sdmt_correct=as.numeric(sdmt$sdmt_correct)
summary(sdmt$sdmt_correct)
sdmt=sdmt %>% filter(is.na(sdmt_correct)==F) 
data=left_join(pls,sdmt,by="id_participant")
data=data %>% filter(is.na(sdmt_correct)==F)
data=data %>% mutate(diff=(as.Date(sdmt_date)-as.Date(blood_date))/365.25)
data=data %>% filter(diff>=-1) 
data=data %>% group_by(id_participant) %>% mutate(time=diff-min(diff))

data=data %>%
  mutate(type=factor(subtype, levels = c("RRMS","PMS")))

tmp=data [!duplicated(data$id_participant),c("id_participant","pTau181","pTau217","Total.Tau","GFAP","NFL")]
tmp$pTau181_scale=as.numeric(scale(tmp$pTau181))
tmp$pTau217_scale=as.numeric(scale(tmp$pTau217))
tmp$Total.Tau_scale=as.numeric(scale(tmp$Total.Tau))
tmp$gfap_scale=as.numeric(scale(tmp$GFAP))
tmp$nfl_scale=as.numeric(scale(tmp$NFL))
tmp$ptau181_tert=ntile(tmp$pTau181,2)
tmp$ptau217_tert=ntile(tmp$pTau217,2)
tmp$ttau_tert=ntile(tmp$Total.Tau,2)
tmp$nfl_tert=ntile(tmp$NFL,2)
tmp$gfap_tert=ntile(tmp$GFAP,2)
data=left_join(data[,c("id_participant", "blood_age","sex","race_ethnicity","duration","blood_date","sdmt_date","diff","time","type","sdmt_correct")],
                tmp, by="id_participant")
data=data %>% filter(!is.na(duration))
data$id_participant=as.factor(data$id_participant)
length(unique(data$id_participant))
outcome="sdmt_correct"
marker=c("pTau181_scale","pTau217_scale","Total.Tau_scale","nfl_scale","gfap_scale")
covar="type+blood_age+sex+race_ethnicity+duration"
gee_result=data.frame( marker=c(), beta=c(), se=c(), p=c())
i=0
for (m in marker){
    f=as.formula(paste(outcome,paste0(m,"+",covar), sep="~"))
    i=i+1
    gee_result[i,1]=m
    gee<- geeglm(f,
              data = data,
              id = id_participant,
              family = gaussian(),
              corstr = "unstructured"
              )
    gee_result[i,2:4]=summary(gee)$coef[2,c(1,2,4)]
}

#### longitudinal t25fw ####
t25=read.csv("data/FunctionalTest.csv")
t25=t25 %>% filter(id_participant %in% pls$id_participant)
summary(as.numeric(t25$t25w1));summary(as.numeric(t25$t25w2))
t25=t25 %>% mutate(t25w1=as.numeric(t25w1),
                   t25w2=as.numeric(t25w2))
t25=t25 %>% mutate(t25w_avg=ifelse(is.na(t25w1)==T & is.na(t25w2)==T, NA,
                                   ifelse(is.na(t25w1)==T & is.na(t25w2)==F, t25w2,
                                          ifelse(is.na(t25w1)==F & is.na(t25w2)==T, t25w1, (t25w1+t25w2)/2))))
t25=t25 %>% filter(is.na(t25w_avg)==F) 
data=left_join(pls,t25,by="id_participant")
data=data %>% filter(is.na(t25w_avg)==F)
data=data %>% mutate(diff=(as.Date(footwalk_date)-as.Date(blood_date))/365.25)
data=data %>% filter(diff >= (-1))
length(unique(data$id_participant))
data=data %>% group_by(id_participant) %>% mutate(time=diff-min(diff))

data=data  %>%
  mutate(type=factor(subtype, levels = c("RRMS","PMS")))

tmp=data [!duplicated(data$id_participant),c("id_participant","pTau181","pTau217","Total.Tau","GFAP","NFL")]
tmp$pTau181_scale=as.numeric(scale(tmp$pTau181))
tmp$pTau217_scale=as.numeric(scale(tmp$pTau217))
tmp$Total.Tau_scale=as.numeric(scale(tmp$Total.Tau))
tmp$gfap_scale=as.numeric(scale(tmp$GFAP))
tmp$nfl_scale=as.numeric(scale(tmp$NFL))
tmp$ptau181_tert=ntile(tmp$pTau181,2)
tmp$ptau217_tert=ntile(tmp$pTau217,2)
tmp$ttau_tert=ntile(tmp$Total.Tau,2)
tmp$nfl_tert=ntile(tmp$NFL,2)
tmp$gfap_tert=ntile(tmp$GFAP,2)
data=left_join(data[,c("id_participant", "blood_age","sex","race_ethnicity","duration","blood_date","footwalk_date","diff","time","type","t25w_avg")],
               tmp, by="id_participant")
data=data %>% filter(!is.na(duration))
data$id_participant=as.factor(data$id_participant)
data$lt25w_avg=log(data$t25w_avg)
length(unique(data$id_participant))
outcome="t25w_avg"
marker=c("pTau181_scale","pTau217_scale","Total.Tau_scale","nfl_scale","gfap_scale")
covar="type+blood_age+sex+race_ethnicity+duration"
gee_result=data.frame( marker=c(), beta=c(), se=c(), p=c())
i=0
for (m in marker){
  f=as.formula(paste(outcome,paste0(m,"+",covar), sep="~"))
  i=i+1
  gee_result[i,1]=m
  gee<- gee(f,
               data = data,
               id = id_participant,
               family = gaussian(),
               corstr = "unstructured"
  )
  gee_result[i,2:4]=c(gee$coefficients[[2]], sqrt(gee$naive.variance[2,2]),2*pnorm(abs(gee$coefficients[[2]]/sqrt(gee$naive.variance[2,2])),lower.tail = F))
  
}

#### longitudinal peg ####
peg=read.csv("data/FunctionalTest.csv")
peg=peg %>% filter(id_participant %in% pls$id_participant)
summary(as.numeric(peg$peg_dom1));summary(as.numeric(peg$peg_dom2))
peg=peg %>% mutate(peg_dom1=as.numeric(peg_dom1),
                   peg_dom2=as.numeric(peg_dom2))

peg=peg %>% mutate(peg_avg=ifelse(is.na(peg_dom1)==T & is.na(peg_dom2)==T, NA,
                                  ifelse(is.na(peg_dom1)==T & is.na(peg_dom2)==F, peg_dom2,
                                         ifelse(is.na(peg_dom1)==F & is.na(peg_dom2)==T, peg_dom1, (peg_dom1+peg_dom2)/2))))
peg=peg %>% filter(is.na(peg_avg)==F) 
data=left_join(pls,peg,by="id_participant")
data=data %>% filter(is.na(peg_avg)==F)
data=data %>% mutate(diff=(as.Date(hpt_date)-as.Date(blood_date))/365.25)
data=data %>% filter(diff >= (-1))
data=data %>% group_by(id_participant) %>% mutate(time=diff-min(diff))
data=data %>% 
  mutate(type=factor(subtype, levels = c("RRMS","PMS")))

tmp=data [!duplicated(data$id_participant),c("id_participant","pTau181","pTau217","Total.Tau","GFAP","NFL")]
tmp$pTau181_scale=as.numeric(scale(tmp$pTau181))
tmp$pTau217_scale=as.numeric(scale(tmp$pTau217))
tmp$Total.Tau_scale=as.numeric(scale(tmp$Total.Tau))
tmp$gfap_scale=as.numeric(scale(tmp$GFAP))
tmp$nfl_scale=as.numeric(scale(tmp$NFL))
tmp$ptau181_tert=ntile(tmp$pTau181,2)
tmp$ptau217_tert=ntile(tmp$pTau217,2)
tmp$ttau_tert=ntile(tmp$Total.Tau,2)
tmp$nfl_tert=ntile(tmp$NFL,2)
tmp$gfap_tert=ntile(tmp$GFAP,2)
data=left_join(data[,c("id_participant", "blood_age","sex","race_ethnicity","duration","blood_date","hpt_date","diff","time","type","peg_avg")],
               tmp, by="id_participant")
data=data %>% filter(!is.na(duration))
data$id_participant=as.factor(data$id_participant)

outcome="peg_avg"
marker=c("pTau181_scale","pTau217_scale","Total.Tau_scale","nfl_scale","gfap_scale")
covar="type+blood_age+sex+race_ethnicity+duration"
gee_result=data.frame( marker=c(), beta=c(), se=c(), p=c())
i=0
for (m in marker){
  f=as.formula(paste(outcome,paste0(m,"+",covar), sep="~"))
  i=i+1
  gee_result[i,1]=m
  gee<- gee(f,
            data = data,
            id = id_participant,
            family = gaussian(),
            corstr = "unstructured"
  )
  gee_result[i,2:4]=c(gee$coefficients[[2]], sqrt(gee$naive.variance[2,2]),2*pnorm(abs(gee$coefficients[[2]]/sqrt(gee$naive.variance[2,2])),lower.tail = F))
  
}




##### longitudinal PDDS ####
pdds=read.csv("data/pdds.csv")
pdds=pdds %>% filter(id_participant %in% pls$id_participant)
pdds=pdds %>% filter(is.na(PDDS)==F) 
data=left_join(pls,pdds,by="id_participant")
data=data %>% filter(is.na(PDDS)==F)
data=data %>% mutate(diff=(as.Date(pdds_date)-as.Date(blood_date))/365.25)


data=data %>% filter(diff >= (-1))
data=data %>% group_by(id_participant) %>% mutate(time=diff-min(diff))
data=data %>% 
  mutate(type=factor(subtype, levels = c("RRMS","PMS")))

tmp=data [!duplicated(data$id_participant),c("id_participant","pTau181","pTau217","Total.Tau","GFAP","NFL")]
tmp$pTau181_scale=as.numeric(scale(tmp$pTau181))
tmp$pTau217_scale=as.numeric(scale(tmp$pTau217))
tmp$Total.Tau_scale=as.numeric(scale(tmp$Total.Tau))
tmp$gfap_scale=as.numeric(scale(tmp$GFAP))
tmp$nfl_scale=as.numeric(scale(tmp$NFL))
tmp$ptau181_tert=ntile(tmp$pTau181,2)
tmp$ptau217_tert=ntile(tmp$pTau217,2)
tmp$ttau_tert=ntile(tmp$Total.Tau,2)
tmp$nfl_tert=ntile(tmp$NFL,2)
tmp$gfap_tert=ntile(tmp$GFAP,2)
data=left_join(data[,c("id_participant", "blood_age","sex","race_ethnicity","duration","blood_date","pdds_date","diff","time","type","PDDS")],
               tmp, by="id_participant")
data=data %>% filter(!is.na(duration))
data$id_participant=as.factor(data$id_participant)

outcome="PDDS"
marker=c("pTau181_scale","pTau217_scale","Total.Tau_scale","nfl_scale","gfap_scale")
covar="type+blood_age+sex+race_ethnicity+duration"
gee_result=data.frame( marker=c(), beta=c(), se=c(), p=c())
i=0
for (m in marker){
  f=as.formula(paste(outcome,paste0(m,"+",covar), sep="~"))
  i=i+1
  gee_result[i,1]=m
  gee<- gee(f,
            data = data,
            id = id_participant,
            family = gaussian(),
            corstr = "exchangeable"
  )
  gee_result[i,2:4]=c(gee$coefficients[[2]], sqrt(gee$naive.variance[2,2]),2*pnorm(abs(gee$coefficients[[2]]/sqrt(gee$naive.variance[2,2])),lower.tail = F))
  
}
gee_result_pdds=gee_result



setwd("~/OneDrive - University of Pittsburgh/XiaLab/ADBiomarker/manuscript/20240903/")
pls=read.csv("data/plasma_CleanData_MS_231206.csv")
#### longitudinal MRI ####
mri=read.csv("data/neuroquant.csv")
colnames(mri)[2]="id_participant"
mri=mri %>% filter(id_participant %in% pls$id_participant)
length(unique(mri$id_participant)) #63
data=left_join(pls,mri,by="id_participant")
data=data %>% mutate(diff=as.numeric((as.Date(mri_date)-as.Date(blood_date))/365.25))
summary(data$diff)
data=data %>% filter(!is.na(diff))
table(data$subtype,useNA="ifany")
data=data %>% mutate(type=factor(subtype, level=c("RRMS","PMS")))
data=data %>% filter(diff>=(0)) 
length(unique(data$id_participant))
data=data %>% group_by(id_participant) %>% mutate(time=diff-min(diff))

tmp=data [!duplicated(data$id_participant),c("id_participant","pTau181","pTau217","Total.Tau","GFAP","NFL")]
tmp$pTau181_scale=as.numeric(scale(tmp$pTau181))
tmp$pTau217_scale=as.numeric(scale(tmp$pTau217))
tmp$Total.Tau_scale=as.numeric(scale(tmp$Total.Tau))
tmp$gfap_scale=as.numeric(scale(tmp$GFAP))
tmp$nfl_scale=as.numeric(scale(tmp$NFL))
tmp$ptau181_tert=ntile(tmp$pTau181,3)
tmp$ptau217_tert=ntile(tmp$pTau217,3)
tmp$ttau_tert=ntile(tmp$Total.Tau,3)
tmp$nfl_tert=ntile(tmp$NFL,3)
tmp$gfap_tert=ntile(tmp$GFAP,3)

data2=left_join(data[,c("id_participant", "blood_age","sex","race_ethnicity","duration","blood_date","mri_date","time","diff","type")],
                tmp, by="id_participant")

data_long=left_join(data[,c("id_participant","mri_date","ICV","TBV","GM","TM","TBV_per","GM_per","TM_per")], 
                    data2[,c("id_participant","blood_date","mri_date","time","diff","blood_age","sex","type","race_ethnicity","duration",
                                       "pTau181_scale","pTau217_scale","Total.Tau_scale","nfl_scale"  ,"gfap_scale",
                                       "ptau181_tert","ptau217_tert","ttau_tert","nfl_tert","gfap_tert")],
                    by=c("id_participant","mri_date"))

data_long=data_long %>% arrange(id_participant, diff)

data_long$id_participant=factor(data_long$id_participant)

outcome=c("TBV","GM")
marker=c("pTau181_scale","pTau217_scale","Total.Tau_scale","nfl_scale","gfap_scale")
gee_result=data.frame(outcome=c(), marker=c(), beta=c(), se=c(), p=c())
i=0
for (o in outcome){
  for (m in marker){
    f=as.formula(paste(o,m, sep="~"))
    i=i+1
    gee_result[i,1]=o
    gee_result[i,2]=m
    gee<- gee(f,
              data = data_long,
              id = id_participant,
              family = gaussian(),
              corstr = "exchangeable")
    gee_result[i,3:5]=c(gee$coefficients[[2]], sqrt(gee$naive.variance[2,2]),2*pnorm(abs(gee$coefficients[[2]]/sqrt(gee$naive.variance[2,2])),lower.tail = F))
    
  }
}
cbind(gee_result$V3, gee_result$V3-1.96*sqrt(gee_result$V4), gee_result$V3+1.96*sqrt(gee_result$V4))/10
covar="type+blood_age+sex+race_ethnicity+duration+ICV"
gee_result=data.frame(outcome=c(), marker=c(), beta=c(), se=c(), p=c())
i=0
for (o in outcome){
  for (m in marker){
    f=as.formula(paste(o,paste0(m,"+",covar), sep="~"))
    i=i+1
    gee_result[i,1]=o
    gee_result[i,2]=m
    gee<- gee(f,
              data = data_long,
              id = id_participant,
              family = gaussian(),
              corstr = "exchangeable")
    gee_result[i,3:5]=c(gee$coefficients[[2]], sqrt(gee$naive.variance[2,2]),2*pnorm(abs(gee$coefficients[[2]]/sqrt(gee$naive.variance[2,2])),lower.tail = F))
    
  }
}
cbind(gee_result$V3, gee_result$V3-1.96*sqrt(gee_result$V4), gee_result$V3+1.96*sqrt(gee_result$V4))/10

outcome=c("TBV","GM")
marker=c("ptau181_tert","ptau217_tert","ttau_tert","nfl_tert","gfap_tert")
covar="type+blood_age+sex+race_ethnicity+duration+ICV"
gee_result=data.frame(outcome=c(), marker=c(), beta=c(), se=c(), p=c())
i=0
for (o in outcome){
  for (m in marker){
    f=as.formula(paste(o,paste0("as.factor(",m,")","+",covar), sep="~"))
    i=i+1
    gee_result[i,1]=o
    gee_result[i,2]=m
    gee<- gee(f,
              data = data_long,
              id = id_participant,
              family = gaussian(),
              corstr = "exchangeable")
    gee_result[i,3:5]=c(gee$coefficients[[2]], sqrt(gee$naive.variance[2,2]),2*pnorm(abs(gee$coefficients[[2]]/sqrt(gee$naive.variance[2,2])),lower.tail = F))
    
  }
}

cbind(gee_result$V3, gee_result$V3-1.96*sqrt(gee_result$V4), gee_result$V3+1.96*sqrt(gee_result$V4))/10
gee<- gee(edss_score ~ pTau181_scale+pTau217_scale+Total.Tau_scale+nfl_scale+gfap_scale + type + blood_age + sex + race_ethnicity + 
            duration, 
          data = data,
          id = id_participant,
          family = gaussian(),
          corstr = "exchangeable"
)
c(gee$coefficients[[6]], gee$coefficients[[6]]-sqrt(gee$naive.variance[6,6]),gee$coefficients[[6]]+sqrt(gee$naive.variance[6,6]))



##### longitudinal EDSS ####

outcome="edss_score"
marker=c("ptau181_tert","ptau217_tert","ttau_tert","nfl_tert","gfap_tert")
covar="type+blood_age+sex+race_ethnicity+duration"
gee_result=data.frame( marker=c(), beta=c(), se=c(), p=c())
i=0
for (m in marker){
  f=as.formula(paste(outcome,paste0("as.factor(",m,")+",covar), sep="~"))
  i=i+1
  gee_result[i,1]=m
  gee<- gee(f,
            data = data,
            id = id_participant,
            family = gaussian(),
            corstr = "exchangeable"
  )
  gee_result[i,2:4]=c(gee$coefficients[[2]], sqrt(gee$naive.variance[2,2]),2*pnorm(abs(gee$coefficients[[2]]/sqrt(gee$naive.variance[2,2])),lower.tail = F))
  
}
gee_result
cbind(gee_result$V3, gee_result$V3-1.96*sqrt(gee_result$V4), gee_result$V3+1.96*sqrt(gee_result$V4))


