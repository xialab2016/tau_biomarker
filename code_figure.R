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

pls=read.csv("~/OneDrive - University of Pittsburgh/XiaLab/ADBiomarker/plasma_CleanData_230828.csv")
pls=pls %>% filter(subtype_blood !="Related disorder")
pls=pls %>% filter(is.na(pTau217)==F)

#### Table 1 ####
pls=pls %>% mutate(race_ethnicity=ifelse(race=="White" & ethnicity=="Not Hispanic or Latino","Non-Hispanic White","Other" ))
pls=pls %>% mutate(subtype_blood1=ifelse(subtype_blood %in% c("RIS","RRMS"),"RRMS",
                                         ifelse(subtype_blood=="SPMS","SPMS",
                                                ifelse(subtype_blood=="PPMS","PPMS","Control"))))
pls$subtype_blood1=factor(pls$subtype_blood1, levels = c("Control","RRMS","PPMS","SPMS"))
pls$subtype=ifelse(pls$subtype_blood1=="Control","Control","MS")
tab=CreateTableOne(vars = c("blood_age","sex","race","ethnicity","race_ethnicity","pTau181", "pTau217","Total.Tau", "NFL","GFAP"),
                   strata = "subtype",
                   addOverall = T,
                   data=pls,
                   factorVars = c("sex","race","ethnicity","race_ethnicity"),
                   includeNA = F,
                   argsApprox = list(correct = F), #chi-square test no correction
                   argsNormal = list(var.equal = T) # 2-sample t-test with equal var
)
tab=print(tab,
          exact=c("sex","race","ethnicity","race_ethnicity"),  #exact fisher test
          nonnormal = c("blood_age","pTau181", "pTau217","Total.Tau", "NFL","GFAP"),  
          catDigits = 2,
          contDigits = 2,
          pDigits = 3,
          formatOptions = list(big.mark = ","),
          showAllLevels = TRUE)
#### sFigure 1 age, sex, and race wiith marrkers in MS ####
p1=pls %>% filter(subtype=="MS") %>% 
  ggplot(aes(x=blood_age,y=pTau181))+
  geom_point()+
  geom_smooth(method="loess")+
  labs(x="Age (years)",y="p-tau181 (pg/mL)")+
  theme_light()+
  theme(axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12))
cor.test(x=pls$blood_age[pls$subtype=="MS"],y=pls$pTau181[pls$subtype=="MS"],method="spearman") #0.27, p<0.01
p2=pls %>% filter(subtype=="MS") %>% 
  ggplot(aes(x=blood_age,y=pTau217))+
  geom_point()+
  geom_smooth(method="loess")+
  labs(x="Age (years)",y="p-tau217 (pg/mL)")+
  theme_light()+
  theme(axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12))
cor.test(x=pls$blood_age[pls$subtype=="MS"],y=pls$pTau217[pls$subtype=="MS"],method="spearman") #0.38, p<0.01

p3=pls %>% filter(subtype=="MS") %>% 
  ggplot(aes(x=blood_age,y=Total.Tau))+
  geom_point()+
  geom_smooth(method="loess")+
  labs(x="Age (years)",y="t-tau (pg/mL)")+
  theme_light()+
  theme(axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12))
cor.test(x=pls$blood_age[pls$subtype=="MS"],y=pls$Total.Tau[pls$subtype=="MS"],method="spearman") #-0.06, p=0.42

p4=pls %>% filter(subtype=="MS") %>% 
  ggplot(aes(x=blood_age,y=NFL))+
  geom_point()+
  geom_smooth(method="loess")+
  labs(x="Age (years)",y="NfL (pg/mL)")+
  theme_light()+
  theme(axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12))
cor.test(x=pls$blood_age[pls$subtype=="MS"],y=pls$NFL[pls$subtype=="MS"],method="spearman") #0.49, p<0.01

p5=pls %>% filter(subtype=="MS") %>% 
  ggplot(aes(x=blood_age,y=GFAP))+
  geom_point()+
  geom_smooth(method="loess")+
  labs(x="Age (years)",y="GFAP (pg/mL)")+
  theme_light()+
  theme(axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12))
cor.test(x=pls$blood_age[pls$subtype=="MS"],y=pls$GFAP[pls$subtype=="MS"],method="spearman") #0.39, p<0.01

P_age=p1+p2+p3+p4+p5+plot_layout(ncol=5)

pls$sex=ifelse(pls$sex=="M","Male","Female")
p1=pls %>% filter(subtype=="MS") %>% 
  ggplot(aes(x=sex,y=pTau181))+
  geom_boxplot(aes(color=sex))+
  labs(y="p-tau181 (pg/mL)",x="")+
  theme_light()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        legend.position = "none")
wilcox.test(pls$pTau181[pls$subtype=="MS"] ~ pls$sex[pls$subtype=="MS"])  #p=0.26

p2=pls %>% filter(subtype=="MS") %>% 
  ggplot(aes(x=sex,y=pTau217))+
  geom_boxplot(aes(color=sex))+
  labs(y="p-tau217 (pg/mL)",x="")+
  theme_light()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        legend.position = "none")
wilcox.test(pls$pTau217[pls$subtype=="MS"] ~ pls$sex[pls$subtype=="MS"])  #p=0.31
p3=pls %>% filter(subtype=="MS") %>% 
  ggplot(aes(x=sex,y=Total.Tau))+
  geom_boxplot(aes(color=sex))+
  labs(y="t-tau (pg/mL)",x="")+
  theme_light()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        legend.position = "none")
wilcox.test(pls$Total.Tau[pls$subtype=="MS"] ~ pls$sex[pls$subtype=="MS"])  #p<0.01

p4=pls %>% filter(subtype=="MS") %>% 
  ggplot(aes(x=sex,y=NFL))+
  geom_boxplot(aes(color=sex))+
  labs(y="NfL (pg/mL)",x="")+
  theme_light()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        legend.position = "none")
wilcox.test(pls$NFL[pls$subtype=="MS"] ~ pls$sex[pls$subtype=="MS"])  #p=0.14
p5=pls %>% filter(subtype=="MS") %>% 
  ggplot(aes(x=sex,y=GFAP))+
  geom_boxplot(aes(color=sex))+
  labs(y="GFAP (pg/mL)",x="")+
  theme_light()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        legend.title=element_blank(),
        legend.text = element_text(face="bold",size=13))
wilcox.test(pls$GFAP[pls$subtype=="MS"] ~ pls$sex[pls$subtype=="MS"])  #p=0.89

P_sex=p1+p2+p3+p4+p5+plot_layout(ncol=5)

pls$race_ethnicity=ifelse(pls$race_ethnicity=="Non-Hispanic White","Non-hispanic white","Not non-hispanic white")
pls$race_ethnicity=as.factor(pls$race_ethnicity)
p1=pls %>% filter(subtype=="MS") %>% 
  ggplot(aes(x=race_ethnicity,y=pTau181,color=race_ethnicity))+
  geom_boxplot()+
  scale_color_brewer(palette="Set1")+
  labs(x="",y="p-tau181 (pg/mL)")+
  theme_light()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        legend.position = "none")



wilcox.test(pls$pTau181[pls$subtype=="MS"] ~ pls$race_ethnicity[pls$subtype=="MS"])  #p=0.46

p2=pls %>% filter(subtype=="MS") %>% 
  ggplot(aes(x=race_ethnicity,y=pTau217))+
  geom_boxplot(aes(color=race_ethnicity))+
  scale_color_brewer(palette="Set1")+
  labs(x="",y="p-tau217 (pg/mL)")+
  theme_light()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        legend.position = "none")
wilcox.test(pls$pTau217[pls$subtype=="MS"] ~ pls$race_ethnicity[pls$subtype=="MS"])  #p=0.58
p3=pls %>% filter(subtype=="MS") %>% 
  ggplot(aes(x=race_ethnicity,y=Total.Tau))+
  geom_boxplot(aes(color=race_ethnicity))+
  scale_color_brewer(palette="Set1")+
  labs(x="",y="t-tau (pg/mL)")+
  theme_light()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        legend.position = "none")
wilcox.test(pls$Total.Tau[pls$subtype=="MS"] ~ pls$race_ethnicity[pls$subtype=="MS"])  #p=0.76

p4=pls %>% filter(subtype=="MS") %>% 
  ggplot(aes(x=race_ethnicity,y=NFL))+
  scale_color_brewer(palette="Set1")+
  geom_boxplot(aes(color=race_ethnicity))+
  labs(x="",y="NfL (pg/mL)")+
  theme_light()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        legend.position = "none")
wilcox.test(pls$NFL[pls$subtype=="MS"] ~ pls$race_ethnicity[pls$subtype=="MS"])  #p=0.42
p5=pls %>% filter(subtype=="MS") %>% 
  ggplot(aes(x=race_ethnicity,y=GFAP))+
  scale_color_brewer(palette="Set1")+
  geom_boxplot(aes(color=race_ethnicity))+
  labs(x="",y="GFAP (pg/mL)")+
  theme_light()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        legend.title=element_blank(),
        legend.text = element_text(face="bold",size=13))
        
wilcox.test(pls$GFAP[pls$subtype=="MS"] ~ pls$race_ethnicity[pls$subtype=="MS"])  #p=0.73

P_race=p1+p2+p3+p4+p5+plot_layout(ncol=5)

#### sFigure 1 correlation of markers in MS ####
tmp=pls
colnames(tmp)[c(11,12,13,16,17)]=c("p-tau181","p-tau217","t-tau","NfL","GFAP")
ggpairs(tmp[,c(11,12,13,16,17)],
        upper = list(continuous = wrap("cor", method = "spearman")))
#### Figure 2 compared MS and subtype vvs. controls ####
tmp<- pls %>% mutate(lpTau181=log(pTau181),lpTau217=log(pTau217),
                     lTotal.Tau=log(Total.Tau),lGFAP=log(Total.Tau),
                     lNFL=log(NFL))

m=lm(lpTau181 ~subtype+ blood_age+sex,data=tmp); summary(m) #beta=0.12, p=0.38
tmp=tmp%>% mutate(pred_ptau181=exp(predict(m)))
p1=ggplot(tmp, aes(x=subtype, y=pred_ptau181,fill=subtype)) + 
  geom_violin(trim=FALSE)+
  scale_fill_brewer(palette="Reds") + 
  labs(x="", y = "Adjusted p-tau181 (pg/mL)")+
  geom_boxplot(width=0.1,fill="white")+
  theme_bw()+
  ylim(c(1.5,5.1))+
  theme(axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        legend.position = "none")
wilcox.test(tmp$pred_ptau181 ~ tmp$subtype) #0.000884


m=lm(lpTau217 ~subtype+ blood_age+sex,data=tmp); summary(m) #beta=0.15, p=0.27
tmp=tmp%>% mutate(pred_ptau217=exp(predict(m)))
p2=ggplot(tmp, aes(x=subtype, y=pred_ptau217,fill=subtype)) + 
  geom_violin(trim=FALSE)+
  scale_fill_brewer(palette="Greens") + 
  labs(x="", y = "Adjusted p-tau217 (pg/mL)")+
  geom_boxplot(width=0.1,fill="white")+
  theme_bw()+
  ylim(c(0.39, 1.7))+
  theme(axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        legend.position = "none")
wilcox.test(tmp$pred_ptau217 ~ tmp$subtype) #0.0001984


m=lm(lTotal.Tau ~subtype+ blood_age+sex,data=tmp); summary(m) #beta=-0.04, p=0.82
tmp=tmp%>% mutate(pred_ttau=exp(predict(m)))
p3=ggplot(tmp, aes(x=subtype, y=pred_ttau,fill=subtype)) + 
  geom_violin(trim=FALSE)+
  scale_fill_brewer(palette="Greys") + 
  labs(x="", y = "Adjusted t-tau (pg/mL)")+
  geom_boxplot(width=0.1,fill="white")+
  theme_bw()+
  ylim(c(1.5, 3.4))+
  theme(axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        legend.position = "none")
wilcox.test(tmp$pred_ttau ~ tmp$subtype) #0.11

tmp<-tmp %>% mutate(subtype2=ifelse(subtype_blood1%in% c("PPMS","SPMS"),"PMS",ifelse(subtype_blood1=="RRMS","RRMS","Control")),
                    subtype2=factor(subtype2, levels = c("Control","RRMS","PMS")))
m=lm(lpTau181 ~subtype2+ blood_age+sex,data=tmp); summary(m) #beta=0.12, p=0.38
anova(m)#p=0.002469 
l=matrix(c(0,1,0,0,0,
           0,0,1,0,0,
           0,1,-1,0,0),byrow = T,nrow = 3)
summary(multcomp::glht(m, l))

tmp=tmp%>% mutate(pred_ptau181=exp(predict(m)))
p4=ggplot(tmp, aes(x=subtype2, y=pred_ptau181,fill=subtype2)) + 
  geom_violin(trim=FALSE)+
  scale_fill_brewer(palette="Reds") + 
  labs(x="", y = "Adjusted p-tau181 (pg/mL)")+
  geom_boxplot(width=0.1,fill="white")+
  theme_bw()+
  ylim(c(1.8,5.3))+
  theme(axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        legend.position = "none")

kruskal.test(tmp$pred_ptau181 ~tmp$subtype2) #p-value < 2.2e-16
dunn.test(tmp$pred_ptau181,tmp$subtype2,kw=T,method = "bh")


m=lm(lpTau217 ~subtype2+ blood_age+sex,data=tmp); summary(m) 
anova(m)#p=0.0015496
tmp=tmp%>% mutate(pred_ptau217=exp(predict(m)))
p5=ggplot(tmp, aes(x=subtype2, y=pred_ptau217,fill=subtype2)) + 
  geom_violin(trim=FALSE)+
  scale_fill_brewer(palette="Greens") + 
  labs(x="", y = "Adjusted p-tau217 (pg/mL)")+
  geom_boxplot(width=0.1,fill="white")+
  theme_bw()+
  ylim(c(0.4, 1.8))+
  theme(axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        legend.position = "none")

m=lm(lTotal.Tau ~subtype2+ blood_age+sex,data=tmp); summary(m) 
anova(m)#p=0.6009
tmp=tmp%>% mutate(pred_ttau=exp(predict(m)))

p6=ggplot(tmp, aes(x=subtype2, y=pred_ttau,fill=subtype2)) + 
  geom_violin(trim=FALSE)+
  scale_fill_brewer(palette="Greys") + 
  labs(x="", y = "Adjusted t-tau (pg/mL)")+
  geom_boxplot(width=0.1,fill="white")+
  theme_bw()+
  ylim(c(1.4, 3.6))+
  theme(axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        legend.position = "none")
tiff("Fig2.tiff", units="in", width=10, height=8, res=300)
p1+p2+p3+p4+p5+p6+plot_layout(ncol = 3)
dev.off()

#### Figure 3A ####
pls_PMS=pls %>% filter(subtype_blood1 %in% c("SPMS","PPMS","RRMS")) %>% 
  mutate(outcome=ifelse(subtype_blood1 =="RRMS",0,1))
pls_PMS $duration=pls_PMS$blood_age -pls_PMS$age_ms_dx
pls_PMS $duration[which(pls_PMS$duration<0)]=0
summary(pls_PMS$duration)
pls_PMS=pls_PMS %>% mutate(spTau181=as.numeric(scale(pTau181)),
                           spTau217=as.numeric(scale(pTau217)),
                           sTotal.Tau=as.numeric(scale(Total.Tau)),
                           sGFAP=as.numeric(scale(GFAP)),
                           sNFL=as.numeric(scale(NFL)))

pls_PMS=pls %>% filter(subtype_blood1 %in% c("SPMS","PPMS","RRMS")) %>% 
  mutate(outcome=ifelse(subtype_blood1 =="RRMS",0,1))
pls_PMS $duration=pls_PMS$blood_age -pls_PMS$age_ms_dx
pls_PMS $duration[which(pls_PMS$duration<0)]=0
summary(pls_PMS$duration)
pls_PMS=pls_PMS %>% mutate(spTau181=as.numeric(scale(pTau181)),
                           spTau217=as.numeric(scale(pTau217)),
                           sTotal.Tau=as.numeric(scale(Total.Tau)),
                           sGFAP=as.numeric(scale(GFAP)),
                           sNFL=as.numeric(scale(NFL)))
pls_PMS=pls_PMS %>% mutate(duration_c=ifelse(duration <=5,"<=5 yrs", ifelse(duration<=15,"5-15 yrs",">15 yrs")),
                           age_c=ifelse(blood_age<=45,"<=45 yrs",ifelse(blood_age<=60,"45-60 yrs",">60 yrs")),
                           duration_c=factor(duration_c,levels = c("<=5 yrs","5-15 yrs",">15 yrs")),
                           age_c=factor(age_c,levels = c("<=45 yrs","45-60 yrs",">60 yrs")))
match_obj <- matchit(outcome ~ sex + age_c+race_ethnicity+duration_c,
                     data = pls_PMS, method="full",distance = "glm",
                     ratio =1,
                     replace = F)
datanew=match.data(match_obj)
table(datanew$outcome)
m=glm(outcome ~duration+I(duration^2)+GFAP+sNFL+spTau181+spTau217+sTotal.Tau,data=datanew,  family = binomial())
summary(m)

matched=data.frame(marker=rep(c("p-tau181","p-tau217","t-tau","NfL","GFAP"),2),
                   model=c(rep("Model 1: Separately",5),rep("Model 2: Simultaneously",5)),
                   or=c(2.30,3.03, 1.03, 1.10, 1.05,
                        exp(0.3097051),exp(0.5031172),exp( -0.1033199),exp(0.0896585),exp(0.0031262)),
                   low=c( 1.40,1.75, 0.65, 0.74, 0.66,
                          exp(0.3097051-1.96*0.1406883),exp(0.5031172-1.96*0.2352299),exp( -0.1033199-1.96*0.2193829),exp(0.0896585-1.96*0.2943977),exp(0.0031262-1.96*0.0041077  )),
                   high=c( 4.13, 5.67,1.61, 1.64, 1.68,
                           exp(0.3097051+1.96*0.1406883),exp(0.5031172+1.96*0.2352299),exp( -0.1033199+1.96*0.2193829),exp(0.0896585+1.96*0.2943977),exp(0.0031262+1.96*0.0041077 )))
matched$marker=factor(matched$marker,levels = c("NfL","GFAP","p-tau181","p-tau217","t-tau"))
tiff("Fig3a.tiff", units="in", width=9, height=6, res=300)
matched %>% ggplot(aes(x=marker,y=or,fill=marker))+
  geom_hline(yintercept = 1,col=2,linetype="dashed",linewidth=0.8)+
  geom_errorbar(aes(ymin=low, ymax=high), width=.1,linewidth=1) +
  geom_point(aes(x=marker,y=or),size=2)+
  ylim(c(0,6))+
  labs(x="",y="OR")+
  theme_bw()+
  facet_wrap(~model,ncol = 1)+
  theme(legend.position = "none",
        axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        strip.text = element_text(size = 12, face="bold")
)

dev.off()
#### Figure 3B ####
m_0=glm(outcome ~race_ethnicity+duration,data=pls_PMS,  family = binomial())
rocobj=roc(pls_PMS$outcome==1, predict(m_0,type="response"),ci = T)
auc(rocobj) #auc=0.6372
m_1=glm(outcome ~blood_age+race_ethnicity+duration+sGFAP+sNFL,data=pls_PMS,  family = binomial())
rocobj=roc(pls_PMS$outcome==1, predict(m_1,type="response"),ci = T)
auc(rocobj) #auc=0.8189

m_2=glm(outcome ~blood_age+race_ethnicity+duration+spTau181+spTau217+sTotal.Tau,data=pls_PMS,  family = binomial())
rocobj=roc(pls_PMS$outcome==1, predict(m_2,type="response"),ci = T)
auc(rocobj) #auc=0.8275

m_3=glm(outcome ~blood_age+sex+race_ethnicity+duration+sGFAP+sNFL+spTau181+spTau217+sTotal.Tau,data=pls_PMS,  family = binomial())
rocobj=roc(pls_PMS$outcome==1, predict(m_3,type="response"),ci = T)
auc(rocobj) #0.8912

rocs <- list()
rocs[["Clinical features only (AUC=0.64)"]] <- roc(pls_PMS$outcome, predict(m_0),ci = T)
rocs[["Clinical features+NfL+GFAP (AUC=0.82)"]] <- roc(pls_PMS$outcome, predict(m_1),ci=T)
rocs[["Clinical features+Tau (AUC=0.83)"]] <- roc(pls_PMS$outcome, predict(m_2),ci=T)
rocs[["Clinical features+NfL+GFAP+Tau (AUC=0.89)"]] <- roc(pls_PMS$outcome, predict(m_3),ci=T)

cis=list()

cc <- ci.se(rocs$`Clinical features only (AUC=0.64)`, specificities=seq(0, 1, l=25))
cis[["Clinical features only (AUC=0.64)"]] <- data.frame(x = as.numeric(rownames(cc)),
                                                         upper = cc[, 1],
                                                         lower=  cc[, 3])

cc <- ci.se(rocs$`Clinical features+NfL+GFAP (AUC=0.82)`, specificities=seq(0, 1, l=25))
cis[["Clinical features+NfL+GFAP (AUC=0.82)"]] <- data.frame(x = as.numeric(rownames(cc)),
                                                             upper = cc[, 1],
                                                             lower=  cc[, 3])

cc <- ci.se(rocs$`Clinical features+Tau (AUC=0.83)`, specificities=seq(0, 1, l=25))
cis[["Clinical features+Tau (AUC=0.83)"]] <- data.frame(x = as.numeric(rownames(cc)),
                                                             upper = cc[, 1],
                                                             lower=  cc[, 3])
cc <- ci.se(rocs$`Clinical features+NfL+GFAP+Tau (AUC=0.89)`, specificities=seq(0, 1, l=25))
cis[["Clinical features+NfL+GFAP+Tau (AUC=0.89)"]] <- data.frame(x = as.numeric(rownames(cc)),
                                                        upper = cc[, 1],
                                                        lower=  cc[, 3])


pp=ggroc(rocs)+
  theme_bw()+
  theme(axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        legend.position = c(0.75, 0.39),
        legend.title = element_blank(),
        legend.text = element_text(size=11, face="bold"))
pp

p=pp+geom_ribbon(
    data = cis[[1]],
    aes(x = x, ymin = lower, ymax = upper),
    fill = 1 + 1,
    alpha = 0.2,
    inherit.aes = F) 
p=p+geom_ribbon(
  data = cis[[2]],
  aes(x = x, ymin = lower, ymax = upper),
  fill = 2 + 1,
  alpha = 0.2,
  inherit.aes = F) 
p=p+geom_ribbon(
  data = cis[[3]],
  aes(x = x, ymin = lower, ymax = upper),
  fill = 3 + 1,
  alpha = 0.2,
  inherit.aes = F)
p=p+geom_ribbon(
  data = cis[[4]],
  aes(x = x, ymin = lower, ymax = upper),
  fill = "plum2",
  alpha = 0.2,
  inherit.aes = F) 
p

net_benefit <- function(actual,prediction,thresholds){
  N <- length(actual)
  sapply(thresholds,function(threshold){
    prediction_thresholded <- as.integer(prediction > threshold)
    tpr <- sum(prediction_thresholded*actual) / N
    fpr <- sum(prediction_thresholded*(1-actual)) / N
    tpr - fpr * threshold/(1-threshold)
  })
}

actual=pls_PMS$outcome
thresholds=seq(0,1,l=100)
nb_0=net_benefit(actual = actual, prediction = predict(m_0,type = "response"), thresholds)
nb_1=net_benefit(actual = actual, prediction = predict(m_1,type = "response"), thresholds)
nb_2=net_benefit(actual = actual, prediction = predict(m_2,type = "response"), thresholds)
nb_3=net_benefit(actual = actual, prediction = predict(m_3,type = "response"), thresholds)
nb=data.frame(thresholds=rep(thresholds,4), nb=c(nb_0, nb_1, nb_2, nb_3),model=c(rep("Clinical features only",100),rep("Clinical features only+NfL+GFAP",100), rep("Clinical features only+Tau",100), rep("Clinical features only+NfL+GFAP+Tau",100) ))
nb$model=factor(nb$model, levels = c("Clinical features only","Clinical features only+NfL+GFAP", "Clinical features only+Tau","Clinical features only+NfL+GFAP+Tau"))
nb$nb[nb$nb<0| is.na(nb$nb)==T]=0

ggplot(data=nb, aes(x=thresholds, y=nb, group=model)) +
  geom_smooth(method = "loess",aes(color=model), se = F)+
  labs(x="Threshold probability", y="Net benefit")+ 
  theme_bw()+
  theme(axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        legend.position = c(0.75, 0.49),
        legend.title = element_blank(),
        legend.text = element_text(size=11, face="bold"))
  

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

model0 <- train(outcome ~ race_ethnicity+duration, data=train_df,
                method="glm", 
                family=binomial, 
                trControl=ctrlspecs)
predictions <- predict(model0, newdata=test_df)
confusionMatrix(data=predictions, test_df$outcome,positive = "PMS")

model1 <- train(outcome ~ blood_age + race_ethnicity + duration +   sGFAP + sNFL,
                data=train_df,
                method="glm", 
                family=binomial, 
                trControl=ctrlspecs)
predictions <- predict(model1, newdata=test_df)
confusionMatrix(data=predictions, test_df$outcome,positive = "PMS")

model2 <- train(outcome ~ blood_age + race_ethnicity + duration + 
                  spTau181 + spTau217 + sTotal.Tau,
                data=train_df,
                method="glm", 
                family=binomial, 
                trControl=ctrlspecs)
predictions <- predict(model2, newdata=test_df)
confusionMatrix(data=predictions, test_df$outcome,positive = "PMS")

model3 <- train(outcome ~ blood_age + sex + race_ethnicity + duration + 
                  sGFAP + sNFL + spTau181 + spTau217 + sTotal.Tau,
                data=train_df,
                method="glm", 
                family=binomial, 
                trControl=ctrlspecs)
predictions <- predict(model3, newdata=test_df)
confusionMatrix(data=predictions, test_df$outcome,positive = "PMS")



  #### blood data ####
bld=read.csv("~/OneDrive - University of Pittsburgh/XiaLab/PROMOTEDATA/Clean_20231206//BloodSample.csv")
bld=bld %>% dplyr::select(id_participant,blood1_date,height_feet_blood,height_inches_blood,weight_blood,
                          blood1_dmt,blood1_dmt_detail,blood1_pddsv2,blood1_relapse_number)
bld<- bld %>% mutate(bmi=weight_blood/(height_feet_blood*12+height_inches_blood)^2*703,
                     bmi=ifelse(bmi<10 ,NA, bmi))

bld=bld %>% group_by(id_participant) %>% mutate(bmi_avg=mean(bmi, na.rm=T))
bld$bmi[which(is.na(bld$bmi)==T)]=bld$bmi_avg[which(is.na(bld$bmi)==T)]
bld<- bld %>% mutate(obese=ifelse(bmi>=30, 1,0))
table(bld$blood1_dmt_detail,useNA = "ifany")
bld=bld %>% mutate(dmt=ifelse(blood1_dmt==0 | blood1_dmt_detail==99| blood1_dmt_detail==11, "None", 
                              ifelse(blood1_dmt_detail %in% c(1,2,3,4,6,8,9,14,19),"Standard","High")))
table(bld$blood1_dmt,bld$dmt,useNA = "ifany") 
#high:Tysabri (Natalizumab),Ocrevus (Ocrelizumab),Rituxan (Rituximab),Mavenclad (Cladribine),Kesimpta (Ofatumumab)
#standard: Aubagio (Teriflunomide),	Avonex (Interferon Beta-1a),Betaseron (Interferon Beta-1b),Copaxone/Glatopa (Glatiramer acetate),
#Gilenya (Fingolimod),Rebif (Interferon Beta-1a),Tecfidera (Dimethyl fumarate/BG12),Plegridy (Pegylated interferon Beta-1a)
bld=bld %>% mutate(relapse_1y=ifelse(blood1_relapse_number ==99,0,1))
bld=bld %>% dplyr::select(id_participant,blood1_date,bmi, obese,dmt, relapse_1y,blood1_pddsv2)
colnames(bld)[2]="blood_date"
pls=read.csv("~/OneDrive - University of Pittsburgh/XiaLab/ADBiomarker/plasma_CleanData_230828.csv")
pls=pls %>% filter(subtype_blood !="Related disorder")
pls=pls %>% filter(is.na(pTau217)==F)
pls=pls %>% mutate(race_ethnicity=ifelse(race=="White" & ethnicity=="Not Hispanic or Latino","Non-Hispanic White","Other" ))
pls=pls %>% mutate(subtype_blood1=ifelse(subtype_blood %in% c("RIS","RRMS"),"RRMS",
                                         ifelse(subtype_blood=="SPMS","SPMS",
                                                ifelse(subtype_blood=="PPMS","PPMS","Control"))))
pls$subtype_blood1=factor(pls$subtype_blood1, levels = c("Control","RRMS","PPMS","SPMS"))
pls$subtype=ifelse(pls$subtype_blood1=="Control","Control","MS")
pls_MS=pls %>% filter(subtype=="MS")
pls_MS=pls_MS %>% left_join(bld, by=c("id_participant","blood_date"))
pls_MS=pls_MS [-128,]
pls_MS=pls_MS %>% mutate(duration=blood_age-age_sx_onset,
                         subtype=ifelse(subtype_blood%in%c("RIS","RRMS"),"RRMS","PMS"))
pls_MS=pls_MS %>% dplyr::select(id_participant,blood_date,blood_age,sex,race_ethnicity,subtype,duration,
                         bmi, obese,dmt, relapse_1y,blood1_pddsv2.y,pTau181,pTau217,Total.Tau,GFAP,NFL)
colnames(pls_MS)[12]="pdds_bl"
write.csv(pls_MS,"~/OneDrive - University of Pittsburgh/XiaLab/ADBiomarker/plasma_CleanData_MS_231206.csv")

#### Forest plot 3 months####
result=readxl::read_xlsx("~/OneDrive - University of Pittsburgh/XiaLab/ADBiomarker/manuscript/20240903//association_3month.xlsx")
result$Model=factor(result$Model, levels = c("Separate","Simultaneous"))
result$marker=factor(result$marker, levels = c("NfL","GFAP","p-tau181","p-tau217","t-tau"))
p1=result %>% filter(outcome=="PDDS") %>%
  ggplot(aes(x=marker,y=est))+
  geom_hline(yintercept = 0,col=2,linetype="dashed",linewidth=0.8)+
  geom_errorbar(aes(ymin=low, ymax=up,color=Model), width=.1,linewidth=1, position = position_dodge(0.6)) +
  geom_point(aes(x=marker,y=est, shape=Model,color=Model),size=2.5, position = position_dodge(0.6))+
  labs(x="",y="Beta (95%CI)")+
  ggtitle("A. PDDS (n=160)")+
  scale_color_manual(values=c('darkorange','darkblue'))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        plot.title = element_text(size=14, face="bold")
       
)

p2=result %>% filter(outcome=="Normalized age-related MS severity") %>%
  ggplot(aes(x=marker,y=est))+
  geom_hline(yintercept = 0,col=2,linetype="dashed",linewidth=0.8)+
  geom_errorbar(aes(ymin=low, ymax=up,color=Model), width=.1,linewidth=1, position = position_dodge(0.6)) +
  geom_point(aes(x=marker,y=est, shape=Model,color=Model),size=2.5, position = position_dodge(0.6))+
  labs(x="",y="Beta (95%CI)")+
  ggtitle("B. Normalized age-related MS severity (n=83)")+
  scale_color_manual(values=c('darkorange','darkblue'))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        plot.title = element_text(size=14, face="bold")
        
  )
p3=result %>% filter(outcome=="t25fw") %>%
  ggplot(aes(x=marker,y=est))+
  geom_hline(yintercept = 0,col=2,linetype="dashed",linewidth=0.8)+
  geom_errorbar(aes(ymin=low, ymax=up,color=Model), width=.1,linewidth=1, position = position_dodge(0.6)) +
  geom_point(aes(x=marker,y=est, shape=Model,color=Model),size=2.5, position = position_dodge(0.6))+
  labs(x="",y="Beta (95%CI)")+
  ggtitle("C. Timed 25-foot walk (n=103)")+
  scale_color_manual(values=c('darkorange','darkblue'))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        plot.title = element_text(size=14, face="bold")
        
  )

p4=result %>% filter(outcome=="9peg") %>%
  ggplot(aes(x=marker,y=est))+
  geom_hline(yintercept = 0,col=2,linetype="dashed",linewidth=0.8)+
  geom_errorbar(aes(ymin=low, ymax=up,color=Model), width=.1,linewidth=1, position = position_dodge(0.6)) +
  geom_point(aes(x=marker,y=est, shape=Model,color=Model),size=2.5, position = position_dodge(0.6))+
  labs(x="",y="Beta (95%CI)")+
  ggtitle("D. Nine-hole peg (n=99)")+
  scale_color_manual(values=c('darkorange','darkblue'))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        plot.title = element_text(size=14, face="bold")
        
  )

p5=result %>% filter(outcome=="SDMT") %>%
  ggplot(aes(x=marker,y=est))+
  geom_hline(yintercept = 0,col=2,linetype="dashed",linewidth=0.8)+
  geom_errorbar(aes(ymin=low, ymax=up,color=Model), width=.1,linewidth=1, position = position_dodge(0.6)) +
  geom_point(aes(x=marker,y=est, shape=Model,color=Model),size=2.5, position = position_dodge(0.6))+
  labs(x="",y="Beta (95%CI)")+
  ggtitle("E. Symbol digit modalities test (n=95)")+
  scale_color_manual(values=c('darkorange','darkblue'))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        plot.title = element_text(size=14, face="bold")
        
  )


p6=result %>% filter(outcome=="RNFL") %>%
  ggplot(aes(x=marker,y=est))+
  geom_hline(yintercept = 0,col=2,linetype="dashed",linewidth=0.8)+
  geom_errorbar(aes(ymin=low, ymax=up,color=Model), width=.1,linewidth=1, position = position_dodge(0.6)) +
  geom_point(aes(x=marker,y=est, shape=Model,color=Model),size=2.5, position = position_dodge(0.6))+
  labs(x="",y="Beta (95%CI)")+
  ggtitle("F. Retinal nerve fiber layer thickness (n=84)")+
  scale_color_manual(values=c('darkorange','darkblue'))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        plot.title = element_text(size=14, face="bold")
        
  )


p7=result %>% filter(outcome=="Percentage of total brain volume loss") %>%
  ggplot(aes(x=marker,y=est))+
  geom_hline(yintercept = 0,col=2,linetype="dashed",linewidth=0.8)+
  geom_errorbar(aes(ymin=low, ymax=up,color=Model), width=.1,linewidth=1, position = position_dodge(0.6)) +
  geom_point(aes(x=marker,y=est, shape=Model,color=Model),size=2.5, position = position_dodge(0.6))+
  labs(x="",y="Beta (95%CI)")+
  ggtitle("G. Percentage of total brain volume loss (n=69) ")+
  scale_color_manual(values=c('darkorange','darkblue'))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        plot.title = element_text(size=14, face="bold")
        
  )


p8=result %>% filter(outcome=="Percentage of gray matter volume loss") %>%
  ggplot(aes(x=marker,y=est))+
  geom_hline(yintercept = 0,col=2,linetype="dashed",linewidth=0.8)+
  geom_errorbar(aes(ymin=low, ymax=up,color=Model), width=.1,linewidth=1, position = position_dodge(0.6)) +
  geom_point(aes(x=marker,y=est, shape=Model,color=Model),size=2.5, position = position_dodge(0.6))+
  labs(x="",y="Beta (95%CI)")+
  ggtitle("H. Percentage of gray matter volume loss (n=69) ")+
  scale_color_manual(values=c('darkorange','darkblue'))+
  theme_bw()+
  theme(
        axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        plot.title = element_text(size=14, face="bold")
        
  )
tiff("~/OneDrive - University of Pittsburgh/XiaLab/ADBiomarker/manuscript/20240903/Assoc_3month.tiff", units="in", width=12, height=12, res=300)
p1+p2+p3+p4+p5+p6+p7+p8+plot_layout(ncol=2)
dev.off()

#### Forest plot 6 months####
result=readxl::read_xlsx("~/OneDrive - University of Pittsburgh/XiaLab/ADBiomarker/manuscript/20240903//association_6month.xlsx")
result$Model=factor(result$Model, levels = c("Separate","Simultaneous"))
result$marker=factor(result$marker, levels = c("NfL","GFAP","p-tau181","p-tau217","t-tau"))
p1=result %>% filter(outcome=="PDDS") %>%
  ggplot(aes(x=marker,y=est))+
  geom_hline(yintercept = 0,col=2,linetype="dashed",linewidth=0.8)+
  geom_errorbar(aes(ymin=low, ymax=up,color=Model), width=.1,linewidth=1, position = position_dodge(0.6)) +
  geom_point(aes(x=marker,y=est, shape=Model,color=Model),size=2.5, position = position_dodge(0.6))+
  labs(x="",y="Beta (95%CI)")+
  ggtitle("A. PDDS (n=142)")+
  scale_color_manual(values=c('darkorange','darkblue'))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        plot.title = element_text(size=14, face="bold")
        
  )

p2=result %>% filter(outcome=="Normalized age-related MS severity") %>%
  ggplot(aes(x=marker,y=est))+
  geom_hline(yintercept = 0,col=2,linetype="dashed",linewidth=0.8)+
  geom_errorbar(aes(ymin=low, ymax=up,color=Model), width=.1,linewidth=1, position = position_dodge(0.6)) +
  geom_point(aes(x=marker,y=est, shape=Model,color=Model),size=2.5, position = position_dodge(0.6))+
  labs(x="",y="Beta (95%CI)")+
  ggtitle("B. Normalized age-related MS severity (n=66)")+
  scale_color_manual(values=c('darkorange','darkblue'))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        plot.title = element_text(size=14, face="bold")
        
  )
p3=result %>% filter(outcome=="t25fw") %>%
  ggplot(aes(x=marker,y=est))+
  geom_hline(yintercept = 0,col=2,linetype="dashed",linewidth=0.8)+
  geom_errorbar(aes(ymin=low, ymax=up,color=Model), width=.1,linewidth=1, position = position_dodge(0.6)) +
  geom_point(aes(x=marker,y=est, shape=Model,color=Model),size=2.5, position = position_dodge(0.6))+
  labs(x="",y="Beta (95%CI)")+
  ggtitle("C. Timed 25-foot walk (n=84)")+
  scale_color_manual(values=c('darkorange','darkblue'))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        plot.title = element_text(size=14, face="bold")
        
  )

p4=result %>% filter(outcome=="9peg") %>%
  ggplot(aes(x=marker,y=est))+
  geom_hline(yintercept = 0,col=2,linetype="dashed",linewidth=0.8)+
  geom_errorbar(aes(ymin=low, ymax=up,color=Model), width=.1,linewidth=1, position = position_dodge(0.6)) +
  geom_point(aes(x=marker,y=est, shape=Model,color=Model),size=2.5, position = position_dodge(0.6))+
  labs(x="",y="Beta (95%CI)")+
  ggtitle("D. Nine-hole peg (n=82)")+
  scale_color_manual(values=c('darkorange','darkblue'))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        plot.title = element_text(size=14, face="bold")
        
  )

p5=result %>% filter(outcome=="SDMT") %>%
  ggplot(aes(x=marker,y=est))+
  geom_hline(yintercept = 0,col=2,linetype="dashed",linewidth=0.8)+
  geom_errorbar(aes(ymin=low, ymax=up,color=Model), width=.1,linewidth=1, position = position_dodge(0.6)) +
  geom_point(aes(x=marker,y=est, shape=Model,color=Model),size=2.5, position = position_dodge(0.6))+
  labs(x="",y="Beta (95%CI)")+
  ggtitle("E. Symbol digit modalities test (n=77)")+
  scale_color_manual(values=c('darkorange','darkblue'))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        plot.title = element_text(size=14, face="bold")
        
  )


p6=result %>% filter(outcome=="RNFL") %>%
  ggplot(aes(x=marker,y=est))+
  geom_hline(yintercept = 0,col=2,linetype="dashed",linewidth=0.8)+
  geom_errorbar(aes(ymin=low, ymax=up,color=Model), width=.1,linewidth=1, position = position_dodge(0.6)) +
  geom_point(aes(x=marker,y=est, shape=Model,color=Model),size=2.5, position = position_dodge(0.6))+
  labs(x="",y="Beta (95%CI)")+
  ggtitle("F. Retinal nerve fiber layer thickness (n=71)")+
  scale_color_manual(values=c('darkorange','darkblue'))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        plot.title = element_text(size=14, face="bold")
        
  )


p7=result %>% filter(outcome=="Percentage of total brain volume loss") %>%
  ggplot(aes(x=marker,y=est))+
  geom_hline(yintercept = 0,col=2,linetype="dashed",linewidth=0.8)+
  geom_errorbar(aes(ymin=low, ymax=up,color=Model), width=.1,linewidth=1, position = position_dodge(0.6)) +
  geom_point(aes(x=marker,y=est, shape=Model,color=Model),size=2.5, position = position_dodge(0.6))+
  labs(x="",y="Beta (95%CI)")+
  ggtitle("G. Percentage of total brain volume loss (n=57) ")+
  scale_color_manual(values=c('darkorange','darkblue'))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(face="bold",size=11),
        axis.text.y = element_text(face="bold",size=11),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        plot.title = element_text(size=14, face="bold")
        
  )


p8=result %>% filter(outcome=="Percentage of gray matter volume loss") %>%
  ggplot(aes(x=marker,y=est))+
  geom_hline(yintercept = 0,col=2,linetype="dashed",linewidth=0.8)+
  geom_errorbar(aes(ymin=low, ymax=up,color=Model), width=.1,linewidth=1, position = position_dodge(0.6)) +
  geom_point(aes(x=marker,y=est, shape=Model,color=Model),size=2.5, position = position_dodge(0.6))+
  labs(x="",y="Beta (95%CI)")+
  ggtitle("H. Percentage of gray matter volume loss (n=57) ")+
  scale_color_manual(values=c('darkorange','darkblue'))+
  theme_bw()+
  theme(
    axis.text.x = element_text(face="bold",size=11),
    axis.text.y = element_text(face="bold",size=11),
    axis.title.x = element_text(face="bold",size=12),
    axis.title.y = element_text(face="bold",size=12),
    plot.title = element_text(size=14, face="bold")
    
  )
tiff("~/OneDrive - University of Pittsburgh/XiaLab/ADBiomarker/manuscript/20240903/Assoc_6month.tiff", units="in", width=12, height=12, res=300)
p1+p2+p3+p4+p5+p6+p7+p8+plot_layout(ncol=2)
dev.off()
