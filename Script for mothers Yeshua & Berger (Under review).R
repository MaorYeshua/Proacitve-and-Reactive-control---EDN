library(tidyverse)
library(car)
library(MKinfer)
library(lmerTest)
library(performance)
library(parameters)
library(ggeffects)
library(patchwork)
library(emmeans)

HLM_PMC_M <- read.csv("HLM_PMC_M.csv") 
parameters_M <- read.csv("parameters_psM.csv")

HLM_PMC_M$delayed_disinhibition <- factor(HLM_PMC_M$delayed_disinhibition, levels = c("NDD","DD"))
HLM_PMC_M <- HLM_PMC_M %>% mutate(delayed_disinhibition = ifelse(delayed_disinhibition == "NDD","NPI","PI"))

HLM_PMC_M$Sex_C <- factor(HLM_PMC_M$Sex_C, levels = c("Male","Female"))


SES <- parameters_M %>% select("Education_Years_1","Income","Home_2","Home_3")
SES <- as.data.frame(apply(SES, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x)))


factor_model <- factanal(
  SES,
  factors = 1,       
  method = "ml",  
  rotation = "varimax",
  scores = "regression")
print(factor_model) 
SES$SES <- factor_model$scores
parameters_M$SES <- factor_model$scores 

parameters_M <- parameters_M %>% select("ID","SES","Sex_C","Child_Age","Mother_Age","Raven":"RT.diff.PES")


TrialsProportion <- HLM_PMC_M %>%
  mutate(DD_dummy = delayed_disinhibition == "PI",
         PES_dummy = PES == "PE") %>%
  group_by(ID) %>%
  mutate(PI_prop = mean(DD_dummy),
         PE_prop = mean(PES_dummy)) %>%
  select("ID","PI_prop","PE_prop") %>%
  distinct() %>%
  ungroup() 


parameters_M <- merge(TrialsProportion,parameters_M) %>% select("SES","Sex_C","Child_Age","Mother_Age","Raven",
                                                                "PI_prop","PE_prop","mean_RT":"RT.diff.PES")

### Descriptives

parameters_M %>% summarise(n = n(),
                           Age = mean(Mother_Age),
                           SD = sd(Mother_Age),
                           min = min(Mother_Age),
                           max = max(Mother_Age))

parameters_M %>% group_by(Sex_C) %>% summarise(n = n(),
                                               Age = mean(Mother_Age),
                                               SD = sd(Mother_Age),
                                               min = min(Mother_Age),
                                               max = max(Mother_Age))

HLM_PMC_M %>% group_by(ID) %>%
  summarise(n = n()) %>% 
  ungroup() %>%
  summarise(mean = mean(n),
            sd = sd(n),
            min = min(n),
            max = max(n))

### Reliability test

#######################################################

# Delayed Disinhibition

DDRe <- HLM_PMC_M %>%
  group_by(ID,delayed_disinhibition) %>%
  mutate(Trial = 1:n()) 

DDEVEN <- DDRe[DDRe$Trial %% 2 == 0, ]
DDODD <- DDRe[DDRe$Trial %% 2 == 1, ]

DDEVEN <- DDEVEN %>%
  group_by(ID,delayed_disinhibition) %>%
  mutate(RT.DDEVEN = mean(RT)) %>%
  select("ID", "delayed_disinhibition", "RT.DDEVEN") %>%
  ungroup() 

DDEVEN <-  DDEVEN %>%
  group_by(ID) %>%
  distinct() %>%
  ungroup()

DDEVEN$ID <- as.character(DDEVEN$ID)
DDEVEN <-  DDEVEN %>%
  pivot_wider(
    names_from = delayed_disinhibition,
    values_from = RT.DDEVEN
  ) 

DDEVEN <-  DDEVEN %>%
  mutate(RT.diff.DD.EVEN = PI - NPI)

DDODD <- DDODD %>%
  group_by(ID,delayed_disinhibition) %>%
  mutate(RT.DDODD = mean(RT)) %>%
  select("ID", "delayed_disinhibition", "RT.DDODD") %>%
  ungroup() 

DDODD <-  DDODD %>%
  group_by(ID) %>%
  distinct() %>%
  ungroup()

DDODD$ID <- as.character(DDODD$ID)
DDODD <-  DDODD %>%
  pivot_wider(
    names_from = delayed_disinhibition,
    values_from = RT.DDODD
  ) 

DDODD <-  DDODD %>%
  mutate(RT.diff.DD.ODD = PI - NPI)

delayed_disinhibition <- merge(DDODD,DDEVEN,by="ID")


cor.test(delayed_disinhibition$RT.diff.DD.ODD,delayed_disinhibition$RT.diff.DD.EVEN)

#######################################################

# Post Error Slowing

PESRe <- HLM_PMC_M %>%
  group_by(ID,PES) %>%
  mutate(Trial = 1:n()) 

PESEVEN <- PESRe[PESRe$Trial %% 2 == 0, ]
PESODD <- PESRe[PESRe$Trial %% 2 == 1, ]

PESEVEN <- PESEVEN %>%
  group_by(ID,PES) %>%
  mutate(RT.PESEVEN = mean(RT)) %>%
  select("ID", "PES", "RT.PESEVEN") %>%
  ungroup() 

PESEVEN <-  PESEVEN %>%
  group_by(ID) %>%
  distinct() %>%
  ungroup()

PESEVEN$ID <- as.character(PESEVEN$ID)
PESEVEN <-  PESEVEN %>%
  pivot_wider(
    names_from = PES,
    values_from = RT.PESEVEN
  ) 

PESEVEN <-  PESEVEN %>%
  mutate(RT.diff.PES.EVEN = PE - NPE)

PESODD <- PESODD %>%
  group_by(ID,PES) %>%
  mutate(RT.PESODD = mean(RT)) %>%
  select("ID", "PES", "RT.PESODD") %>%
  ungroup() 

PESODD <-  PESODD %>%
  group_by(ID) %>%
  distinct() %>%
  ungroup()

PESODD$ID <- as.character(PESODD$ID)
PESODD <-  PESODD %>%
  pivot_wider(
    names_from = PES,
    values_from = RT.PESODD
  ) 

PESODD <-  PESODD %>%
  mutate(RT.diff.PES.ODD = PE - NPE)

PES <- merge(PESODD,PESEVEN,by="ID")

cor.test(PES$RT.diff.PES.ODD,PES$RT.diff.PES.EVEN)


################################################################################################################

## Zero order correlation

parameters_M <- parameters_M %>% select("SES","Mother_Age":"RT.diff.PES")

Hmisc::rcorr(as.matrix(as.data.frame(lapply(parameters_M, as.numeric))))

p.adjust(c(0.926,0,0.253,0.4115,0.0619,0.1258,0.1939,
           0.1931,0.4988,0.9436,0.8589,0.4355,0.6237,
           0.1913,0,0,0.828,0.8774,
           0.5481,0.3049,0.1566,0.9626,
           0.035,0.13,0.0004,
           0.6946,0,
           0.9314), method = "BH")

# p value after BH adjustment

# 0.9626000 0.0000000 0.5060000 0.7172941 0.2476000 0.4044444 0.4176308
# 0.4176308 0.7759111 0.9626000 0.9626000 0.7172941 0.8731800
# 0.4176308 0.0000000 0.0000000 0.9626000 0.9626000
# 0.8077263 0.5691467 0.4176308 0.9626000
# 0.1633333 0.4044444 0.0022400
# 0.9261333 0.0000000
# 0.9626000

Hmisc::rcorr(as.matrix(as.data.frame(lapply(parameters_M, as.numeric))),type = "spearman")

## specific cor test

cor.test(parameters_M$RT.diff.DD,parameters_M$RT.diff.PES)
cor.test(parameters_M$PI_prop,parameters_M$RT.diff.DD)
cor.test(parameters_M$PE_prop,parameters_M$RT.diff.PES)

## Descriptive statistics

summary(parameters_M)

parameters_M %>% summarise(sdSES = sd(SES),
                           sdAge = sd(Mother_Age),
                           sdRaven = sd(Raven),
                           sdPI_prop = sd(PI_prop),
                           sdPE_prop = sd(PE_prop),
                           sdRT = sd(mean_RT),
                           sdDD = sd(RT.diff.DD),
                           sdPES = sd(RT.diff.PES,na.rm = T))


###################################################################################################################################

## HLM mother

## Multiple regression analyses

# DD

lm1 <- lm(RT.diff.DD ~
            1
          ,data = parameters_M)
summary(lm1)


lm2 <- lm(RT.diff.DD ~
            scale(PI_prop) +
            scale(Raven) +
            scale(Mother_Age) 
          ,data = parameters_M)
summary(lm2)

anova(lm1,lm2)


# PES

lm1 <- lm(RT.diff.PES ~
            1
          ,data = parameters_M)
summary(lm1)


lm2 <- lm(RT.diff.PES ~
            scale(PE_prop) +
            scale(Raven) +
            scale(Mother_Age) 
          ,data = parameters_M)
summary(lm2)

anova(lm1,lm2)

