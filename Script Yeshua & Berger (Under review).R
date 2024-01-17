library(tidyverse)
library(car)
library(MKinfer)
library(lmerTest)
library(performance)
library(parameters)
library(ggeffects)
library(patchwork)
library(emmeans)

HLM_PMC <- read.csv("HLM_PMC.csv")
parameters_C <- read.csv("parameters_ps.csv") 

HLM_PMC$delayed_disinhibition <- factor(HLM_PMC$delayed_disinhibition, levels = c("NDD","DD"))
HLM_PMC <- HLM_PMC %>% mutate(delayed_disinhibition = ifelse(delayed_disinhibition == "NDD","NPI","PI"))

HLM_PMC$Sex_C <- factor(HLM_PMC$Sex_C, levels = c("Male","Female"))

SES <- parameters_C %>% select("Education_Years_1","Income","Home_2","Home_3")
SES <- as.data.frame(apply(SES, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x)))


factor_model <- factanal(
  SES,
  factors = 1,       
  method = "ml",  
  rotation = "varimax",
  scores = "regression")
print(factor_model) 
SES$SES <- factor_model$scores
parameters_C$SES <- factor_model$scores 

parameters_C <- parameters_C %>% select("SES","Sex_C","Child_Age","RavenScore":"RT.diff.PES")

TrialsProportion <- HLM_PMC %>%
  mutate(DD_dummy = delayed_disinhibition == "PI",
         PES_dummy = PES == "PE") %>%
  group_by(ID) %>%
  mutate(PI_prop = mean(DD_dummy),
         PE_prop = mean(PES_dummy)) %>%
  select("ID","PI_prop","PE_prop","Sex_C","Child_Age","RavenScore","SES") %>%
  distinct() %>%
  ungroup() 


parameters_C <- merge(TrialsProportion,parameters_C) %>% select("SES","Sex_C","Child_Age","RavenScore",
                                                                "PI_prop","PE_prop","mean_RT":"RT.diff.PES")

### Descriptives

parameters_C %>% summarise(n = n(),
                          Age = mean(Child_Age),
                          SD = sd(Child_Age),
                          min = min(Child_Age),
                          max = max(Child_Age))

parameters_C %>% group_by(Sex_C) %>% summarise(n = n(),
                                               Age = mean(Child_Age),
                                               SD = sd(Child_Age),
                                               min = min(Child_Age),
                                               max = max(Child_Age))

HLM_PMC %>% group_by(ID) %>%
  summarise(n = n()) %>% 
  ungroup() %>%
  summarise(mean = mean(n),
            sd = sd(n),
            min = min(n),
            max = max(n))

### Reliability test

#######################################################

# Delayed Disinhibition

DDRe <- HLM_PMC %>%
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

PESRe <- HLM_PMC %>%
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

parameters_C$Sex_C <- as.factor(parameters_C$Sex_C)
parameters_C$Sex_C <- factor(parameters_C$Sex_C, levels = c("Male","Female"))

Hmisc::rcorr(as.matrix(as.data.frame(lapply(parameters_C, as.numeric))))

p.adjust(c(0.1856,0.2825,0.7354,0.2036,0.9869,0.5299,0.3687,0.1615,
           0.743,0.0233,0.6245,0.8828,0.066,0.0027,0.6085,
           0,0.8498,0.0387,0.0043,0.6334,0.0381,
           0.9989,0.0057,0,0.0817,0.019,
           0.7051,0.6759,0,0.7521,
           0.1327,0.5657,0.0003,
           0.0031,0.141,
           0.011), method = "BH")

# p value after BH adjustment

# 0.35166316 0.48428571 0.84611250 0.36648000 0.99890000 0.82940870 0.60332727 0.32300000 
# 0.84611250 0.07625455 0.84453333 0.93472941 0.16971429 0.01860000 0.84453333
# 0.00000000 0.92705455 0.10716923 0.02211429 0.84453333 0.10716923
# 0.99890000 0.02565000 0.00000000 0.19608000 0.06840000
# 0.84611250 0.84611250 0.00000000 0.84611250
# 0.29857500 0.84453333 0.00270000
# 0.01860000 0.29858824
# 0.04400000


Hmisc::rcorr(as.matrix(as.data.frame(lapply(parameters_C, as.numeric))),type = "spearman")


## specific cor test

cor.test(parameters_C$RT.diff.DD,parameters_C$RT.diff.PES)
cor.test(parameters_C$PI_prop,parameters_C$RT.diff.DD)
cor.test(parameters_C$PE_prop,parameters_C$RT.diff.PES)

## Descriptive statistics

summary(parameters_C)

parameters_C %>% summarise(sdSES = sd(SES),
                           sdAge = sd(Child_Age),
                           sdRaven = sd(RavenScore),
                           sdPIprop = sd(PI_prop),
                           sdPEprop = sd(PE_prop),
                           sdRT = sd(mean_RT),
                           sdDD = sd(RT.diff.DD),
                           sdPES = sd(RT.diff.PES,na.rm = T))


###################################################################################################################################
## Multiple regression analyses

# DD

lm1 <- lm(RT.diff.DD ~
           1
           ,data = parameters_C)
summary(lm1)


lm2 <- lm(RT.diff.DD ~
           scale(PI_prop) +
           scale(RavenScore) +
           Sex_C +
           scale(Child_Age) 
         ,data = parameters_C)
summary(lm2)

lm3 <- lm(RT.diff.DD ~
            scale(PI_prop) +
            scale(RavenScore) *
           Sex_C +
          scale(Child_Age) *
           Sex_C
         ,data = parameters_C)

summary(lm3)

anova(lm1,lm2,lm3)

ggemmeans(lm3, c("Child_Age","Sex_C"))
G1 <- ggemmeans(lm3, c("Child_Age","Sex_C"))
plot(G1) + labs(
  title = "",
  x = "Age",
  y = "DD effect",
  colour = "Sex") +
  scale_color_manual(values = c("gray", "black")) +  # Set custom colors for Sex variable
  scale_fill_manual(values = c("lightgray", "darkgray")) + # Set custom colors for CI
  theme_classic()

emmeans(lm3, c("Child_Age","Sex_C"),
                  at = list(
                    Child_Age = seq(min(parameters_C$Child_Age),max(parameters_C$Child_Age), len = 30),
                    Sex_C = c("Male","Female")))

joint_tests(lm3,by = "Sex_C")


# PES


lm1 <- lm(RT.diff.PES ~
            1
          ,data = parameters_C)
summary(lm1)


lm2 <- lm(RT.diff.PES ~
            scale(PE_prop) +
            scale(RavenScore) +
            Sex_C +
            scale(Child_Age) 
          ,data = parameters_C)
summary(lm2)

lm3 <- lm(RT.diff.PES ~
            scale(PE_prop) +
            scale(RavenScore) *
            Sex_C +
            scale(Child_Age) *
            Sex_C
          ,data = parameters_C)

summary(lm3)

anova(lm1,lm2,lm3)

ggemmeans(lm2, c("Child_Age"))
G1 <- ggemmeans(lm2, c("Child_Age"))
plot(G1) + labs(
  title = "",
  x = "Age",
  y = "RT") +
  scale_color_grey() +
  theme_classic() 
