library(tidyverse)
library(prepdat)
library(psych)
library(lme4)
library(misty)
library(lmerTest)
library(sjPlot)
library(ggplot2)
library(reghelper)
library(parameters)
library(bayestestR)
library(car)

parameters_CFILTER <- read.csv("parameters_CFILTER.csv")

Hmisc::rcorr(as.matrix(as.data.frame(lapply(parameters_CFILTER, as.numeric))))

ggplot(parameters_CFILTER, aes(RT.diff.DD,RT.diff.PES)) +
  geom_point() +
  stat_smooth(method = "lm", color = "black") +
  labs(x="DD standerdized effect",y="PES standerdized effect") +
  scale_color_grey() +
  theme_classic()

parameters_CFILTER$Sex_C <- as.factor(parameters_CFILTER$Sex_C)

ggplot(parameters_CFILTER, aes(Sex_C,RT.diff.DD)) +
  geom_boxplot() +
  labs(x="Sex",y="DD standerdized effect") +
  scale_x_discrete(labels=c("0" = "Male", "1" = "Female")) +
  scale_color_grey() +
  theme_classic()

ggplot(parameters_CFILTER, aes(Child_Age,RT.diff.PES)) +
  geom_point() +
  stat_smooth(method = "lm", color = "black") +
  labs(x="Children's age",y="PES standerdized effect") +
  scale_color_grey() +
  theme_classic()

ggplot(parameters_CFILTER, aes(RavenScore,RT.diff.PES)) +
  geom_point() +
  stat_smooth(method = "lm", color = "black") +
  labs(x="Raven score",y="PES standerdized effect") +
  scale_color_grey() +
  theme_classic()

ggplot(parameters_CFILTER, aes(RavenScore,Child_Age)) +
  geom_point() +
  stat_smooth(method = "lm", color = "black") +
  labs(x="Raven score",y="Children's age") +
  scale_color_grey() +
  theme_classic()


summary(parameters_CFILTER)
sd(parameters_CFILTER$RT.diff.DD) 
sd(parameters_CFILTER$RT.diff.PES) 
sd(parameters_CFILTER$Child_Age) 
sd(parameters_CFILTER$RavenScore) 

cor.test(parameters_CFILTER$RT.diff.DD,
         parameters_CFILTER$RT.diff.PES)

cor.test(parameters_CFILTER$RT.diff.PES,
         parameters_CFILTER$Child_Age)

cor.test(parameters_CFILTER$Child_Age,
         parameters_CFILTER$RavenScore)

leveneTest(RT.diff.DD ~ Sex_C, parameters_CFILTER, center = "mean")

## 0 = Male 1 = Female

t.test(parameters_CFILTER$RT.diff.DD[parameters_CFILTER$Sex_C == "0"],
       parameters_CFILTER$RT.diff.DD[parameters_CFILTER$Sex_C == "1"],
       var.equal = T)


leveneTest(RT.diff.PES ~ Sex_C, parameters_CFILTER, center = "mean")

t.test(parameters_CFILTER$RT.diff.PES[parameters_CFILTER$Sex_C == "0"],
       parameters_CFILTER$RT.diff.PES[parameters_CFILTER$Sex_C == "1"],
       var.equal = T)

PivoteData <- parameters_CFILTER %>% pivot_longer(names_to = "Conds",
                                                  values_to = "RT",
                                                  cols = c(DD:NPE)) 

leveneTest(RT ~ Conds, PivoteData, center = "mean")

t.test(PivoteData$RT[PivoteData$Conds == "DD"],
       PivoteData$RT[PivoteData$Conds == "NDD"],
       var.equal = T)

t.test(PivoteData$RT[PivoteData$Conds == "PE"],
       PivoteData$RT[PivoteData$Conds == "NPE"],
       var.equal = T)


### HLM

HLM2 <- read.csv("HLM2.csv")

HLM2 %>%
  group_by(ID) %>% 
  summarise(n = n()) %>%
  ungroup() %>% print(n = 200) ## Sample size = 67


HLM2 %>% group_by(ID,PES) %>% 
  summarise(n = n()) %>%
  ungroup() %>% print(n = 200) 

HLM2 %>% group_by(ID,delayed_disinhibition) %>% 
  summarise(n = n()) %>%
  ungroup() %>% print(n = 200) 

HLM2$delayed_disinhibition <- factor(HLM2$delayed_disinhibition, 
                                     levels = c("NDD", "DD"))

hist(HLM2$Child_Age)
min(HLM2$Child_Age, na.rm = T)
max(HLM2$Child_Age, na.rm = T)
mean(HLM2$Child_Age, na.rm = T)
sd(HLM2$Child_Age, na.rm = T)

model0 <- lmer(RT ~ PES +
                 delayed_disinhibition +
                 RavenScore +
                 Sex_C +
                 Child_Age +
                 (PES + delayed_disinhibition | ID),
               data = HLM2)
summary(model0)
model_parameters(model0)

library(lmerTest)
library(performance)
r2_nakagawa(model0)

model1 <- lmer(RT ~ (PES +
                       delayed_disinhibition +
                       Sex_C +
                       Child_Age)^2 +
                 RavenScore +
                 (PES + delayed_disinhibition | ID),
               data = HLM2)
summary(model1)
model_parameters(model1)
r2_nakagawa(model1, tolerance = 0)

model1.2 <- lmer(RT ~ PES *
                   delayed_disinhibition *
                   Sex_C *
                   Child_Age +
                   RavenScore +
                   (PES + delayed_disinhibition | ID),
                 data = HLM2)
summary(model1.2)
model_parameters(model1.2)
r2_nakagawa(model1.2, tolerance = 0)

anova(model0,model1,model1.2)

library(emmeans)


HLM2$Sex_C <- factor(HLM2$Sex_C, 
                     levels = c("Female", "Male"))

em_DDPE <- emmeans(model1.2, c("Sex_C", "Child_Age", "PES", "delayed_disinhibition"),
                   at = list(
                     Child_Age = seq(min(HLM2$Child_Age),max(HLM2$Child_Age), len = 20),
                     PES = "PE",
                     delayed_disinhibition = "DD"
                   ))
c_DDPE <- contrast(em_DDPE, method = "pairwise", by = "Child_Age") %>% 
  update(by = NULL, adjust = "none")

emmip(c_DDPE, ~ Child_Age, CIs = TRUE,
      xlab = "Age (Z scores)",
      ylab = "Linear prediction of high cognitive load effect") + 
  geom_hline(yintercept = 0) +
  scale_color_grey() +
  theme_classic() 

library(ggplot2)


HLM2 <- HLM2 %>% mutate(Condition = ifelse(PES == "PE" &
                                             delayed_disinhibition == "DD", "PE + PI",
                                           ifelse(PES == "PE" &
                                                    delayed_disinhibition == "NDD", "PE + NPI",
                                                  ifelse(PES == "NPE" &
                                                           delayed_disinhibition == "DD", "NPE + PI",
                                                         "NPE + NPI"))))

HLM2 <- HLM2 %>% mutate(delayed_disinhibition = ifelse(delayed_disinhibition == "DD","PI","NPI"))
HLM2 <- HLM2 %>% rename(`Delayed Disinhibition` = "delayed_disinhibition") 


ggplot(HLM2, 
       aes(x=Child_Age,y=RT, color = `Delayed Disinhibition`,
           group = `Delayed Disinhibition`)) +
  stat_smooth(method = "lm") +
  facet_grid(PES~Sex_C, scales = 'free') +
  xlab("Age (Z scores)") +
  ylab("Response Time (ms)") +
  scale_color_grey() +
  theme_classic()

ggplot(HLM2, 
       aes(x=Child_Age,y=RT, color = PES,
           group = PES)) +
  stat_smooth(method = "lm") +
  facet_grid(`Delayed Disinhibition`~Sex_C, scales = 'free') +
  xlab("Age (Z scores)") +
  ylab("Response Time (ms)") +
  scale_color_grey() +
  theme_classic()

HLM2 <- HLM2 %>% mutate(delayed_disinhibition = ifelse(delayed_disinhibition == "PI","DD","NDD"))
HLM2 <- HLM2 %>% rename("delayed_disinhibition" = `Delayed Disinhibition`) 


library(emmeans)
emt <- emmeans::emmeans(model1.2, ~ PES + delayed_disinhibition + Child_Age + Sex_C,
                        at = list(Child_Age = seq(-2.47, 2.27, len = 20)),
                        lmer.df = "S")
emt[1:5]

pe_contrasts <- contrast(emt, method = "revpairwise", by = c("Child_Age", "Sex_C",
                                                             "delayed_disinhibition"),
                         adjust = "none") %>% 
  update(by = "Sex_C")

emmip(pe_contrasts, delayed_disinhibition ~ Child_Age | Sex_C, CIs = TRUE,
      xlab = "Age (Z scores)",
      ylab = "Linear prediction of PES effect") + 
  guides(color = guide_legend(title = "DD")) +
  geom_hline(yintercept = 0) +
  scale_color_grey() +
  theme_classic() 


pe_contrasts <- contrast(emt, method = "revpairwise", by = c("Child_Age", "Sex_C",
                                                             "PES"),
                         adjust = "none") %>% 
  update(by = "Sex_C")

emmip(pe_contrasts, PES ~ Child_Age | Sex_C, CIs = TRUE,
      xlab = "Age (Z scores)",
      ylab = "Linear prediction of DD effect") + 
  geom_hline(yintercept = 0) +
  scale_color_grey() +
  theme_classic()

emt <- emmeans::emmeans(model1.2, ~ delayed_disinhibition + Child_Age + Sex_C,
                        at = list(Child_Age = seq(-2.3, 2.16, len = 20)),
                        lmer.df = "S")
emt[1:5]

pe_contrasts <- contrast(emt, method = "revpairwise", by = c("Child_Age", "Sex_C"),
                         adjust = "none") %>% 
  update(by = "Sex_C") 

emmip(pe_contrasts, Sex_C  ~ Child_Age, CIs = TRUE,
      xlab = "Age (Z scores)",
      ylab = "Linear prediction of DD effect") +
  guides(color = guide_legend(title = "Sex")) +
  geom_hline(yintercept = 0) +
  scale_color_grey() +
  theme_classic()