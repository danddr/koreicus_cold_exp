setwd("/home/dared/GitHub/koreicus_cold_exp/")
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(ggeffects)
library(broom)
library(broom.mixed)

 # overdispersion function
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type="pearson")
  Pearson.chisq <- sum(rp^2)
  ratio <- Pearson.chisq / rdf
  p <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  return(c(ChiSq=Pearson.chisq, Ratio=ratio, df=rdf, p=p))
}

#set out directory
outdir<-"outputs/"
#set color palette
myCols <- c("Diapausing"= "#003366", "Non-Diapausing"= "#FF8E00")

#---- 1. Load Control and Treatment observations ----
df_treat <- read_excel(path="data/koreicus_fluct_cold_eggs_data.xlsx",sheet=1) %>% 
  select(Type, Temperature, Period_length, Cup_N, N_EGGS, TOTAL) %>% 
  rename(TypeEggs = Type, replicate = Cup_N, start_eggs = N_EGGS, hatched = TOTAL) %>% 
  mutate(
    ExpSetting = "Treatment",
    TypeEggs = ifelse(TypeEggs == "D", "Diapausing", "Non-Diapausing"),
    replicate = as.factor(replicate),
    # treat.hatchRate = hatched / start_eggs,
    join.id = paste(TypeEggs, abs(Temperature), Period_length, replicate, sep="_")  
  )

df_contr <- read_excel(path="data/koreicus_fluct_cold_eggs_data.xlsx",sheet=2) %>% 
  select(Type, Confront_Temperature, Period_length, Cup_N, N_EGGS, TOTAL) %>% 
  rename(TypeEggs = Type, Temperature = Confront_Temperature, replicate = Cup_N, start_eggs = N_EGGS, hatched = TOTAL) %>% 
  mutate(
    ExpSetting = "Control",
    TypeEggs = ifelse(TypeEggs == "D", "Diapausing", "Non-Diapausing"),
    replicate = as.factor(replicate),
    # control.hatchRate = hatched / start_eggs,
    join.id = paste(TypeEggs, abs(Temperature), Period_length, replicate, sep="_")
  ) 
# rename(control_hatched = hatched, control_start_eggs = start_eggs)

# put together
df_stat <- bind_rows(df_treat, df_contr) %>%
  mutate(
    Period_length = paste(Period_length, "days"),
    Period_length = factor(Period_length, levels = c("2 days", "5 days", "10 days", "20 days", "30 days")),
    TypeEggs = factor(TypeEggs, levels = c("Non-Diapausing", "Diapausing"))
  )
outname <- paste0(outdir, "Aedes_koreicus_eggsHatching_coldExposure.csv")
write.csv(df_stat, outname)

##---- 1.1 summary table ----
# aggregated summary 
df_stat %>%
  mutate(prop_hatched = hatched / start_eggs) %>%
  group_by(TypeEggs, ExpSetting) %>%
  summarise(
    mean_prop = mean(prop_hatched, na.rm = TRUE),
    se_prop = sd(prop_hatched, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    mean_se = sprintf("%.3f ± %.3f", mean_prop, se_prop)
  ) 

# brokem dowm summary
df_summary <- df_stat %>%
  mutate(prop_hatched = hatched / start_eggs) %>%
  group_by(TypeEggs, Temperature, Period_length, ExpSetting) %>%
  summarise(
    mean_prop = mean(prop_hatched, na.rm = TRUE),
    se_prop = sd(prop_hatched, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    mean_se = sprintf("%.3f ± %.3f", mean_prop, se_prop)
  ) 

df_summary_wide <- df_summary %>% 
  # ensure all combinations are represented, even if missing
  complete(TypeEggs, Temperature, ExpSetting, Period_length) %>%
  select(TypeEggs, Temperature, ExpSetting, Period_length, mean_se) %>%
  pivot_wider(
    names_from = Period_length,
    values_from = mean_se
  )

df_summary_wide

##---- 1.2 summary plot ----
df_summary %>% 
  ggplot()+
  geom_point(aes(x = Temperature, y = mean_prop, col = ExpSetting), size = 2.8 ) +
  geom_errorbar(aes(x = Temperature, ymin = mean_prop - se_prop, ymax = mean_prop + se_prop ,col = ExpSetting), 
                width = 0.2, size = 0.8)+
    ylim(0,1)+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  scale_color_manual(  values = c("#ff99be","#4e8bc4" ),
                       labels = c("Control (constant +5°C)", "Treatment (fluctuating < 0°C)"))+
  labs(x = "Temperature (°C)", y="Hatched eggs (%)", col="")+
  facet_grid(TypeEggs~Period_length)+
  theme_classic()+
  coord_fixed(ratio=20)+
  theme(legend.position = "bottom",
        text = element_text(size=14),
        legend.key.size = unit(1.5, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )

##---- 1.3 Treatment Vs Control model (Eq.1) ----
###---- 1.3.1 Modelling ----
bin0<-glmmTMB(cbind(hatched, start_eggs-hatched) ~ ExpSetting + TypeEggs,
              data=df_stat, family=binomial)
#check for overdispersion in the binomial GLM
overdisp_fun(bin0) # some overdisperison

# Simpler model (check if interaction is needed)
betabin0 <- glmmTMB(cbind(hatched, start_eggs-hatched) ~ ExpSetting + TypeEggs,
                    data=df_stat, family=betabinomial)
# summary(betabin0)
overdisp_fun(betabin0)
AIC(bin0, betabin0)

betabin1<-glmmTMB(cbind(hatched, start_eggs-hatched)~ExpSetting*TypeEggs + (1|replicate), 
                  data=df_stat, family=betabinomial)
summary(betabin1)
AIC(betabin0, betabin1)
anova(betabin0, betabin1) # ok we keep the interaction

betabin2<-glmmTMB(cbind(hatched, start_eggs-hatched)~ExpSetting*TypeEggs, 
                  data=df_stat, family=betabinomial)
AIC(betabin1, betabin2)
anova(betabin1, betabin2)
summary(betabin2)

betabin3<-glmmTMB(cbind(hatched, start_eggs-hatched)~ExpSetting*TypeEggs*Period_length, 
                  data=df_stat, family=betabinomial)
AIC(betabin2, betabin3)
anova(betabin2, betabin3)
summary(betabin3) #too complex, we stay with betabin2

# model diagnostic
res <- simulateResiduals(betabin2)
plot(res)
DHARMa::testDispersion(res)
DHARMa::testZeroInflation(res)
testOutliers(res, type="bootstrap")
testUniformity(res) 

###---- 1.3.2 Model coefficients table ----
# Extract tidy coefficients from glmmTMB model
summary(betabin2)
tab_betabin2 <- tidy(betabin2, effects = "fixed") %>%
  mutate(
    # Compute 95% confidence intervals (on β scale)
    conf.low = estimate - 1.96 * std.error,
    conf.high = estimate + 1.96 * std.error
  ) %>%
  mutate(
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    conf.low = round(conf.low, 3),
    conf.high = round(conf.high, 3),
    p.value = ifelse(p.value < 0.001, "<0.001", round(p.value, 3))
  ) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high)

# Print the table
tab_betabin2

###---- 1.3.3 Model predictions ----
# Get model predictions
preds <- ggpredict(betabin2, terms = c("ExpSetting", "TypeEggs"))
preds 

# Plot
ggplot(preds, aes(x = x, y = predicted, color = group, group = group)) +
  geom_point(position = position_dodge(0.2), size = 3) +
  # geom_line(position = position_dodge(0.2), linewidth = 1) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                position = position_dodge(0.2), width = 0.1) +
  scale_color_manual(values = c("#ff99be","#4e8bc4")) +
  labs(x = "Experimental setting",
       y = "Predicted hatching probability",
       color = "Egg type",
       # title = "Predicted hatching probability by egg type and treatment"
       ) +
  theme_classic(base_size = 16)+
  theme(legend.position = "bottom")

#---- 2. Temperature, Eggs and cold Exposure model (Eq. 2) ----
# Load data
df_treat <- read_excel(path="data/koreicus_fluct_cold_eggs_data.xlsx",sheet=1) %>% 
  select(Type, Temperature, Period_length, Cup_N, N_EGGS, TOTAL) %>% 
  rename(TypeEggs = Type, replicate = Cup_N, start_eggs = N_EGGS, hatched = TOTAL) %>% 
  mutate(
    ExpSetting = "Treatment",
    TypeEggs = ifelse(TypeEggs == "D", "Diapausing", "Non-Diapausing"),
    replicate = as.factor(replicate),
    # treat.hatchRate = hatched / start_eggs,
    join.id = paste(TypeEggs, abs(Temperature), Period_length, replicate, sep="_")  
  )


p <- df_treat %>% 
  mutate(Period_length = paste(Period_length, "days"), 
         Period_length = factor(Period_length, c("2 days", "5 days", "10 days", "20 days", "30 days"))) %>%
  ggplot()+
  stat_summary(aes(x = Temperature, y = hatched / start_eggs, col = TypeEggs),
               fun = mean, geom = "point", size = 1.8 ) +
  stat_summary(aes(x = Temperature, y = hatched / start_eggs, col = TypeEggs),
    fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8) +
  ylim(0,1)+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  scale_color_manual(values = myCols)+
  labs(x = "Temperature (°C)", y="Hatched eggs (%)", col="Eggs")+
  facet_wrap(~Period_length, ncol=5)+
  theme_classic()+
  coord_fixed(ratio=20)+
  theme(legend.position = "bottom",
        text = element_text(size=16),
        # legend.key.size = unit(1.5, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )
p
outname <- paste0(outdir, "/hatch_observed",Sys.Date(), ".png")
ggsave(p, filename = outname, width = 10, height = 10, device='png', dpi=320)

##---- 2.1  Additive Modelling ----
# additive model with interactions
bin1 <-glmmTMB(cbind(hatched, start_eggs-hatched)~Temperature*TypeEggs + TypeEggs*Period_length  + Temperature*Period_length, 
               data=df_treat, family=binomial)
#check for overdispersion in the binomial GLM
overdisp_fun(bin1) # some overdisperison
betabin<-glmmTMB(cbind(hatched, start_eggs-hatched)~Temperature*TypeEggs + TypeEggs*Period_length + Temperature*Period_length, data=df_treat, family=betabinomial)
betabin2<-glmmTMB(cbind(hatched, start_eggs-hatched)~Temperature*TypeEggs + TypeEggs*Period_length  + Temperature*Period_length + (1|replicate), data=df_treat, family=betabinomial)
overdisp_fun(betabin)
overdisp_fun(betabin2)
AIC(bin1, betabin, betabin2)

# Likelihood Ratio Test
anova(bin1, betabin) 
anova( betabin, betabin2)
# betabinomial without random effect seems the best model
summary(betabin)

#try quadratic effect on temperature
betabin3<-glmmTMB(cbind(hatched, start_eggs-hatched)~poly(Temperature,2)*TypeEggs + TypeEggs*Period_length + poly(Temperature,2)*Period_length, data=df_treat, family=betabinomial)
AIC(betabin, betabin3) 
anova( betabin, betabin3)

#betabin check
sim_fmp <- DHARMa::simulateResiduals(betabin3) 
DHARMa::testZeroInflation(sim_fmp)
testOutliers(sim_fmp, type="bootstrap")
testDispersion(sim_fmp)     # Dispersion test
testUniformity(sim_fmp) 
plot(sim_fmp)
plot(df_treat$hatched/df_treat$start_eggs, sim_fmp$scaledResiduals)
abline(0,1)

##---- 2.2  Three-ways interaction Additive Modelling ----
three.inter_bin1<-glmmTMB(cbind(hatched, start_eggs-hatched)~Temperature*TypeEggs*Period_length, data=df_treat, family=binomial)
#check for overdispersion in the binomial GLM
overdisp_fun(three.inter_bin1) # some overdisperison
three.inter_betabin<-glmmTMB(cbind(hatched, start_eggs-hatched)~Temperature*TypeEggs*Period_length, data=df_treat, family=betabinomial)
three.inter_betabin2<-glmmTMB(cbind(hatched, start_eggs-hatched)~Temperature*TypeEggs*Period_length + (1|replicate), data=df_treat, family=betabinomial)
overdisp_fun(three.inter_betabin)
overdisp_fun(three.inter_betabin2)
AIC(three.inter_bin1, three.inter_betabin, three.inter_betabin2)

# Likelihood Ratio Test
anova(three.inter_bin1, three.inter_betabin)  
anova( three.inter_betabin, three.inter_betabin2)
summary(betabin)
#again, simple betabinomial model is better

# try quadratic effect on temperature
three.inter_betabin3<-glmmTMB(cbind(hatched, start_eggs-hatched)~poly(Temperature,2)*TypeEggs*Period_length, data=df_treat, family=betabinomial)
AIC(three.inter_bin1, three.inter_betabin, three.inter_betabin2, three.inter_betabin3)
anova(three.inter_betabin, three.inter_betabin3) # quadratic betabin better 

summary(three.inter_betabin3)

#betabin check
sim_fmp <- DHARMa::simulateResiduals(three.inter_betabin3) 
DHARMa::testZeroInflation(sim_fmp)
testOutliers(sim_fmp, type="bootstrap")
testDispersion(sim_fmp)     # Dispersion test
testUniformity(sim_fmp) 
plot(sim_fmp)
plot(df_treat$hatched/df_treat$start_eggs, sim_fmp$scaledResiduals)

##---- 2.3  Best models comparisons ----
AIC(betabin3, three.inter_betabin3)
anova(betabin3, three.inter_betabin3)

sim_fmp <- DHARMa::simulateResiduals(betabin3) 
DHARMa::testZeroInflation(sim_fmp)
testOutliers(sim_fmp, type="bootstrap")
testDispersion(sim_fmp)     # Dispersion test
testUniformity(sim_fmp) 
plot(sim_fmp)
plot(df_treat$hatched/df_treat$start_eggs, sim_fmp$scaledResiduals)


###---- 1.3.2 Model coefficients table ----
# Extract tidy coefficients from glmmTMB model
summary(betabin3)
tab_betabin3 <- tidy(betabin3, effects = "fixed") %>%
  mutate(
    # Compute 95% confidence intervals (on β scale)
    conf.low = estimate - 1.96 * std.error,
    conf.high = estimate + 1.96 * std.error
  ) %>%
  mutate(
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    conf.low = round(conf.low, 3),
    conf.high = round(conf.high, 3),
    p.value = ifelse(p.value < 0.001, "<0.001", round(p.value, 3))
  ) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high)

# Print the table
tab_betabin3

##---- 2.4 get marginal means ----
# Choose specific values to evaluate at
ref_grid <- ref_grid(betabin3, 
                     at = list(Temperature = c(-20, -15, -10, -5, 0), 
                               Period_length = c(2, 5, 10, 20, 30)))

# Get estimated marginal means by TypeEggs
emm <- emmeans(ref_grid, ~ TypeEggs | Temperature + Period_length, type = "response")
summary(emm)

# get two heatmaps
emm_df <- as.data.frame(emm)
p <- emm_df %>% 
  ggplot(aes(x = Temperature, y = as.factor(Period_length), fill = prob)) +
  geom_tile() +
  scale_fill_viridis_c(option = "mako", name = "Proportion of hatched eggs", limits=c(0,1)) +
  labs(
       x = "Temperature (°C)", y = "Cold exposure (days)") +
  facet_wrap(~TypeEggs)+
  theme_classic()+
  theme(legend.position = "bottom",
        text = element_text(size=16),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )
outname <- paste0(outdir, "/heatmap",Sys.Date(), ".png")
ggsave(p, filename = outname, width = 15, height = 10, device='png', dpi=320)


##---- 2.5 Model predictions ----
#predict to intepred the model output
pred <- expand.grid(Temperature=seq(min(df_treat$Temperature),max(df_treat$Temperature)),
                    Period_length=unique(df_treat$Period_length),
                    TypeEggs=c("Diapausing","Non-Diapausing"))
pred$hatched <- predict(betabin3, pred, type = "response")
pred<- pred %>% mutate(se.fit = predict(betabin3, pred, se.fit = TRUE)$se.fit)

pred %>%
    mutate(Period_length = paste(Period_length, "days"), 
       Period_length = factor(Period_length, c("2 days", "5 days", "10 days", "20 days", "30 days"))) %>% 
ggplot()+
  geom_line( aes(x=Temperature, y=hatched, col=TypeEggs))+
  geom_ribbon(aes(x=Temperature, ymin = hatched - se.fit, ymax = hatched + se.fit, fill = TypeEggs), alpha = 0.2, color = NA)+
  geom_point(data= df_treat%>%
               mutate(Period_length = paste(Period_length, "days"), 
                      Period_length = factor(Period_length, c("2 days", "5 days", "10 days", "20 days", "30 days"))),
             aes(x=Temperature, y=hatched/start_eggs, col=TypeEggs))+
  ylim(-0.25,1.1)+
  geom_hline(yintercept = 0.5, linetype="dashed")+
  scale_color_manual(values = myCols)+
  scale_fill_manual(values = myCols)+
  labs(x = "Temperature (°C)", y="Hatched eggs (%)", col="Eggs", fill="Eggs")+
  facet_wrap(~Period_length, ncol=5)+
  theme_classic()+
  coord_fixed(ratio=20)+
  theme(legend.position = "bottom",
        text = element_text(size=14),
        legend.key.size = unit(1.5, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )

#---- 3. Ae. koreicus Vs. Ae. albopictus (temperate strain from Tippelt et al 2020) comparison ----
albo <- read_excel(path="data/Tippelt2020_Fig2.xlsx",sheet=1) %>% 
  select(`life stage`, Temperature, `cold exposure`, replicate, Value, `Initial number of individuals`, species) %>% 
  rename(TypeEggs = `life stage`, Period_length = `cold exposure`, 
         start_eggs = `Initial number of individuals` , prop_hatched = Value) %>% 
  mutate(
    ExpSetting = "Treatment",
    # Period_length = as.factor(paste(Period_length, "days")),
    TypeEggs = ifelse(TypeEggs == "Non-diapausing egg", "Non-Diapausing", "Diapausing"),
    replicate = as.factor(replicate), 
    hatched = round(start_eggs*prop_hatched),
  )

# koreicus
kor <- df_stat %>%
  filter(ExpSetting == "Treatment") %>% 
  mutate(prop_hatched = hatched / start_eggs,
         species = "Aedes koreicus", 
         Period_length =as.numeric(gsub(" days", "", Period_length))) %>% 
  select(all_of(names(albo)))

# combine
cfr.df <- bind_rows(kor, albo)
cfr.df
##---- 3.1 summary stat ----
cfr.df %>% 
  group_by(species, TypeEggs) %>%
  summarise(
    mean_prop = mean(prop_hatched, na.rm = TRUE),
    se_prop = sd(prop_hatched, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    mean_se = sprintf("%.3f ± %.3f", mean_prop, se_prop)
  )

##---- 3.2 modelling ----
###---- 3.2.1 Model selection ----
bin0<-glmmTMB(cbind(hatched, start_eggs-hatched) ~ poly(Temperature, 2)*species +
                poly(Temperature, 2)*TypeEggs +  TypeEggs*species + Period_length*species ,
              data=cfr.df, family=binomial)

#check for overdispersion in the binomial GLM
overdisp_fun(bin0) # some overdisperison

# Simpler model (check if interaction is needed)
betabin0 <- glmmTMB(cbind(hatched, start_eggs-hatched) ~ poly(Temperature, 2)*species + 
                      poly(Temperature, 2)*TypeEggs +  TypeEggs*species+ Period_length*species ,
                    data=cfr.df, family=betabinomial)
# summary(betabin0)
overdisp_fun(betabin0)
AIC(bin0, betabin0)

betabin1<-glmmTMB(cbind(hatched, start_eggs-hatched) ~ poly(Temperature, 2)*species + 
                    poly(Temperature, 2)*TypeEggs +  TypeEggs*species+ Period_length*species + (1|replicate),
                  data=cfr.df, family=betabinomial)
# summary(betabin1)
AIC(betabin0, betabin1)
anova(betabin0, betabin1) # ok without the random eff
summary(betabin0)

# model diagnostic
res <- simulateResiduals(betabin0)
plot(res)
DHARMa::testDispersion(res)
DHARMa::testZeroInflation(res)
testOutliers(res, type="bootstrap")
testUniformity(res) 

###---- 3.2.2 Model coefficients table ----
# Extract tidy coefficients from glmmTMB model
summary(betabin0)
tab_betabin0 <- tidy(betabin0, effects = "fixed") %>%
  mutate(
    # Compute 95% confidence intervals (on β scale)
    conf.low = estimate - 1.96 * std.error,
    conf.high = estimate + 1.96 * std.error
  ) %>%
  mutate(
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    conf.low = round(conf.low, 3),
    conf.high = round(conf.high, 3),
    p.value = ifelse(p.value < 0.001, "<0.001", round(p.value, 3))
  ) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high)

# Print the table
tab_betabin0

##---- 3.3 get marginal means ----
# Choose specific values to evaluate at
ref_grid <- ref_grid(
  betabin0,
  at = list(
    Temperature    = c(-20, -15, -10, -5, 0),
    Period_length  = c(2,5,10,20,30), #c("2 days", "5 days", "10 days", "20 days", "30 days"),
    species        = c("Aedes albopictus", "Aedes koreicus"),
    TypeEggs       = c("Diapausing", "Non-Diapausing")
  )
)


# Get estimated marginal means by TypeEggs
emm <- emmeans(ref_grid,
               ~ species * TypeEggs | Temperature + Period_length,
               type = "response")

# get two heatmaps
emm_df <- as.data.frame(emm)
species_labs <- c(
  "Aedes albopictus" = "italic('Aedes albopictus')",
  "Aedes koreicus"   = "italic('Aedes koreicus')"
)
p <- emm_df %>% 
  ggplot(aes(x = Temperature, y = as.factor(Period_length), fill = prob)) +
  geom_tile() +
  scale_fill_viridis_c(option = "mako", name = "Proportion of hatched eggs", limits=c(0,1)) +
  labs(
    x = "Temperature (°C)", y = "Cold exposure (days)") +
  facet_grid(species ~ TypeEggs)+
  theme_classic()+
  theme(legend.position = "bottom",
        text = element_text(size=16),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size=16),
        strip.text.y =  element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )
p
outname <- paste0(outdir, "/heatmap_cfr_alboXkor_",Sys.Date(), ".png")
ggsave(p, filename = outname, width = 15, height = 10, device='png', dpi=320)

p <- emm_df %>%
  filter(species=="Aedes albopictus") %>% 
  ggplot(aes(x = Temperature, y = as.factor(Period_length), fill = prob)) +
  geom_tile() +
  scale_fill_viridis_c(option = "mako", name = "Proportion of hatched eggs", limits=c(0,1)) +
  labs(
    x = "Temperature (°C)", y = "Cold exposure (days)") +
  facet_wrap( ~ TypeEggs)+
  theme_classic()+
  theme(legend.position = "bottom",
        text = element_text(size=16),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size=16),
        strip.text.y =  element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )
p
outname <- paste0(outdir, "/heatmap_albo_",Sys.Date(), ".png")
ggsave(p, filename = outname, width = 15, height = 10, device='png', dpi=320)
