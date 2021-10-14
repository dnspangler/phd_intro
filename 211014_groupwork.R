


# Remember to stick everything in your working directory! use getwd()


####################### Setup

# Load required packages - Uncomment the next lines to install:
#install.packages("tidyverse")
#install.packages("readxl")
#install.packages("lme4")
#install.packages("lmeresampler")
#install.packages("mice")

library(tidyverse)
library(readxl)
library(lme4)
library(lmeresampler)
library(mice)

#Define ids with missing data
missing_ids <- c("16","35")






####################### Load data

# Load weight data from excel weights sheet
weight <- read_excel("Anticancer Agent Study - Data set_clean.xlsx",sheet = "weights") %>%
  mutate(id = factor(id), # Convert numeric id to factor
         group = relevel(factor(group),ref = "Ctrl"), # Set Ctrl as reference category
         new_drug = factor(grepl("_227-3",group)), # Recorde group variable as before
         old_drug = relevel(factor(gsub("_227-3","",group)),ref = "Ctrl"))

# read in data from excel (from the volumes sheet) 
raw_vol <- read_excel("Anticancer Agent Study - Data set_clean.xlsx",
                      sheet = "volumes") %>%
  mutate(id = factor(id))


vol <- raw_vol %>% # Assign the following to the "vol" variable
  pivot_longer(-c(id,group)) %>% # Pivot all data except id and group to long format
  mutate(day = as.numeric(name), # convert day to numeric variable
         log_value = log(value), # Generate log transformed volume value
         id = factor(id), # convert id to factor (categorical) variable
         group = relevel(factor(group),ref = "Ctrl"),
         new_drug = relevel(factor(ifelse(grepl("_227-3",group), # generate new_drug variable based on the presence of the string "_227-3" in the group variable and convert to factor
                                          "227-3",
                                          "Ctrl")),
                            ref = "Ctrl"), 
         old_drug = relevel(factor(gsub("_227-3","",group)), # Recode group variable to remove "_227-3" and convert to factor with "ctrl" as reference category
                            ref = "Ctrl"),
         miss = ifelse(id %in% missing_ids,"solid","dashed"))


# Generate line graph using ggplot package
ggplot(data = vol, # use vol as underlying data
       aes(x=day, # use day as x variable
           y=value, # use volume as y variable
           color=old_drug, # color by old drug
           linetype=new_drug)) + # Set line type (solid/dashed) by new drug
  geom_line(aes(group = id), alpha = 0.6) + # add lines to figure (group by id and set transparency to 0.4)
  geom_smooth(se = F,size = 2) # Add smoothed moving average (loess) to figure (don't draw standard errors and make the line bigger)












####################### Analyze weight

weight %>% # Take weight and stick it in as the data parameter in the following functions:
  group_by(old_drug,new_drug,group) %>% # Group by group for summarization
  summarise(mean = mean(weight,na.rm = T), # Calculate the mean value
            se = sd(weight,na.rm = T)/sqrt(n())) %>% # Calculate standard error of mean
  ggplot(aes(y = mean, # generate plot (using above as underlying data)
             x=group,
             fill = old_drug)) +
  geom_bar(stat = "identity") + # Add bars to plot Using the actual value (identity) instead of some statistic
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se),width = 0.5) # Add error bars at +- SEM

# Fit linear regression model to estimate effects
weight_fit <- lm(weight ~ old_drug + new_drug, data = weight)

summary(weight_fit)


#Include interaction?
weight_fit <- lm(weight ~ old_drug + new_drug,data = weight)
weight_fit_int <- lm(weight ~ old_drug * new_drug,data = weight)

AIC(weight_fit)
AIC(weight_fit_int)
logLik(weight_fit)
logLik(weight_fit_int)

anova(weight_fit,weight_fit_int)

# No.
# ... But apparently average treatment effects are too confusing, so generate some variables to estimate models at different reference levels.

weight <- weight %>%
  mutate(old_drug_CAT = relevel(old_drug,ref="CAT"),
         old_drug_CMF = relevel(old_drug,ref="CMF"))

# Handle missing data

mi_weight <- select(weight,-group,-old_drug,-new_drug) %>% # Create new dataset for imputation
  left_join(filter(vol,day==6),by = "id") %>% # Get tumor volumes at day 6 (last day with full data)
  select(id,weight,group,value,old_drug,old_drug_CAT,old_drug_CMF,new_drug) %>% # Include only needed variables
  mice(m=10, method ="pmm", seed = 101,maxit = 30) # Perform imputation

plot(mi_weight)



mi_weight_complete <- complete(mi_weight, 'long') %>% # make dataset for visualizing imputed data
  mutate(imp = factor(.imp),
         miss = ifelse(id %in% missing_ids,1,0))

# Generate plot of weight estimates
ggplot(data = mi_weight_complete,
       aes(x=group,
           y=value,
           color=old_drug,
           size = miss,
           shape = factor(miss)), # Make imputed points trianges
       size=1) +
  geom_point() +
  guides(size=F)


# Build regression models on imputed datasets (Note: Could also apply rubins rules to get SE based CIs fom a single model)
weight_fit_mi <- with(data = mi_weight, exp = lm(weight ~ old_drug * new_drug)) 

weight_fit_mi_CAT <- with(data = mi_weight, exp = lm(weight ~ old_drug_CAT * new_drug))

weight_fit_mi_CMF <- with(data = mi_weight, exp = lm(weight ~ old_drug_CMF * new_drug))

# Combine results 
pooled <- pool(weight_fit_mi)
pooled_CAT <- pool(weight_fit_mi_CAT)
pooled_CMF <- pool(weight_fit_mi_CMF)

### Get final results (with different contrasts)
weight_ci_table <- summary(pooled, conf.int = TRUE) %>%
  select(term,estimate, lower = `2.5 %`,upper = `97.5 %`)

weight_ci_table_CAT <- summary(pooled_CAT, conf.int = TRUE) %>%
  select(term,estimate, lower = `2.5 %`,upper = `97.5 %`)

weight_ci_table_CMF <- summary(pooled_CMF, conf.int = TRUE) %>%
  select(term,estimate, lower = `2.5 %`,upper = `97.5 %`)

weight_ci_table
weight_ci_table_CAT
weight_ci_table_CMF

####################### Analyze volumes

# Add variables for different contrasts

vol <- vol %>%
  mutate(old_drug_CAT = relevel(old_drug,ref="CAT"),
         old_drug_CMF = relevel(old_drug,ref="CMF"))

# try a simple model using only the final tumor volumes
vol_final <- vol %>%
  filter(day == 10) # Filter to only include the final values on day 10

vol_fit_final <- lm(value ~ old_drug * new_drug, # Fit linear regression model for additive effects of old/new drugs on final tumor volume
              data = vol_final)

summary(vol_fit_final)


# Estimate fancy models

# This model is wrong because it doesn't account for repeated measures
vol_fit <- lm(log_value ~ 1 + day + old_drug:day + new_drug:day, 
                     data = vol)

summary(vol_fit)

# This model accounts for random individual level effects (note the substantially reduced significance of the coefs)
vol_fit_ml <- lmer(log_value ~ day + old_drug:day + new_drug:day +# new_drug:old_drug + # define fixed effects
                   (1 + day | id), # define random effects
                   data = vol)

vol_fit_ml_int <- lmer(log_value ~ day + old_drug:day * new_drug:day + # define fixed effects
                       (1 + day | id), # define random effects
                       data = vol)

# Compare models with / without interaction effect with old drugs
anova(vol_fit_ml,vol_fit_ml_int)

# Again, no reason to report different results between CAT/CMF, but whatever.

vol_fit_ml_int_CAT <- lmer(log_value ~ day + old_drug_CAT:day * new_drug:day + # define fixed effects
                         (1 + day | id), # define random effects
                       data = vol)

vol_fit_ml_int_CMF <- lmer(log_value ~ day + old_drug_CMF:day * new_drug:day + # define fixed effects
                         (1 + day | id), # define random effects
                       data = vol)


# Print model summary
summary(vol_fit_ml)
summary(vol_fit_ml_int)


# Add predicted values from model to dataset
vol$pred = predict(vol_fit_ml,newdata = vol,type = "response")
vol$pred_fe = predict(vol_fit_ml,newdata = vol,type = "response",re.form = NA)

vol$pred_int = predict(vol_fit_ml_int,newdata = vol,type = "response")
vol$pred_int_fe = predict(vol_fit_ml_int,newdata = vol,type = "response",re.form = NA)

vol$resid[!is.na(vol$log_value)] = resid(vol_fit_ml)

# Plot predicted vs. observed values
ggplot(vol,
       aes(x=day,
           y=value,
           color=old_drug)) +
  geom_line(aes(group = id),size=1) + # Add lines for observations grouped by individual
  #geom_smooth(se = F,size = 1) + # uncomment this to compare with smoothed moving average from before
  geom_line(aes(y=exp(pred),
                group = id,
                linetype = miss), alpha = 0.4,size = 2) +
  geom_line(aes(y=exp(pred_fe)),size = 2,color = "black") +
  guides(color = FALSE,
         linetype = FALSE) + # Donät show legend for colors
  facet_grid(old_drug~new_drug) + # divide up observations for old/new drugs into seperate plots next to eachother
  coord_cartesian(ylim = c(0,3))

# Make plot of residuals per day (evaluate deviations from exponential curve)
vol %>%
ggplot(aes(x = day,
           y=resid,
           color = old_drug,
           linetype=new_drug)) +
  geom_point() +
  geom_smooth(se=F,size=1)


# Generate line graph including predicted values
ggplot(data = vol, # use vol as underlying data
       aes(x=day, # use day as x variable
           y=value, # use volume as y variable
           color=old_drug, # color by old drug
           linetype=new_drug)) + # Set line type (solid/dashed) by new drug
  geom_line(aes(group = id)) + # add lines to figure (group by id and set transparency to 0.4)
  geom_line(aes(y=exp(pred),group = id),size = 2, alpha = 0.6) # Add smoothed moving average (loess) to figure (don't draw standard errors and make the line bigger)

# check fitted vs. observed values
ggplot(aes(y=pred,x=log_value,color = old_drug,shape = new_drug),data = vol) +
  geom_point()


# Generate boostrapped replicates
vol_fit_boot <- bootstrap(vol_fit_ml_int, .f = fixef, type = "case", B = 100,resample = c(T,T))
vol_fit_boot_CAT <- bootstrap(vol_fit_ml_int_CAT, .f = fixef, type = "case", B = 100,resample = c(T,T))
vol_fit_boot_CMF <- bootstrap(vol_fit_ml_int_CMF, .f = fixef, type = "case", B = 100,resample = c(T,T))


### Get final results for various contrasts of old_drug
vol_ci_table <- confint(vol_fit_boot,type = "perc") %>% # Extract percentile confidence intervals
  select(-type,-level) %>%
  mutate(across(c(estimate,lower,upper),function(x)round(exp(x),3))) # Back transform to original scale (interpret as percentage change per day)


vol_ci_table_CAT <- confint(vol_fit_boot_CAT,type = "perc") %>% 
  select(-type,-level) %>%  
  mutate(across(c(estimate,lower,upper),function(x)round(exp(x),5))) 


vol_ci_table_CMF <- confint(vol_fit_boot_CMF,type = "perc") %>% 
  select(-type,-level) %>%
  mutate(across(c(estimate,lower,upper),function(x)round(exp(x),4)))

vol_ci_table
vol_ci_table_CAT
vol_ci_table_CMF


