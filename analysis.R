


# Remember to stick everything in your working directory! use getwd()


# Setup -------------

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






# Load data -----------

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





# Analyze weight -----------

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
weight_fit <- lm(weight ~ group, data = weight)
weight_fit_CAT <- lm(weight ~ relevel(group,ref = "CAT"), data = weight)
weight_fit_CMF <- lm(weight ~ relevel(group,ref = "CMF"), data = weight)

easy_weight_coefs <- rbind(cbind(summary(weight_fit)$coef,
                                 as.data.frame(confint(weight_fit))),
                           cbind(summary(weight_fit_CAT)$coef,
                                  as.data.frame(confint(weight_fit_CAT))),
                           cbind(summary(weight_fit_CMF)$coef,
                                  as.data.frame(confint(weight_fit_CMF))))


#Include interaction?
weight_fit_noint <- lm(weight ~ old_drug + new_drug,data = weight)
weight_fit_int <- lm(weight ~ old_drug * new_drug,data = weight)

AIC(weight_fit_noint)
AIC(weight_fit_int)
logLik(weight_fit_noint)
logLik(weight_fit_int)

anova(weight_fit_noint,weight_fit_int)

# No, but apparently average treatment effects are confusing...

# Handle missing data - Tumor weight seems like it would be highly correlated with volume, 
# so should be able to estimate what the tumor sizes would have been if the mouse had lived.. 
# We'll impute weight based on tumor volumes on the last day all rats were alive (day 6).

mi_weight <- select(weight,-group,-old_drug,-new_drug) %>% # Create new dataset for imputation
  left_join(filter(vol,day==6),by = "id") %>% # Get tumor volumes at day 6
  select(id,weight,group,value,old_drug,new_drug) %>% # Include only needed variables
  mice(m=50, method ="pmm", seed = 101,maxit = 20) # Perform imputation

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
  guides(size="none")


# Build regression models on imputed datasets (Note: Could also apply rubins rules to get SE based CIs fom a single model)

# A bit of bad form repeating code - Should be done with functions and apply(), but I'm lazy

weight_fit_mi <- with(data = mi_weight, 
                      exp = lm(weight ~ group))
weight_fit_mi_CAT <- with(data = mi_weight, 
                          exp = lm(weight ~ relevel(group,ref = "CAT")))
weight_fit_mi_CMF <- with(data = mi_weight, 
                          exp = lm(weight ~ relevel(group,ref = "CMF")))

# Combine results 
pooled <- pool(weight_fit_mi)
pooled_CAT <- pool(weight_fit_mi_CAT)
pooled_CMF <- pool(weight_fit_mi_CMF)

### Get final results (with different contrasts)
weight_ci_table <- summary(pooled, conf.int = TRUE) %>%
  select(term,estimate, lower = `2.5 %`,upper = `97.5 %`) %>%
  mutate(ref = rep("Ctrl",nrow(.)))

weight_ci_table_CAT <- summary(pooled_CAT, conf.int = TRUE) %>%
  select(term,estimate, lower = `2.5 %`,upper = `97.5 %`) %>%
  mutate(ref = rep("CAT",nrow(.)))

weight_ci_table_CMF <- summary(pooled_CMF, conf.int = TRUE) %>%
  select(term,estimate, lower = `2.5 %`,upper = `97.5 %`) %>%
  mutate(ref = rep("CMF",nrow(.)))

weight_ci <- rbind(weight_ci_table,weight_ci_table_CAT,weight_ci_table_CMF)


# Analyze volumes --------------


# Generate line graph using ggplot package
ggplot(data = vol, # use vol as underlying data
       aes(x=day, # use day as x variable
           y=value, # use volume as y variable
           color=old_drug, # color by old drug
           linetype=new_drug)) + # Set line type (solid/dashed) by new drug
  geom_line(aes(group = id), alpha = 0.6) + # add lines to figure (group by id and set transparency to 0.4)
  geom_smooth(se = F,size = 2) # Add smoothed moving average (loess) to figure (don't draw standard errors and make the line bigger)

#Add variables for different model contrasts

vol <- vol %>%
  mutate(group_CAT = relevel(group,ref="CAT"),
         group_CMF = relevel(group,ref="CMF"))

# try a simple model using only the final tumor volumes
vol_final <- vol %>%
  filter(day == 10) # Filter to only include the final values on day 10

vol_fit_final <- lm(value ~ group, # Fit linear regression model for additive effects of old/new drugs on final tumor volume
              data = vol_final)

# Fit seperate models to get contrasts - I think there are packages out there to get contrasts from a single model... But this should work.
vol_fit_final_CAT <- lm(value ~ group_CAT, data = vol_final)
vol_fit_final_CMF <- lm(value ~ group_CMF, data = vol_final)


easy_vol_coefs <- rbind(cbind(summary(vol_fit_final)$coef,
                                 as.data.frame(confint(vol_fit_final))),
                           cbind(summary(vol_fit_final_CAT)$coef,
                                 as.data.frame(confint(vol_fit_final_CAT))),
                           cbind(summary(vol_fit_final_CMF)$coef,
                                 as.data.frame(confint(vol_fit_final_CMF))))

# On to daily volume data... This model is wrong because it doesn't account for repeated measures
vol_fit <- lm(log_value ~ day + group:day, 
              data = vol)

summary(vol_fit)

# Estimate fancy models ---------

# This model accounts for random individual level effects (note the substantially reduced significance of the coefs)
vol_fit_ml_noint <- lmer(log_value ~ day + old_drug:day + new_drug:day + # define fixed effects
                   (1 + day | id), # define random effects
                   data = vol)

vol_fit_ml_int <- lmer(log_value ~ day + old_drug:day * new_drug:day + # add multiplicative term (*) here
                       (1 + day | id),
                       data = vol) # note a term gets dropped (missing Ctrl_227-3 factor I believe)

# Compare models with / without interaction effect with old drugs
anova(vol_fit_ml_noint,vol_fit_ml_int)

# Again, no reason to report different results between CAT/CMF, but whatever. Use groups instead then.

vol_fit_ml <- lmer(log_value ~ day + group:day + # define fixed effects
                         (1 + day | id), # define random effects
                       data = vol)

vol_fit_ml_CAT <- lmer(log_value ~ day + group_CAT:day + # define fixed effects
                             (1 + day | id), # define random effects
                           data = vol)

vol_fit_ml_CMF <- lmer(log_value ~ day + group_CMF:day + # define fixed effects
                         (1 + day | id), # define random effects
                       data = vol)


# Print model summary
summary(vol_fit_ml)


# Add predicted values from model to dataset
vol$pred = predict(vol_fit_ml,newdata = vol,type = "response")
vol$pred_fe = predict(vol_fit_ml,newdata = vol,type = "response",re.form = NA)

vol$pred_int = predict(vol_fit_ml,newdata = vol,type = "response")
vol$pred_int_fe = predict(vol_fit_ml,newdata = vol,type = "response",re.form = NA)

vol$resid[!is.na(vol$log_value)] = resid(vol_fit_ml)

# Plot predicted vs. observed values
ggplot(vol,
       aes(x=day,
           y=log_value,
           color=old_drug)) +
  geom_line(aes(group = id),size=1) + # Add lines for observations grouped by individual
  #geom_smooth(se = F,size = 1) + # uncomment this to compare with smoothed moving average from before
  geom_line(aes(y=pred,
                group = id,
                linetype = miss), alpha = 0.4,size = 2) +
  geom_line(aes(y=pred_fe),size = 2,color = "black") +
  guides(color = FALSE,
         linetype = FALSE) + # Donät show legend for colors
  facet_grid(old_drug~new_drug) #+ # divide up observations for old/new drugs into seperate plots next to eachother
  #coord_cartesian(ylim = c(0,3))

# Make plot of residuals per day (evaluate deviations from exponential curve)
vol %>%
ggplot(aes(x = day,
           y=resid,))+
           #color = old_drug,
           #linetype=new_drug)) +
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
vol_fit_boot <- bootstrap(vol_fit_ml, .f = fixef, type = "case", B = 100,resample = c(T,T))
vol_fit_boot_CAT <- bootstrap(vol_fit_ml_CAT, .f = fixef, type = "case", B = 100,resample = c(T,T))
vol_fit_boot_CMF <- bootstrap(vol_fit_ml_CMF, .f = fixef, type = "case", B = 100,resample = c(T,T))


### Get final results for various contrasts of old_drug
vol_ci_table <- confint(vol_fit_boot,type = "perc") %>% # Extract percentile confidence intervals
  select(-type,-level) %>% # remove unneeded variables
  mutate(across(c(estimate,lower,upper),function(x)round(exp(x),3)),# Back transform coefs to original scale
         ref = rep("Ctrl",nrow(.)))

vol_ci_table_CAT <- confint(vol_fit_boot_CAT,type = "perc") %>% 
  select(-type,-level) %>%  
  mutate(across(c(estimate,lower,upper),function(x)round(exp(x),3)),
         ref = rep("CAT",nrow(.)))


vol_ci_table_CMF <- confint(vol_fit_boot_CMF,type = "perc") %>% 
  select(-type,-level) %>%
  mutate(across(c(estimate,lower,upper),function(x)round(exp(x),3)),
         ref = rep("CMF",nrow(.)))


volume_ci <- rbind(vol_ci_table,
                   vol_ci_table_CAT,
                   vol_ci_table_CMF)
