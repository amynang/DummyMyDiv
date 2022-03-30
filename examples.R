library(tidyverse)
library(glmmTMB)
library(performance)
library(sjPlot)
set.seed(2022)
##### Bait lamina strips #####

# creating a dummy version of MyDiv

# 2 blocks, 40 plots each
dum = data.frame(block = rep(c("1","2"),40),
                 rich = NA,
                 myco = NA)
# in each block 10 replicates of monocultures, 15 of 2 and 4 species
dum$rich[dum$block=="1"] = c(rep(1,10),rep(c(2,4),15))
dum$rich[dum$block=="2"] = c(rep(1,10),rep(c(2,4),15))

# in each block and tree richness level replicates of mycorhiza "treatment"
dum$myco[dum$block=="1" & dum$rich==1] = rep(c("AM","EM"), 5)
dum$myco[dum$block=="1" & dum$rich==2] = rep(c("AM","EM","Both"), 5)
dum$myco[dum$block=="1" & dum$rich==4] = rep(c("AM","EM","Both"), 5)

dum$myco[dum$block=="2" & dum$rich==1] = rep(c("AM","EM"), 5)
dum$myco[dum$block=="2" & dum$rich==2] = rep(c("AM","EM","Both"), 5)
dum$myco[dum$block=="2" & dum$rich==4] = rep(c("AM","EM","Both"), 5)

# check numbers
table(dum$rich,dum$myco)
table(dum$block,dum$myco)
table(dum$block,dum$rich)

# rearrange dataframe
dum = dum %>% arrange(block,rich,myco)

# re-leveling mycorhiza treatment making ectomycorhizas the reference point
dum$myco = factor(dum$myco, levels = c("EM","AM","Both"))

# random scenario 1 - no effect of predictors on response
# here is the distribution of the random data
hist(rbeta(80,5,5))
plot(density(rbeta(80,5,5)))
dum$BLS = rbeta(80,5,5)

# a beta regression mixed model with the two predictors and their interaction
m.bls = glmmTMB(BLS ~ scale(rich) + myco + myco:scale(rich) + (1|block),
                family=beta_family(),
                data = dum)

summary(m.bls)
# model diagnostics
check_model(m.bls)
# model estimates are weird because are given as odds and I am not sure how to change that
tab_model(m.bls)

# but the model predictions look like they are in the right scale
ggplot2::ggplot(dum, (aes(x=rich, y=BLS, color=myco))) +  
  geom_point(position = position_jitter(width = 0.10), alpha = 0.9) + 
  geom_line(aes(y=predict(m.bls, type = "response", dum))) +
  theme(panel.background=element_blank(),
        legend.position="right",
        text = element_text(size=16),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.line.y = element_line(colour = "black", size=.8, lineend="round"),
        axis.line.x = element_line(colour = "black", size=.8, lineend="round")
        )

plot_model(m.bls, 
           show.values = TRUE, 
           value.offset = .3) 

# random scenario 2 - bls dependent on richness and both mycorhiza 
# types together lead to higher activity than individually
dum$BLS = case_when(dum$myco == "AM" ~ dum$rich/10,
                    dum$myco == "EM" ~ dum$rich/10,
                    dum$myco == "Both" ~ (dum$rich+2)/10)
dum$BLS = dum$BLS + rnorm(80,0,.05)


m.bls = glmmTMB(BLS ~ scale(rich) + myco + myco:scale(rich) + (1|block),
                family=beta_family(),
                data = dum)

summary(m.bls)
check_model(m.bls)
tab_model(m.bls)


ggplot2::ggplot(dum, (aes(x=rich, y=BLS, color=myco))) +  
  geom_point(position = position_jitter(width = 0.10), alpha = 0.9) + 
  geom_line(aes(y=predict(m.bls, type = "response", dum))) +
  theme(panel.background=element_blank(),
        legend.position="right",
        text = element_text(size=16),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.line.y = element_line(colour = "black", size=.8, lineend="round"),
        axis.line.x = element_line(colour = "black", size=.8, lineend="round")
  )

plot_model(m.bls, 
           #transform = "plogis",
           show.values = TRUE, 
           value.offset = .3) 

##### Prey dummies #####

# now for the prey dummies
# same as before but 16 times everything so we have one row per dummy
dum = data.frame(block = rep(c("1","2"),40*16),
                 rich = NA,
                 myco = NA)
dum$rich[dum$block=="1"] = c(rep(1,10*16),rep(c(2,4),15*16))
dum$rich[dum$block=="2"] = c(rep(1,10*16),rep(c(2,4),15*16))

dum$myco[dum$block=="1" & dum$rich==1] = rep(c("AM","EM"), 5*16)
dum$myco[dum$block=="1" & dum$rich==2] = rep(c("AM","EM","Both"), 5*16)
dum$myco[dum$block=="1" & dum$rich==4] = rep(c("AM","EM","Both"), 5*16)

dum$myco[dum$block=="2" & dum$rich==1] = rep(c("AM","EM"), 5*16)
dum$myco[dum$block=="2" & dum$rich==2] = rep(c("AM","EM","Both"), 5*16)
dum$myco[dum$block=="2" & dum$rich==4] = rep(c("AM","EM","Both"), 5*16)

table(dum$rich,dum$myco)
table(dum$block,dum$myco)
table(dum$block,dum$rich)

dum = dum %>% arrange(block,rich,myco)

# now we need a plot column
dum$plot = rep(c(1:80), each = 16)
# each plot gets 16 dummies
dum$preydummy = rep(1:16, 80)
# each dummy has two chances to get attacked, from our perspective
# as a random example, let's say 50% probability for each trial
dum$date1 = as.integer(rbernoulli(1280,.5))
dum$date2 = as.integer(rbernoulli(1280,.5))

# let's say 5% of the data are NAs
dum = dum %>% mutate(.data = as_tibble(dum),
                     across(
                       .cols = all_of(c("date1", "date2")),
                       .fns = ~ ifelse(row_number(.x) %in% sample(1:n(), 
                                                                  size = (5 * n(
                       ) / 100)), NA, .x)
                     ))
# so now I can count the number of attack events, "bites" and the trials, removing
# the NAs
# So this means out of the times that we know of, this many times where the dummies attacked
# You are still assuming missingness is random, while in reality missingness is an indication 
# of an attack (bird took the dummy and left) but accounting for that would 
# complicate things a lot...
# but this way at least you don't throw away information
dum.sum = dum %>% group_by(block,rich,myco,plot) %>% 
                  summarise(# number of total attacks per plot minus NAs
                            bites = sum(c(date1,date2), na.rm = TRUE),
                            # 32 trials minus NAs
                            trials = n()*2 - sum(is.na(c(date1,date2))))

# get the bites/trials ratio, as a probability of attack
dum.sum$bites.p = dum.sum$bites/dum.sum$trials

# a binomial model, the response is probability of attack and the weights 
# argument gives the number of trials
# you could also model this as (bites/trials) ~ ... results would be the same
m.prey = glmmTMB(bites.p ~ scale(rich) + myco + myco:scale(rich) + (1|block),
                family="binomial",
                weights = trials,
                data = dum.sum)
summary(m.prey)
check_model(m.prey)
tab_model(m.prey)


ggplot2::ggplot(dum.sum, (aes(x=rich, y=bites.p, color=myco))) +  
  geom_point(position = position_jitter(width = 0.10), alpha = 0.9) + 
  geom_line(aes(y=predict(m.prey, type = "response", dum.sum))) +
  theme(panel.background=element_blank(),
        legend.position="right",
        text = element_text(size=16),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.line.y = element_line(colour = "black", size=.8, lineend="round"),
        axis.line.x = element_line(colour = "black", size=.8, lineend="round"))

plot_model(m.bls, 
           #transform = "plogis",
           show.values = TRUE, 
           value.offset = .3) 
