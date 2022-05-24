### SD Creek Campus Dr rating curve development pre and post tmdl

install.packages("rstatix")

#load libraries
library("ggplot2")
library("dplyr")
library(rstatix)

### read in sediment data
sed <- read.csv("./data/01_input/CampusDr_sedimentratingdata.csv") %>% 
  mutate(log.Q.MCM = log10(AnnualQ_MCM),
         log.sed.tons = log10(Annual_sediment_tons))

#add column for TMDL
sed$tmdl <- "Pre-TMDL (1983-1999)"
sed$tmdl[sed$Year > 1999] <- "Post-TMDL (2000-2018)"
#set levels of tmdl
sed$tmdl <- factor(sed$tmdl, levels = c("Pre-TMDL (1983-1999)", "Post-TMDL (2000-2018)"))

sed$AnnualQ_acreft <- as.numeric(sed$AnnualQ_acreft)



#plot sediment rating curve
sed.rat <- ggplot(sed, aes(x=AnnualQ_MCM, y=Annual_sediment_tons, color = tmdl)) +
  geom_point() +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  geom_smooth(method = 'lm', , se=F) + #
  labs(x = "Annual Discharge (MCM)", y = "Annual Sediment Load (tons)", colour = "Time Period") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position="top")
  
#save jpg
ggsave(sed.rat, file="./data/01_output/SdCkCampus_ratings.jpg", dpi=300, width = 5, height= 5)

#fit linear model to log-log transformed to get a and b
sed.pre <- sed %>% 
  filter(Year < 2000)
#fit linear model to log log
pre.fit <- lm(log.sed.tons ~ log.Q.MCM, sed.pre)
#pre tmdl av annual flow
pre.mean.ann.flow.MCM <- mean(sed.pre$AnnualQ_MCM)
pre.mean.ann.sed.tons <- mean(sed.pre$Annual_sediment_tons)
pre.mean.ann.sed.tons.km2 <- pre.mean.ann.sed.tons/306


#fit linear model to post tmdl
sed.post <- sed %>% 
  filter(Year > 1999)
#fit linear model to log log
post.fit <- lm(log.sed.tons ~ log.Q.MCM, sed.post)
#post tmdl av annual flow
post.mean.ann.flow.MCM <- mean(sed.post$AnnualQ_MCM)
post.mean.ann.sed.tons <- mean(sed.post$Annual_sediment_tons)
post.mean.ann.sed.tons.km2 <- post.mean.ann.sed.tons/306

#ancova test sig diff between two models post and pre tmdl
anova.diff <- anova_test(sed, log.sed.tons ~ log.Q.MCM + tmdl, type=3, detailed = TRUE) #type 3 is ancova


#2005-2019
sed.20052019 <- sed %>% 
  filter(Year > 2004)
#post tmdl av annual flow
post.mean.ann.flow.MCM.20052019 <- mean(sed.20052019$AnnualQ_MCM)
post.mean.ann.sed.tons.20052019 <- mean(sed.20052019$Annual_sediment_tons)

#############################################################################
####read in all long-term gages
sed.all <- read.csv("./data/01_input/Campus_Peters_Culver_sedimentratingdata.csv") %>% 
  mutate(log.Q.MCM = log10(AnnualQ_MCM),
         log.sed.tons = log10(Annual_sediment_tons))

#add column for TMDL
sed.all$tmdl <- "Pre-TMDL (1983-1999)"
sed.all$tmdl[sed.all$Year > 1999] <- "Post-TMDL (2000-2018)"
#set levels of tmdl
sed.all$tmdl <- factor(sed.all$tmdl, levels = c("Pre-TMDL (1983-1999)", "Post-TMDL (2000-2018)"))

sed.all$AnnualQ_acreft <- as.numeric(sed.all$AnnualQ_acreft)

#plot sediment rating curve
sed.rat.all <- ggplot(sed.all, aes(x=AnnualQ_MCM, y=Annual_sediment_tons, color = tmdl)) +
  geom_point() +
  facet_wrap(~site, scales = "free") +
  #scale_x_continuous(trans='log10') +
  #scale_y_continuous(trans='log10') +
  geom_smooth(method = 'lm', se=F) +
  labs(x = "Annual Discharge (MCM)", y = "Annual Sediment Load (tons)", colour = "Time Period") +
  theme(legend.position="top")

#save jpg
ggsave(sed.rat.all, file="./data/01_output/SdCkCampus_ratings.jpg", dpi=300, width = 5, height= 5)

