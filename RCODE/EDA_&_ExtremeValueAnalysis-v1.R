
library(foreign)
library(quantreg)
library(dplyr)
library(ggplot2)
library(extRemes)



## Load data -----------------------------------------------------------------

fPath <- "./DATA/AA_1990_2015_ICNF_v1.dbf"

firedb <- read.dbf(fPath, as.is = TRUE)
firedb <- firedb[,-4]
colnames(firedb) <- c("year","area")


# Aux funs -------------------------------------------------------------------

countIf <- function(x,thresh=10000) length(x[x > thresh])


## Preliminary analyses ------------------------------------------------------

fire_qhigh <- firedb %>% 
                group_by(year) %>%
                summarise(q90 = quantile(area, probs = 0.90),
                          q95 = quantile(area, probs = 0.95),
                          q99 = quantile(area, probs = 0.99),
                          maxf = max(area)) %>%
                as.data.frame

# Median 99% quantile
median(fire_qhigh$q99)
# Median of maxima
median(fire_qhigh$maxf)


# Linear trend on selected quantiles / maxima

lmq99 <- lm(q99~year, data = fire_qhigh %>% filter(year<=2005))

lmq95 <- lm(q95~year, data = fire_qhigh)

lmq90 <- lm(q90~year, data = fire_qhigh)

# Quantile plots

g <- ggplot(fire_qhigh) + 
     # Quantile 90%
     geom_point(aes(x=year, y=q90),color="dark green") + 
     geom_line(aes(x=year, y=q90),color="dark green") + 
     # Quantile 95%
     geom_point(aes(x=year, y=q95),color="blue") + 
     geom_line(aes(x=year, y=q95),color="blue") + 
     # Quantile 99%
     geom_point(aes(x=year, y=q99),color="red") + 
     geom_line(aes(x=year, y=q99),color="red")

plot(g)     


# Log area plot (2005 - 2015) with quantile regression for (1%, 50% and 99% quantiles)
# Perhaps not the best option for this kind of analyses...
#
g <-  ggplot(firedb %>% filter(year>2005),aes(x=year, y=log10(area))) + 
      geom_point(alpha=0.5) + 
      geom_quantile(quantiles = c(0.01,0.5,0.99))

plot(g)



## Calculate the frequency of wildfires above certain area thresholds ---------



fire_freq <- firedb %>% 
  group_by(year) %>%
  summarise(f100 = countIf(area, thresh=100),
            f250 = countIf(area, thresh=250),
            f500 = countIf(area, thresh=500),
            f1k = countIf(area, thresh=1000),
            f2k = countIf(area, thresh=2000),
            f5k = countIf(area, thresh=5000),
            f7_5k = countIf(area, thresh=7500),
            f10k = countIf(area, thresh=10000)) %>%
  as.data.frame



## Extreme value analysis using an incremental moving window approach --------- 
##
## Uses the block maxima approach (with Generalized Extreme Value distribution) 
## and pre-filters for a minimum area (to consider only very large fires, i.e., 
## fires above 2000ha)
##
## The starting moving window has ten years and it follows until it includes all 
## the available time-series (hence incremental moving window ;-)
##
##

mws <- 10                 # Window size
yr_st <- min(firedb$year) # Start year
yr_en <- max(firedb$year) # End year


i<-0
yr1 <- yr_st

repeat{
  
  # Set the end year
  yr2 <- (yr_st - 1) + mws + i
  
  if(yr2 > yr_en) break

  # Perform GEV fitting (block maxima approach) for fires above 2000ha
   z <- fevd(x = area ~ year, type = "GEV",
           data = firedb %>% filter(area > 2000 & (year >= yr1 & year <= yr2)))
  
  # Perform POT approach with Generalized Pareto Distribution 
  # with a threshold above 2000 ha
  #z <- fevd(x = area ~ year, type = "GP", threshold = 2000,
  #          data = firedb %>% filter(year >= yr1 & year <= yr2))
  
  #summary(z)
  
  rl <- return.level(z,return.period = c(5,10,20,50,100))
  print(rl)
  
  if(i==0)
    TMP <- as.numeric(rl)
  else
    TMP <- rbind(TMP, as.numeric(rl))
  
  #plot(z)
  i <- i + 1

}


## Plot the change in the maximum fire size for the 10-yr return time 

colnames(TMP) <- paste("rp_",c(5,10,20,50,100),sep="")

g<-ggplot(as.data.frame(TMP), aes(x=mws:(mws+nrow(TMP)-1), y= rp_10)) + 
  geom_point(size=3, alpha=0.5) + 
  geom_smooth(color="red") + 
  ylab("Fire size in hectares (10-yrs return time)") + 
  xlab("Analysis window-size (from 1990; in years)")

plot(g)

#ggsave(filename = "./OUT/FireSize10yrReturnTime-POT-GP-v1.png", plot = g)

ggsave(filename = "./OUT/FireSize10yrReturnTime-BlockMax-GEV-v1.png", plot = g)





