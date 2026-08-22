Sys.setLanguage("en")

library(gnm)
library(lubridate)
library(MatchIt)

# load data
rdat <- readRDS("station_hi-forecast_mort_2024.rds")

# remove stations with lowest data
lowd <- c("Basco","Cuyo")
dat <- rdat[!rdat$station%in%lowd,]

# add month and day of week
dat$month <- month(dat$date)
dat$dow <- wday(dat$date)



################ FALSE NEGATIVES #######################
# create data
d2 <- dat[dat$catobs=="Danger" & dat$catfrc=="Extreme Caution" & abs(dat$hiobs-dat$hifrc)>=2,]
d2$type <- 1 # false negative
d1 <- dat[dat$catobs=="Danger" & dat$catfrc=="Danger" & abs(dat$hiobs-dat$hifrc)<2 & dat$station %in% d2$station,]
d1$type <- 0
d3 <- rbind(d1,d2)

# match
d4 <- matchit(type ~ hiobs + month + station, data=d3, distance = "glm", caliper=0.2, exact = c("station"), link="logit")
d5 <- match.data(d4,subclass="sub")
d5 <- d5[order(d5$sub),]
d5$station <- factor(d5$station)

# subset data with less than 3 days differences
dif <- tapply(d5$date,d5$sub,function(z)abs(interval(z[1],z[2])/days(1)))
sum(dif<=3)
summary(dif)
d6 <- d5[d5$sub%in% which(dif>3),]

# model using GLM
mod1 <- glm(exmort ~ factor(type) + station + hiobs, data=d6, family=quasipoisson)
summary(mod1)
exp(coef(mod1)[2])
exp(confint(mod1)[2,])



################ FALSE POSITIVES #######################
# create data
d1 <- dat[dat$catobs=="Extreme Caution" & dat$catfrc=="Extreme Caution" & abs(dat$hiobs-dat$hifrc)<2,]
d1$type <- 0
d2 <- dat[dat$catobs=="Extreme Caution" & dat$catfrc=="Danger" & abs(dat$hiobs-dat$hifrc)>=2,]
d2$type <- 1 # false positives
d3 <- rbind(d1,d2)

# match
d4 <- matchit(type ~ hiobs + month + station, data=d3, distance = "glm", caliper=0.2, exact = c("station"), link="logit")
d5 <- match.data(d4,subclass="sub")
d5 <- d5[order(d5$sub),]
d5$station <- factor(d5$station)

# check date differences
dif <- tapply(d5$date,d5$sub,function(z)abs(interval(z[1],z[2])/days(1)))
sum(dif<=3)
summary(dif)
d6 <- d5[d5$sub%in% which(dif>3),]

# model using GLM
mod1 <- glm(exmort ~ factor(type) + station + hiobs, data=d6, family=quasipoisson)
summary(mod1)
exp(coef(mod1)[2])
exp(confint(mod1)[2,])



################ HIGHER WARNING LEVEL #######################
# create data
d2 <- dat[dat$catobs=="Danger" & dat$catfrc=="Danger" & abs(dat$hiobs-dat$hifrc)<2,]
d2$type <- 1 
d1 <- dat[dat$catobs=="Extreme Caution" & dat$catfrc=="Extreme Caution" & abs(dat$hiobs-dat$hifrc)<2 & dat$station %in% d2$station,]
d1$type <- 0
d3 <- rbind(d1,d2)

# match
d4 <- matchit(type ~ month + station, data=d3, distance = "glm", caliper=0.2, exact = c("station"), link="logit")
d5 <- match.data(d4,subclass="sub")
d5 <- d5[order(d5$sub),]
d5$station <- factor(d5$station)

# check date differences
dif <- tapply(d5$date,d5$sub,function(z)abs(interval(z[1],z[2])/days(1)))
sum(dif<=3)
summary(dif)
d6 <- d5[d5$sub%in% which(dif>3),]

# model using GLM
mod1 <- glm(exmort ~ factor(type) + station, data=d6, family=quasipoisson)
summary(mod1)
exp(coef(mod1)[2])
exp(confint(mod1)[2,])

#rm(list=ls());gc()