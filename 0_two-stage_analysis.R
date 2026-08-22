
Sys.setLanguage("en")

library(splines)
library(dlnm)
library(mixmeta)
library(gnm)
library(lubridate)

source("findmmt.R") # function for minimum mortality temperature from Aurelio Tobias

# load data
dat0 <- readRDS("mort_station-muni_daily_2006-2024.rds") # CANNOT PUBLICLY SHARE
mpred0 <- readRDS("48station_metapred.rds")

# remove Cuyo and Basco
dat0 <- dat0[!dat0$station%in%c("Basco","Cuyo"),]
mpred0 <- mpred0[!mpred0$station%in%c("Basco","Cuyo"),]

# add stratum
dat0$year <- year(dat0$date)
dat0$month <- month(dat0$date)
dat0$dow <- wday(dat0$date)

# select years and months
yr <- 2015:2024
mo <- 3:5
dat <- dat0[year(dat0$date) %in% yr,]
mpred <- mpred0[,c("station","lat","lon","elev")]
ag1 <- aggregate(heatindex~station,data=dat[dat$month%in%mo,],FUN="mean")
mpred$hi_mam <- ag1$heatindex[match(mpred$station,ag1$station)]
ag1 <- aggregate(rain~year+station,data=dat,FUN="sum")
ag2 <- aggregate(rain~station,data=ag1,FUN="mean")
mpred$rain <- ag2$rain[match(mpred$station,ag2$station)]
ag1 <- aggregate(all~year+station,data=dat[dat$month%in%mo,],FUN="sum")
ag2 <- aggregate(all~station,data=ag1,FUN="mean")
mpred$mort <- ag2$all[match(mpred$station,ag2$station)]

# model specs
tkn <- 0.2
lagd <- 1 
fnc <- "ns"; plus <- 1

# objects
muni <- mpred$station
coef1 <- matrix(data=NA,nrow=length(muni),ncol=length(tkn)+plus,dimnames=list(muni))
vcov1 <- vector("list",length(muni)); names(vcov1) <- muni

# FIRST STAGE
par(mfrow=c(3,3),mar=c(4,4,3,1),oma=c(1,1,1,1))
for (i in 1:nrow(mpred)) {
  cat(i," ")
  s1 <- mpred$station[i]
  
  # subset full
  d1 <- dat[dat$station==s1,]
  d1$strat <- factor(factor(d1$year):factor(d1$month):factor(d1$dow))
  
  # lag mortality
  d1$lagall <- zoo::rollapply(data=d1$all,width=lagd+1,FUN="mean",align="right",fill=NA,partial=TRUE)
  
  # subset MAM
  sdat <- d1[d1$month%in%mo,]
  
  # keep
  sdat$keep <- FALSE
  str1 <- aggregate(all~strat,data=sdat,FUN="sum")
  sdat$keep[sdat$strat%in%str1$strat[str1$all>0]] <- TRUE
  
  # crossbasis
  cbhi1 <- crossbasis(sdat$heatindex,lag=lagd,
                      argvar=list(fun=fnc,knots=quantile(sdat$heatindex,tkn,na.rm=TRUE)),
                      arglag=list(fun="integer"),group=sdat$year)
  
  # model and predict
  mod1 <- gnm(all ~ cbhi1, data=sdat, family=quasipoisson, eliminate=strat, subset=keep, na.action=na.exclude)
  cen1 <- findmmt(cbhi1,mod1,by=0.1)
  cp1 <- crossreduce(cbhi1, mod1, cen=cen1, by=0.1)
  plot(cp1,main=s1,xlab="heat index",ylab="RR")
  
  # extract
  coef1[s1,] <- coef(cp1)
  vcov1[[s1]] <- vcov(cp1)
  
}

# SECOND STAGE META-ANALYSIS
meta1 <- mixmeta(formula=coef1~1,S=vcov1,data=mpred)
summary(meta1)$i2[1] 

# PLOT
cn <- 1:(length(tkn)+plus)
obhi1 <- onebasis(dat$heatindex,fun=fnc,knots=quantile(dat$heatindex,tkn,na.rm=TRUE))
cen1 <- findmmt(obhi1, coef=meta1$coefficients[cn], vcov=meta1$vcov[cn,cn],by=0.1)
cp1 <- crosspred(obhi1,coef=meta1$coefficients[cn], vcov=meta1$vcov[cn,cn], model.link="log", cen=cen1, by=0.1)
plot(cp1,xlab="heat index",ylab="Relative risk (95%CI)",main="Pooled risk function (2015-2024, Mar-May)")

# BLUP
blup1 <- blup(meta1,vcov=TRUE)
names(blup1) <- rownames(coef1)
par(mfrow=c(3,3),mar=c(4,4,3,1),oma=c(1,1,1,1))
for (i in 1:length(blup1)) {
  #i=1
  s1 <- names(blup1)[i]
  sdat <- dat[dat$station==s1 & dat$month %in% mo,]
  obhi1 <- onebasis(sdat$heatindex,fun=fnc,knots=quantile(sdat$heatindex,tkn,na.rm=TRUE))
  cf <- blup1[[s1]]
  cen1 <- findmmt(obhi1, coef=cf$blup, vcov=cf$vcov, by=0.1)
  cp1 <- crosspred(obhi1,coef=cf$blup, vcov=cf$vcov, model.link="log", cen=cen1, by=0.1)
  plot(cp1,main=s1,xlab="heat index",ylab="RR")
}

# save
save(meta1,blup1,tkn,lagd,file="cf_allcause_blup_HI_MAM_2015-2024.rda")


#rm(list=ls());gc()
