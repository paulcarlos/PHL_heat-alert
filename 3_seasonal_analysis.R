Sys.setLanguage("en")

library(splines)
library(dlnm)
library(mixmeta)
library(gnm)
library(lubridate)

source("findmmt.R") # function for minimum mortality temperature by Aurelio Tobias

# load data
dat0 <- readRDS("mort_station-muni_daily_2006-2024.rds") # CANNOT PUBLICLY SHARE
mpred0 <- readRDS("48station_metapred.rds")
thr <- read.csv("table_thresholds_46areas_clust.csv")

# remove Basco & Cuyo
rmsta <- c("Basco","Cuyo")
dat0 <- dat0[!dat0$station%in%rmsta,]
mpred0 <- mpred0[!mpred0$station%in%rmsta,]
mpred <- mpred0[,c("station","lat","lon","elev")]

# add stratum
dat0$year <- year(dat0$date)
dat0$month <- month(dat0$date)
dat0$dow <- wday(dat0$date)

# select years and months
yr <- 2015:2024
mam <- 3:5
jja <- 6:8
son <- 9:11
djf <- c(12,1,2)
dat <- dat0[dat0$year %in% yr,]

# model specs
tkn <- 0.2
lagd <- 1 
fnc <- "ns"; plus <- 1

# objects
muni <- mpred$station
cfmam <- cfjja <- cfson <- cfdjf <- matrix(data=NA,nrow=length(muni),ncol=length(tkn)+plus,dimnames=list(muni))
vcmam <- vcjja <- vcson <- vcdjf <- vector("list",length(muni)); names(vcmam) <- names(vcjja) <- names(vcson) <-names(vcdjf) <- muni

# FIRST STAGE
for (i in 1:nrow(mpred)) {
  cat(i," ")
  s1 <- mpred$station[i]
  
  # subset full
  d1 <- dat0[dat0$station==s1 & dat0$year %in% 2014:2024,]
  d1 <- d1[!grepl(paste0("2014-",sprintf("%02d",1:11),collapse="|"),d1$date),]
  d1$strat <- factor(factor(d1$year):factor(d1$month):factor(d1$dow))
  
  # lag mortality
  d1$lagall <- zoo::rollapply(data=d1$all,width=lagd+1,FUN="mean",align="right",fill=NA,partial=TRUE)
  
  # subset seasons
  dmam <- d1[d1$month %in% mam,]
  dmam$dos <- unname(unlist(tapply(dmam$date,dmam$year,function(z)1:length(z))))
  djja <- d1[d1$month %in% jja,]
  djja$dos <- unname(unlist(tapply(djja$date,djja$year,function(z)1:length(z))))
  dson <- d1[d1$month %in% son,]
  dson$dos <- unname(unlist(tapply(dson$date,dson$year,function(z)1:length(z))))
  ddjf <- d1[d1$month %in% djf,]
  ddjf$year[ddjf$month==12] <- ddjf$year[ddjf$month==12]+1
  ddjf$dos <- unname(unlist(tapply(ddjf$date,ddjf$year,function(z)1:length(z))))
  
  # keep
  dmam$keep <- FALSE
  str1 <- aggregate(lagall~strat,data=dmam,FUN="sum")
  dmam$keep[dmam$strat%in%str1$strat[str1$lagall>0]] <- TRUE
  djja$keep <- FALSE
  str1 <- aggregate(lagall~strat,data=djja,FUN="sum")
  djja$keep[djja$strat%in%str1$strat[str1$lagall>0]] <- TRUE
  dson$keep <- FALSE
  str1 <- aggregate(lagall~strat,data=dson,FUN="sum")
  dson$keep[dson$strat%in%str1$strat[str1$lagall>0]] <- TRUE
  ddjf$keep <- FALSE
  str1 <- aggregate(lagall~strat,data=ddjf,FUN="sum")
  ddjf$keep[ddjf$strat%in%str1$strat[str1$lagall>0]] <- TRUE
  
  # crossbasis
  cbh_mam <- crossbasis(dmam$heatindex,lag=lagd,
                        argvar=list(fun=fnc,knots=quantile(dmam$heatindex,tkn,na.rm=TRUE)),
                        arglag=list(fun="integer"),group=dmam$year)
  cbh_jja <- crossbasis(djja$heatindex,lag=lagd,
                        argvar=list(fun=fnc,knots=quantile(djja$heatindex,tkn,na.rm=TRUE)),
                        arglag=list(fun="integer"),group=djja$year)
  cbh_son <- crossbasis(dson$heatindex,lag=lagd,
                        argvar=list(fun=fnc,knots=quantile(dson$heatindex,tkn,na.rm=TRUE)),
                        arglag=list(fun="integer"),group=dson$year)
  cbh_djf <- crossbasis(ddjf$heatindex,lag=lagd,
                        argvar=list(fun=fnc,knots=quantile(ddjf$heatindex,tkn,na.rm=TRUE)),
                        arglag=list(fun="integer"),group=ddjf$year)
  
  # model and predict
  modmam <- gnm(all ~ cbh_mam + ns(rhum,df=3), data=dmam, family=quasipoisson, eliminate=strat, subset=keep, na.action=na.exclude)
  cenmam <- findmmt(cbh_mam,modmam,by=0.1)
  cpmam <- crossreduce(cbh_mam, modmam, cen=cenmam, by=0.1)
  
  modjja <- gnm(all ~ cbh_jja + ns(rhum,df=3), data=djja, family=quasipoisson, eliminate=strat, subset=keep, na.action=na.exclude)
  cenjja <- findmmt(cbh_jja,modjja,by=0.1)
  cpjja <- crossreduce(cbh_jja, modjja, cen=cenjja, by=0.1)
  
  modson <- gnm(all ~ cbh_son + ns(rhum,df=3), data=dson, family=quasipoisson, eliminate=strat, subset=keep, na.action=na.exclude)
  censon <- findmmt(cbh_son,modson,by=0.1)
  cpson <- crossreduce(cbh_son, modson, cen=censon, by=0.1)
  
  moddjf <- gnm(all ~ cbh_djf + ns(rhum,df=3), data=ddjf, family=quasipoisson, eliminate=strat, subset=keep, na.action=na.exclude)
  cendjf <- findmmt(cbh_djf,moddjf,by=0.1)
  cpdjf <- crossreduce(cbh_djf, moddjf, cen=cendjf, by=0.1)
  
  # extract
  cfmam[s1,] <- coef(cpmam)
  cfjja[s1,] <- coef(cpjja)
  cfson[s1,] <- coef(cpson)
  cfdjf[s1,] <- coef(cpdjf)
  vcmam[[s1]] <- vcov(cpmam)
  vcjja[[s1]] <- vcov(cpjja)
  vcson[[s1]] <- vcov(cpson)
  vcdjf[[s1]] <- vcov(cpdjf)
  
}

# save coefficients
save(tkn,lagd,fnc,cfmam,cfjja,cfson,cfdjf,vcmam,vcjja,vcson,vcdjf,file="coef_vcov_season.rda")


# SECOND STAGE META-ANALYSIS
mpred$clust <- factor(mpred$clust)
mtmam <- mixmeta(formula=cfmam~clust,S=vcmam,data=mpred)
summary(mtmam)$i2[1]
mtjja <- mixmeta(formula=cfjja~clust,S=vcjja,data=mpred)
summary(mtjja)$i2[1]
mtson <- mixmeta(formula=cfson~clust,S=vcson,data=mpred)
summary(mtson)$i2[1]
mtdjf <- mixmeta(formula=cfdjf~clust,S=vcdjf,data=mpred)
summary(mtdjf)$i2[1]

# create predicted coef and vcov by cluster
prdmam <- predict(mtmam,newdata=data.frame("clust"=factor(1:4,clust)),vcov=TRUE)
prdjja <- predict(mtjja,newdata=data.frame("clust"=factor(1:4,clust)),vcov=TRUE)
prdson <- predict(mtson,newdata=data.frame("clust"=factor(1:4,clust)),vcov=TRUE)
prddjf <- predict(mtdjf,newdata=data.frame("clust"=factor(1:4,clust)),vcov=TRUE)



#rm(list=ls());gc()