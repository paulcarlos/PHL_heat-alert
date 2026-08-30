Sys.setLanguage("en")

library(splines)
library(dlnm)
library(gnm)
library(lubridate)
library(zoo)
library(earth)
library(partykit)
library(primr)

# load data and risk function
load("cf_allcause_blup_HI_MAM_2015-2024.rda")
rdat <- readRDS("mort_station-muni_daily_2006-2024.rds") # CANNOT PUBLICLY SHARE
rdat <- rdat[rdat$station%in%names(blup1),]

# round heat index
rdat$hi <- round(rdat$heatindex)

# add time variables
rdat$year <- year(rdat$date)
rdat$month <- month(rdat$date)

# subset data
yr <- 2015:2024
mo <- 3:5
dat <- rdat[rdat$year%in%yr & rdat$month%in%mo,]

# MOB rules
a <- c(0.1,0.5,0.8)
t <- c(0.1,0.05,0.02)
mobcon <- data.frame("alpha"=rep(a,each=length(t)),
                     "trim"= rep(t,times=length(a)),
                     "thr"=NA)
rm(a,t)

# MARS rules
e <- c(3,5,10)
th <- c(0.001,0.0001,0)
nk <- c(30,50,100)
deg <- 1:2
marscon <- data.frame("endspan"=rep(e,each=length(th)*length(nk)*length(deg)),
                      "thres"=rep(th,each=length(nk)*length(deg),times=length(e)),
                      "nk"=rep(nk,each=length(deg),times=length(e)*length(th)),
                      "deg"=rep(deg,times=length(e)*length(th)*length(nk)))
sum(duplicated(marscon))

# dataframe to save thresholds
sta <- unique(rdat$station)
df <- data.frame("station"=sta,
                 "thr_95p"=NA,"thr_25pct"=NA,"thr_mob"=NA,"thr_mars"=NA,"thr_prim"=NA)

# LOOP THRESHOLD SEARCH
for (i in 1:nrow(df)) {
  cat(i," ")
  s1 <- df$station[i]
  
  # subset station data
  d1 <- rdat[rdat$station==s1,]
  
  # lag mortality, moving average
  d1$lmort <- zoo::rollapply(data=d1$all,width=lagd+1,FUN="mean",align="right",fill=NA)
  
  # get BLUP coefficients and vcov
  b1 <- blup1[[s1]]
  
  # subset data within March-May
  d2 <- d1[d1$month %in% mo,]
  d2$dos <- unname(unlist(tapply(d2$date,d2$year,function(z)1:length(z))))
  
  # derive attributable deaths
  cen <- min(d2$hi,na.rm=TRUE)
  argvar <- list(fun="ns",knots=quantile(d2$hi,tkn,na.rm=TRUE),Bound=range(d2$hi,na.rm=TRUE))
  bvar <- do.call(onebasis,c(list(x=d2$hi),argvar))
  cp <- crosspred(bvar,coef=b1$blup, vcov=b1$vcov, model.link="log", cen=cen, by=1)
  d2$attrd <- as.vector((1-exp(-bvar%*%b1$blup)))*d2$lmort
  
  # find temperature of 95th percentile heat mortality
  df$thr_95p[df$station==s1] <- round(mean(d2$hi[d2$attrd>=quantile(d2$attrd,0.95,na.rm=TRUE)],na.rm=TRUE))
  
  # find temperature that covers 25% of heat-related mortality
  ag1 <- aggregate(attrd~hi,data=d2,FUN="sum")
  ag1 <- ag1[order(ag1$hi,decreasing = TRUE),]
  ag1$cumul <- cumsum(ag1$attrd)
  ag1$prop <- ag1$cumul/sum(ag1$attrd)
  df$thr_25pct[df$station==s1] <- min(ag1$hi[ag1$prop<=0.25],na.rm=TRUE)
  
  # SUBSET DATA REMOVE NAs
  sdat <- d2[!is.na(d2$hi),]
  
  
  # Model-based Recursive Partitioning (MOB)
  mb <- mobcon
  for (b in 1:nrow(mb)) {
    mobtr1 <- glmtree(all ~ year + ns(dos,df=3) | hi, data=sdat, family = quasipoisson(),
                      alpha=mb$alpha[b], bonferroni=FALSE, trim=mb$trim[b],  restart=FALSE, numsplit="left")
    nodepath1 <- partykit:::.list.rules.party(mobtr1)
    inner_path1 <- unlist(strsplit(nodepath1, " & "))
    rules <- as.numeric(gsub("\\D","",c(inner_path1)))
    mb$thr[b] <- round(max(rules,na.rm=TRUE))
  }
  
  if (sum(mb$thr>0,na.rm=TRUE)>0) {
    df$thr_mob[df$station==s1] <- max(mb$thr,na.rm=TRUE)
  } 
  rm(rules)
  
  
  # Multivariate Adaptive Regression Splines (MARS)
  ms <- marscon
  ms$val <- NA
  for (b in 1:nrow(marscon)) {
    #b=1
    marsfit1 <- earth(all ~ hi + year + ns(dos,df=3), data = sdat, glm = list(family = "quasipoisson"), pmethod="backward",
                      endspan = ms$endspan[b], thresh=ms$thres[b], nk=ms$nk[b], degree=ms$deg[b])
    marsv <- c(marsfit1$cuts[,"hi"])
    ms$val[b] <- round(max(marsv,na.rm=TRUE))
  }
  
  if (sum(ms$val,na.rm=TRUE)>0) {
    df$thr_mars[df$station==s1] <- max(ms$val,na.rm=TRUE)
  }
  
  
  # Patient Rule-Induction Method (PRIM)
  peelres <- peeling(sdat$all, sdat$hi, 
                     beta.stop = 10/nrow(sdat), peeling.side = -1,
                     obj.fun = function(y, x, inbox){
                       y <- y[inbox]
                       x <- x[inbox,]
                       dat <- data.frame(y, x, sdat[inbox, c("dos", "year")])
                       fit <- glm(y ~ ns(dos,df=3) + year, data = dat, family = "quasipoisson")
                       pred <- coef(fit)[2]
                       return(exp(pred))
                     })
  primres <- pasting(peelres, support = 0.05)
  df$thr_prim[df$station==s1] <- round(sapply(extract.box(primres)$limits[[1]], "[", 1))
  
}

# SAVE
saveRDS(df,"table_thresholds_46areas.rds")


#rm(list=ls());gc()
