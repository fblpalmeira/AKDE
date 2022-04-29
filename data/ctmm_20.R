library(ctmm)
library(rgdal)

#########################akde.Rmd############################################

y <- read.csv("vinte.txt")
summary(y)

y1 <- as.telemetry(y, proj=CRS("+init=epsg:31984 +proj=utm +zone=24 +south +ellps=GRS80
+ +towgs84=0,0,0,0,0,0,0 +units=m +no_defs "),drop = FALSE) # for importing data, use drop = FALSE to make sure it's a proper list.
summary(y1)

Vinte <- y1$Vinte
plot(Vinte)
title("Vinte (number of locations = 5590)")
compass(cex=2)
extent(Vinte)

M.IID <- ctmm.fit(Vinte) # no autocorrelation timescales and # in general, you should be running ctmm.select here instead of ctmm.fit
m.ouf <- ctmm.guess(Vinte,interactive=FALSE) # automated model guess or # calculate fit guess object
M.OUF <- ctmm.fit(Vinte,m.ouf)

UD0 <- akde(Vinte,M.IID)# Compute akde object for each model
UD2 <- akde(Vinte,M.OUF)
UD2w <- akde(Vinte,M.OUF,weights=TRUE)
# calculate one extent for all UDs
EXT <- extent(list(UD0,UD2,UD2w),level=0.95)

# Plot data with AKDE contours
plot(Vinte,UD=UD0,xlim=EXT$x,ylim=EXT$y)
title(expression("IID KDE"["C"]))
compass(cex=2)
plot(Vinte,UD=UD2,xlim=EXT$x,ylim=EXT$y)
title(expression("OUF AKDE"["C"]))
compass(cex=2)
plot(Vinte,UD=UD2w,xlim=EXT$x,ylim=EXT$y)
title(expression("weighted OUF AKDE"["C"]))
compass(cex=2)

# sampling intervals in hours
col <- "hr" %#% diff(Vinte$t)
# minimum adjacent sampling interval
col <- pmin(c(Inf,col),c(col,Inf))
# sampling intervals under 1.5 hours
col <- (col < 1.5)
# red (low-frequency) or yellow (high-frequency)
col <- grDevices::rgb(1,col,0)

plot(Vinte,UD=UD0,xlim=EXT$x,ylim=EXT$y,col=col)
title(expression("IID KDE"["C"]))
plot(Vinte,UD=UD2,xlim=EXT$x,ylim=EXT$y,col=col)
title(expression("OUF AKDE"["C"]))
plot(Vinte,UD=UD2w,xlim=EXT$x,ylim=EXT$y,col=col)
title(expression("weighted OUF AKDE"["C"]))

# compare the area estimates and effective sample sizes.
summary(UD0)
summary(UD2w)

###########################variogram.Rmd#######################################

SVF <- variogram(Vinte)
level <- c(0.5,0.95) # 50% and 95% CIs
xlim <- c(0,12 %#% "hour") # 0-12 hour window
plot(SVF,xlim=xlim,level=level)
title("zoomed in")
plot(SVF,fraction=0.65,level=level)
title("zoomed out")

## Variogram Fitting the Hard Way
m.iid <- ctmm(sigma=23 %#% "km^2")
m.ou <- ctmm(sigma=23 %#% "km^2",tau=6 %#% "day")
plot(SVF,CTMM=m.iid,fraction=0.65,level=level,col.CTMM="red")
title("Independent and identically distributed data")
plot(SVF,CTMM=m.ou,fraction=0.65,level=level,col.CTMM="purple")
title("Ornstein-Uhlenbeck movement")

m.ouf <- ctmm(sigma=23 %#% "km^2",tau=c(6 %#% "day",1 %#% "hour"))
plot(SVF,CTMM=m.ou,level=level,col.CTMM="purple",xlim=xlim)
title("Ornstein-Uhlenbeck movement")
plot(SVF,CTMM=m.ouf,level=level,col.CTMM="blue",xlim=xlim)
title("Ornstein-Uhlenbeck-F movement")

## Variogram Fitting the Easy Way
variogram.fit(SVF)

### Variogram Error is Autocorrelated
# simulate fake buffalo with the same sampling schedule
willa <- simulate(m.ouf,t=Vinte$t)
plot(willa)
title("simulation")
# now calculate and plot its variogram
SVF2 <- variogram(willa)
plot(SVF2,CTMM=m.ouf,fraction=0.65,level=level,col.CTMM="blue")
title("simulation")

## Maximum Likelihood Fitting the Hard Way 
M.IID <- ctmm.fit(Vinte,m.iid)# Ind. Ident. Distr. (IID). Null model. The IID model is assumed when no autocorrelation time scales are given. Conventional home range estimation (KDE)
summary(M.IID)
M.OU <- ctmm.fit(Vinte,m.ou)# Ornstein-Uhlenbeck (OU). ?? = ??r. The OU process is therefore appropriate for data that lack evidence of directional persistence, but where restricted space use is apparent. 
summary(M.OU)
M.OUF <- ctmm.fit(Vinte,m.ouf)# Ornstein-Uhlenbeck F (OUF.) ?? = {??r, ??v}. The OUFfeatures both correlated velocities and restricted space use
summary(M.OUF)

FITS <- list(IID=M.IID,OU=M.OU,OUF=M.OUF)
summary(FITS)

## Maximum Likelihood Fitting the Easy Way
# CRAN policy limits us to 2 cores
FITZ <- ctmm.select(Vinte,m.ouf,verbose=TRUE,cores=2) # Maximum Likelihood Fitting the Easy Way
summary(FITZ)

# Variograms as a Disgnostic for Maximum Likelihood
plot(SVF,CTMM=FITS,col.CTMM=c("red","purple","blue"),fraction=0.65,level=0.5)

title("zoomed out")
plot(SVF,CTMM=FITS,col.CTMM=c("red","purple","blue"),xlim=xlim,level=0.5)
title("zoomed in")

##############################################################################

### Home-range estimation
home.vinte<-homerange(Vinte,M.IID,method="AKDE")
plot(home.vinte)
compass(cex=2)
summary(home.vinte)

home.vinte2<-homerange(Vinte,M.OU,method="AKDE")
plot(home.vinte2)
compass(cex=2)
summary(home.vinte2)

home.vinte3<-homerange(Vinte,M.OUF,method="AKDE")
plot(home.vinte3)
compass(cex=2)
summary(home.vinte3)

