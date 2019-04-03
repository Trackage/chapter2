## ----echo=FALSE----------------------------------------------------------
options(width = 60)
options(continue=" ")



## ----print=FALSE,echo=FALSE----------------------------------------------

suppressMessages(library(mgcv))
invisible(capture.output(library(deldir)))
invisible(capture.output(library(spatstat)))
library(foreign)
library(lattice)
library(sp)
library(trip)
suppressMessages(library(rgdal))
suppressMessages(library(spatstat))
invisible(capture.output(library(maptools)))
suppressMessages(library(geosphere))



## ----echo=FALSE,print=FALSE,fig=TRUE-------------------------------------
 load("RawArgos.Rdata")

cols <- bpy.colors(length(cols), cutoff.tails = 0.2)

lmat <- matrix(1, 3, 3)
lmat[1,3] <- 2
layout(lmat)

plot(crds, type = "n", axes = FALSE, xlab = NA, ylab = NA, asp =1)
bg.col <- "white"
usr <- par("usr")

rect(usr[1], usr[3], usr[2], usr[4], col = bg.col, border = NA)

plot(world, add = TRUE, col = "grey")
plot(macq, add = TRUE, col = "black")
plot(grat, add = TRUE)
for (i in 2:nrow(crds)) {
	lines(crds[(i-1):i, ], col = cols[i-1], lwd = 2)
}
aa <- legend(x = -1118299.11849875, y = -566890.365414213,
       legend = format(round(seq(min(tr$gmt), max(tr$gmt), length = 8), "days"), "%d %B"),
      	col = bpy.colors(7), lwd = 2, bg = "white", title = "        date (1999)")

r <- matrix(bpy.colors(7, cutoff.tails = 0.2), 7, 1)

rasterImage(r, aa$rect$left + aa$rect$w/12, aa$rect$top - aa$rect$h + aa$rect$h/11,
            aa$text$x[1] - aa$rect$w/12, aa$rect$top - aa$rect$h/7)


box()


text(-1104281, 950000, "Macquarie Island")
text(177346.2,-1354524, "Ross Sea")
op <- par(mar = c(0, 0, 4.1, 2.1))
bb <- bbox(world) + matrix(9e6 * c(1, 1, -1, -1), 2)

plot(0, type = "n", xlab = NA, ylab = NA, xlim = bb[1,], ylim = bb[2,], axes = FALSE)
usr <- par("usr")

rect(usr[1], usr[3], usr[2], usr[4], col = bg.col, border = NA)

plot(world, axes = FALSE, xlab = NA, ylab = NA, xlim = bb[1,], ylim = bb[2,], col = "lightgrey", add = TRUE)
plot(grat10, add = TRUE, col = "grey")
box(lwd = 2)
for (i in 2:nrow(crds)) {
	lines(crds[(i-1):i, ], col = cols[i-1], lwd = 1)
}
par(op)
## crds1 <- coordinates(tr[tr$seal == unique(tr$seal)[3], ])
## #plot(project(crds1, p4), type = "l")
## plot(crds1, type = "l")

## for (i in 2:nrow(crds1)) {
## 	xy <- crds1[(i-1):i, ]
## 	pg4 <- paste("+proj=", "gnom", " +lon_0=", xy[1,1], " +lat_0=", xy[1,2], " +over", sep = "")
## 	xy <- project(xy, pg4)
## 	xy <- project(cbind(seq(xy[1,1], xy[2,1], length = 20), seq(xy[1,2], xy[2,2], length = 20)), pg4, inv = TRUE)

## 	#xy <- project(xy, p4)
## 	lines(xy, col = "red")
## }


## ----print=FALSE,echo=FALSE,fig=TRUE-------------------------------------

load("RawArgos.Rdata")
tr <- tr[tr$seal == "c026", ]

crds <- coordinates(tr)

sf <- speedfilter(tr, max.speed = 12.5)

crds.f <- crds[sf, ]
w <- spTransform(world, CRS("+proj=longlat"))

tr$lq <- ordered(factor(tr$lq, c("Z", "B", "A", "0", "1", "2", "3")))

## update 2019, overlay was long removed from sp
lf <- is.na(over(SpatialPoints(crds, proj4string = CRS(proj4string(w))), as(w, "SpatialPolygons")))
crds.lf <- crds[lf, ]
qf <- tr$lq > "B"

op <- par(mfcol = c(1,2))
## PLOT the raw data, then colours to show filtered track
plot(crds, type = "n", main = "land filter", xlab = "longitude", ylab = "latitude")
plot(w, add = TRUE, col = "grey")
lines(crds, lty = 2)
text(175, -60, paste(sum(!lf), "points \nremoved"))

cols <- bpy.colors(length(cols), cutoff.tails = 0.2)

for (i in 2:nrow(crds.lf)) {
 lines(crds.lf[(i-1):i, ], col = cols[i-1], lwd = 2)
}
crds.qf <- crds[qf, ]
plot(crds, type = "n", main = "class filter", xlab = "longitude", ylab = "latitude")
plot(w, add = TRUE, col = "grey")
lines(crds.lf, lty = 2)

cols <- bpy.colors(length(cols), cutoff.tails = 0.2)
text(175, -60, paste(sum(!qf), "points \nremoved"))
for (i in 2:nrow(crds.qf)) {
 lines(crds.qf[(i-1):i, ], col = cols[i-1], lwd = 2)
}
par(op)


## ----print=FALSE,echo=FALSE,fig=TRUE-------------------------------------

source("iterationsSpeedFilter.R")
ok <- extra.speedfilter(tr, max.speed = 12.5)
t <- speedfilter(tr, max.speed = 12.5, test = TRUE)
tr$ok <- speedfilter(tr, max.speed = 12.5)
ok <- ok[-c(1, 2)]
layout(matrix(c(1,2, 3, 3), 2, 2))
#plot(tr$gmt, ok[[1]], type = "n", log = "y")
for (i in c(1, 5)) { # 7, 9)) {
 ## plot(tr$gmt, t$rms, log = "y", type = "l", col = "grey", ylab = "RMS speed", main = "first iteration")
 ##    abline(h = 12.5, col = "red")
 ##    points(tr$gmt[t$rms > 12.5], t$rms[t$rms > 12.5], pch = 4, col = "red")
 ##    points(tr$gmt[t <= 12.5], t[t <= 12.5], pch = 16, cex = 0.5)

plot(tr$gmt, t$rms, log = "y", type = "l", col = "grey", xlab = "time", ylab = "RMS speed", main = if(i == 1) "second iteration" else "sixth iteration")
    abline(h = 12.5, col = "black", lwd = 2)
    points(tr$gmt[ok[[i]] > 12.5], ok[[i]][ok[[i]] > 12.5], pch = 4, col = "black")
    points(tr$gmt[ok[[i]] <= 12.5], ok[[i]][ok[[i]] <= 12.5], pch = 16, cex = 0.5, col = rgb(0, 0, 0, 0.3))
}


plot(coordinates(tr), type = "n", main = "tenth iteration", xlab = "longitude", ylab = "latitude")
plot(w, add = TRUE, col = "grey")
lines(crds, lty = 2)
crds.sf <- crds[tr$ok, ]

cols <- bpy.colors(length(cols), cutoff.tails = 0.2)
text(175, -60, paste(sum(!tr$ok), "points \nremoved"))
for (i in 2:nrow(crds.sf)) {
 lines(crds.sf[(i-1):i, ], col = cols[i-1], lwd = 2)
}


## ----print=FALSE,echo=FALSE,fig=TRUE-------------------------------------

 load("RawArgos.Rdata")
## tr <- tr[tr$seal == "c026", ]

## tr0 <- tr
## tr1 <- data.frame(Time = tr0$gmt, Lon = coordinates(tr0)[,1], Lat = coordinates(tr0)[,2])

## tr2 <- filter.penSS(tr1,1,iterlim=2000,print.level=1)
##save(tr, tr1, tr2, file = "E:\\mike\\working\\doc\\PhD\\Ch-1-Tag Location Problem\\Fullc026_PenSS.Rdata")
load("Fullc026_PenSS.Rdata")
plot(tr1$Lon,tr1$Lat, type = "n", xlab = "longitude", ylab = "latitude")
plot(spTransform(world, CRS("+proj=longlat")), add = T, col = "grey")
lines(tr1$Lon,tr1$Lat, type = "l", lty = 2)
lines(tr2$Lon,tr2$Lat, type = "l", col = "blue", lwd = 2)
crds <- cbind(tr2$Lon, tr2$Lat)

cols <- bpy.colors(length(cols), cutoff.tails = 0.2)
for (i in 2:nrow(crds)) {
	lines(crds[(i-1):i, ], col = cols[i-1], lwd = 2)
}


## ----print=FALSE,echo=FALSE,fig=TRUE-------------------------------------
par(mfcol = c(2,1))
plot(tr1$Time, tr1$Lon, type = "l", xlab = "time", ylab = "longitude")
lines(tr$gmt, coordinates(tr)[,1], col = "grey")
for (i in 2:nrow(crds)) {
    lines(tr2$Time[(i-1):i], crds[(i-1):i, 1], col = cols[i-1], lwd = 2)

}
plot(tr1$Time, tr1$Lat, type = "l", xlab = "time", ylab = "latitude")
lines(tr$gmt, coordinates(tr)[,2], col = "grey")
for (i in 2:nrow(crds)) {
    lines(tr2$Time[(i-1):i], crds[(i-1):i, 2], col = cols[i-1], lwd = 2)
}
par(op)


## ----print=FALSE,echo=FALSE,fig=TRUE-------------------------------------

load("smalltrip.Rdata")
 source("makesmalltrip.R")

coarse <- c(0.2, 0.2)
g <- GridTopology(c(-0.2, -0.2) + coarse/2, coarse, c(6, 6))
pgs <- as.SpatialPolygons.GridTopology(g)
#g <- GridTopology(c(-0.1, -0.1), c(0.2, 0.2), c(6, 6))
# pgs <- as.SpatialPolygons.GridTopology(g)


 tripGridEG(exact = TRUE, main = "", cex = 0.8, point.cex = 2)
 text(0.06, -0.03, "START")
 text(0.86, 0.17, "END")
 pgs1 <- as.SpatialPolygons.GridTopology(g)
g1 <- tripGrid(tr, grid = g)




## ----print=FALSE,echo=FALSE,fig=TRUE-------------------------------------
fine <- coarse/5
g <- GridTopology(c(-0.2, -0.2) + fine/2, fine, c(30, 30))
pgs <- as.SpatialPolygons.GridTopology(g)

# g <- GridTopology(c(-0.1, -0.1), c(0.2, 0.2)/5, c(28, 28))
# pgs <- as.SpatialPolygons.GridTopology(g)
 tripGridEG(exact = TRUE, main = "", cex = 0.8, point.cex = 1.2, show.text = FALSE)
 text(0.06, -0.03, "START")
 text(0.86, 0.17, "END")
 plot(pgs1[g1$z > 0, ], add = TRUE, lwd = 2)
g2 <- tripGrid(tr, grid = g, method = "density")


## ----echo=FALSE,print=FALSE,fig=TRUE-------------------------------------

#library(maptools)
#gpclibPermit()
#par(mfcol = c(1,2))
load("Ch1-1.RData")

load("Fullc026_PenSS.Rdata")
load("world2-ArcMin3.Rdata")
plot(tr, col = "white", xlim = bbox(tr)[1,] + c(-1, 1) * diff(bbox(tr)[1,])/4, ylim = bbox(tr)[2,] + c(-1, 1) * diff(bbox(tr)[2,])/4)

#plot(nowrapRecenter(spTransform(world, CRS("+proj=longlat"))), add = T, col = "grey")
plot(world, add = TRUE, col = "grey")
points(tr)
box()

coordinates(tr2) <- ~Lon+Lat
tr2$id <- "1"
tr2 <- trip(tr2, c("Time", "id"))
lines(coordinates(tr2), col = "darkgrey", lwd = 3)
p4 <- " +proj=laea +lon_0=174.5 +lat_0=-65.5 +ellps=WGS84 +over"

#p4 <- " +proj=laea +lon_0=170.123 +lat_0=-72.31095 +ellps=WGS84 +over"


ex2 <- coordinates(tr2)[c(1, which.max(spDistsN1(coordinates(tr2), coordinates(tr2)[1,,drop=FALSE], longlat = TRUE))), ]

gc2 <-greatCircle(ex2[1,], ex2[2,])
gc2[gc2[,1] < 0, 1] <-  gc2[gc2[,1] < 0, 1] + 360
gc2 <- gc2[c(180:360, 1:179), ]
lines(gc2, lty = 2, lwd = 2)

#ex <- ex2
#ex[2,] <- c(182.6450, -75.58298)
ex <- coordinates(tr)[c(1, which.max(spDistsN1(coordinates(tr), coordinates(tr)[1,,drop=FALSE], longlat = TRUE))), ]

gc <-greatCircle(ex[1,], ex[2,])
gc[gc[,1] < 0, 1] <-  gc[gc[,1] < 0, 1] + 360
gc <- gc[c(180:360, 1:179), ]
lines(gc, lty = 1, lwd = 2)


## ----echo=FALSE,print=FALSE,fig=TRUE-------------------------------------
proj4string(tr) <- CRS("+proj=longlat")
tr.proj <- spTransform(tr, CRS(p4))
plot(tr.proj, col = "white", xlim = bbox(tr.proj)[1,] + c(-1, 1) * diff(bbox(tr.proj)[1,])/2, ylim = bbox(tr.proj)[2,] + c(-1, 1) * diff(bbox(tr.proj)[2,])/2)
points(tr.proj)
plot(spTransform(world, CRS(p4)), add = TRUE, col = "grey")
box()

lines(project(gc2, p4), lwd = 2, lty = 2)
lines(project(gc, p4), lwd = 2, lty = 1)


## ----echo=TRUE,print=FALSE,fig=FALSE-------------------------------------

## originally this was the readArgos data, but now is this:
## load("E:\\mike\\working\\doc\\PhD\\Thesis-MDSUMNER\\figuredata\\RawArgos.Rdata")

dat <- read.csv("trackfile.csv")
names(dat)  #[1] "long" "lat"  "seal"  "date"    "local"     "lq"
library(sp)
coordinates(dat) <- c("long", "lat")
dat$gmt <- as.POSIXct(strptime(paste(dat$date, dat$local),
                      "%d-%b-%y %H:%M:%S"), tz = "GMT") - 10 * 3600


## ----echo=TRUE,eval=FALSE------------------------------------------------
## library(trip)

## ----echo=TRUE,eval=FALSE------------------------------------------------
## tr <- trip(dat, c("gmt", "seal"))

## ----echo=FALSE,eval=TRUE,results=verbatim-------------------------------
try_out <- try(tr <- trip(dat, c("gmt", "seal")))
cat(try_out)


## ----echo=TRUE,print=FALSE-----------------------------------------------
library(trip)
## read in some Argos data
argosfiles <- list.files(path = "G:/DATA/tracks/blackBrowed/", pattern = ".dat", full.names = TRUE, ignore.case = TRUE)
argosdata <- readArgos(argosfiles[1:3])
summary(argosdata)


## ------------------------------------------------------------------------
dat <- as.data.frame(dat)


## ------------------------------------------------------------------------
dat <- dat[!duplicated(dat), ]
head(dat)


## ------------------------------------------------------------------------
dat <- dat[order(dat$seal, dat$gmt), ]


## ------------------------------------------------------------------------
dat$gmt <- adjust.duplicateTimes(dat$gmt, dat$seal)


## ------------------------------------------------------------------------
coordinates(dat) <- c("long", "lat")
tr <- trip(dat, c("gmt", "seal"))


## ------------------------------------------------------------------------
tr$class <- ordered(tr$lq, c("Z", "B", "A", "0", "1", "2", "3"))
proj4string(tr) <- CRS("+proj=longlat +ellps=WGS84 +over")


## ------------------------------------------------------------------------
tr$ok <- speedfilter(tr, max.speed = 12)
summary(tr)


## ----echo=TRUE,print=FALSE,fig=TRUE--------------------------------------
plot(tr, axes = TRUE, cex = 0.4)
plot(world, add = TRUE, col = "grey")
lines(tr[tr$ok & tr$class > "B", ], lwd = 2,
      col = bpy.colors()[seq(10, 90, length = 4)])


## ----print=FALSE,eval=TRUE,echo=FALSE------------------------------------
## bearing conflict
detach(package:geosphere)

## ------------------------------------------------------------------------

library(argosfilter)
data(seal)
lat <- seal$lat
lon <- seal$lon
dtime <- seal$dtime
lc <- seal$lc
cfilter <- sdafilter(lat, lon, dtime, lc)

seal$sda <- !(cfilter == "removed")


## ------------------------------------------------------------------------
library(argosfilter)
library(sp)
library(trip)
library(maptools)

trackAngle <- function(xy) {
    angles <- abs(c(trackAzimuth(xy), 0) -
                  c(0, rev(trackAzimuth(xy[nrow(xy):1, ]))))
    angles <- ifelse(angles > 180, 360 - angles, angles)
    angles[is.na(angles)] <- 180
    angles
}


## default values used by sdafilter
vmax <-  2
ang <- c(15, 25)
distlim <- c(2500,5000)


coordinates(seal) <- ~lon+lat
proj4string(seal) <- CRS("+proj=longlat +ellps=WGS84")
seal$id <- "seal"


## ------------------------------------------------------------------------
range(diff(seal$dtime))

which(duplicated(seal$dtime))



## ------------------------------------------------------------------------
seal <- seal[!duplicated(seal$dtime), ]

seal.tr <- trip(seal, c("dtime", "id"))


## ------------------------------------------------------------------------
seal.tr$speed.ok <- speedfilter(seal.tr, max.speed = vmax * 3.6)


## ------------------------------------------------------------------------
## distances in metres
dsts <- trackDistance(coordinates(seal.tr)) * 1000
angs <- trackAngle(coordinates(seal.tr))

dprev <- c(0, dsts)
dnext <- c(dsts, 0)

## speed, lc, max distance
ok <- (seal.tr$speed.ok | dprev <= 5000) &  (seal.tr$lc > -9)



## ------------------------------------------------------------------------
seal.tr$filt.row <- 1:nrow(seal.tr)
seal.tr$ok <- rep(FALSE, nrow(seal.tr))

df <- seal.tr

## first subset
df <- df[ok, ]
## distlim and angles, progressively
for (i in 1:length(distlim)) {
    dsts <- trackDistance(coordinates(df)) * 1000
    angs <- trackAngle(coordinates(df))
    dprev <- c(0, dsts)
    dnext <- c(dsts, 0)
    ok <- (dprev <= distlim[i] | dnext <= distlim[i])  | angs > ang[i]
    ok[c(1:2, (length(ok)-1):length(ok))] <- TRUE
    df <- df[ok, ]
    ok <- rep(TRUE, nrow(df))
}


## ------------------------------------------------------------------------

seal.tr$ok[ match(df$filt.row, seal.tr$filt.row)] <- ok

sum(seal.tr$sda)
table(seal.tr$sda, seal.tr$lc)
sum(seal.tr$ok)
table(seal.tr$ok, seal.tr$lc)



## ----echo=FALSE,print=FALSE,fig=TRUE-------------------------------------
plot(seal.tr[seal.tr$sda, ], cex = 0.4)
lines(seal.tr, col = rgb(0.5, 0.5, 0.5, 0.3))
lines(seal.tr[seal.tr$sda, ], col = "#0000FF22", lwd = 15)
lines(seal.tr[seal.tr$ok, ], col = "#FF758AFF")
degAxis(1)
degAxis(2)
box()



## ----print=FALSE,echo=FALSE,eval=TRUE------------------------------------
  load('tripGridDEMO.Rdata')


## ----echo=TRUE,print=FALSE,fig=TRUE--------------------------------------
trg  <- tripGrid(tr[tr$ok, ])
image(trg, col = oc.colors(100), axes = TRUE)


## ----print=TRUE,eval=FALSE-----------------------------------------------
## names(trg) <- 'grid.filtered'
## gt <- getGridTopology(trg)
## trg$grid.c026  <- tripGrid(tr[tr$ok & tr$seal == "c026", ], grid = gt)$z
## trg$kde1.filtered <- tripGrid(tr[tr$ok, ], grid = gt, method = "density", sigma = 0.1)$z
## trg$kde3.filtered <- tripGrid(tr[tr$ok, ], grid = gt, method = "density", sigma = 0.3)$z
## 
## for (col.name in names(trg)) trg[[col.name]] <- trg[[col.name]]/3600


## ----print=FALSE,echo=FALSE,eval=TRUE------------------------------------
## save(trg, file = 'saved.fourgrids.Rdata')
load('saved.fourgrids.Rdata')

## ----fig=TRUE------------------------------------------------------------
require(lattice)
library(trip)
trellis.par.set("regions", list(col=oc.colors(256)))
print( spplot(trg))


## ----eval=TRUE,print=FALSE,echo=FALSE,fig=FALSE--------------------------
tripTransform <- function(x, p4) {
	tor <- getTORnames(x)
	d <- as(x, "SpatialPointsDataFrame")
	require(rgdal)

	d <- spTransform(d, p4)
	trip(d, tor)
	}
proj4string(tr) <- CRS("+proj=longlat +ellps=WGS84")
p4 <- CRS("+proj=laea +lon_0=174.5 +lat_0=-65.5 +units=km")

ptr <- tripTransform(tr, p4)

gt1 <- makeGridTopology(ptr, cellsize = c(80, 60))

grd2 <- tripGrid(ptr, grid = gt1)
image(grd2, col = oc.colors(256), axes = TRUE)

usr <- par("usr")



## ----echo=TRUE,print=FALSE,fig=TRUE--------------------------------------
proj4string(tr) <- CRS("+proj=longlat +ellps=WGS84")
p4 <- CRS("+proj=laea +lon_0=174.5 +lat_0=-65.5 +units=km")

ptr <- tripTransform(tr, p4)

gt <- makeGridTopology(ptr, c(50, 50))
gt1 <- makeGridTopology(ptr, cellsize = c(80, 60))

grd2 <- tripGrid(ptr, grid = gt1)
image(grd2, col = oc.colors(256), axes = TRUE)


library(maptools)
data(wrld_simpl)

plot(spTransform(wrld_simpl, p4), add = TRUE, col = "grey")



