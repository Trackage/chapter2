extra.speedfilter <- function (x, max.speed = NULL, test = FALSE)
{

extrares <- list()
    if (!is(x, "trip"))
        stop("only trip objects supported")
    projected <- is.projected(x)
    if (is.na(projected)) {
        projected <- FALSE
        warning("coordinate system is NA, assuming longlat . . .")
    }
    FUN <- function(x, aadBUG = FALSE) {
        sqrt(sum((x)^2, na.rm = FALSE)/(if (aadBUG) 1 else length(x)))
    }
    if (is.null(max.speed)) {
        print("no max.speed given, nothing to do here")
        return(x)
    }
    longlat <- !projected
    coords <- coordinates(x)
    time = x[[x@TOR.columns[1]]]
    id = factor(x[[x@TOR.columns[2]]])
    x <- coords[, 1]
    y <- coords[, 2]
    aadBUG = FALSE
    pprm <- 3
    grps <- levels(id)
    if (length(x) != length(y))
        stop("x and y vectors must be of same\nlength")
    if (length(x) != length(time))
        stop("Length of times not equal to number of points")
    okFULL <- NULL


    if (test)
        res <- list(speed = numeric(0), rms = numeric(0))
    for (sub in grps) {
        ind <- id == sub
        xy <- matrix(c(x[ind], y[ind]), ncol = 2)
        tms <- time[ind]
        npts <- nrow(xy)
        if (pprm%%2 == 0 || pprm < 3)
            stop("Points per running mean should be odd and greater than 3, pprm = 3")
        RMS <- rep(max.speed + 1, npts)

        offset <- pprm - 1
        ok <- rep(TRUE, npts)
        if (npts < (pprm + 1)) {
            warning("Not enough points to filter ID: \"", sub,
                "\"\n continuing . . . \n")
            okFULL <- c(okFULL, ok)
            next
        }
        index <- 1:npts
        iter <- 1
        while (any(RMS > max.speed, na.rm = TRUE)) {
            n <- length(which(ok))
            speed1 <- OLDtrackDistance(xy[ok, ], longlat = longlat)/(diff(unclass(tms[ok]))/3600)
            speed2 <- OLDtrackDistance(xy[ok, ], longlat = longlat,
                push = 2)/((unclass(tms[ok][-c(1, 2)]) - unclass(tms[ok][-c(n -
                1, n)]))/3600)
            thisIndex <- index[ok]
            npts <- length(speed1)
            if (npts < pprm) {
                next
            }
            sub1 <- rep(1:2, npts - offset) + rep(1:(npts - offset),
                each = 2)
            sub2 <- rep(c(0, 2), npts - offset) + rep(1:(npts -
                offset), each = 2)
            rmsRows <- cbind(matrix(speed1[sub1], ncol = offset,
                byrow = TRUE), matrix(speed2[sub2], ncol = offset,
                byrow = TRUE))
            RMS <- c(rep(0, offset), apply(rmsRows, 1, FUN, aadBUG = aadBUG))
if (iter == 1) rms <- RMS
            if (test & iter == 1) {
                res$speed <- c(res$speed, 0, speed1)
                res$rms <- c(res$rms, 0, RMS)
                break
            }
            iter <- iter + 1

            bad <- RMS > max.speed
            bad[is.na(bad)] <- FALSE
            segs <- cumsum(c(0, abs(diff(bad))))
            segs[RMS <= max.speed] <- NA
            peaks <- tapply(RMS, segs, which.max)
            for (i in levels(as.factor(segs))) {
                RMS[segs == i & !is.na(segs)][peaks[[i]]] <- NA
            }
            RMS[1] <- 0
            RMS[length(RMS)] <- 0

		rms[thisIndex] <- c(0, RMS)
		rms[-thisIndex] <- 0
extrares[[iter]] <-  rms

            ok[thisIndex][is.na(RMS)] <- FALSE
        }
        okFULL <- c(okFULL, ok)
    }
    if (test)
        return(res)
    filt <- okFULL
    filt
return(extrares)
}


OLDtrackDistance <-
function (track, longlat = FALSE, push = 1)
{

	#print(longlat)
	#track <- coordinates(spData)
    if (!is.matrix(track))
        stop("track must be two-column matrix")
    if (ncol(track) != 2)
        stop("track must be two-column matrix")
    n1 <- nrow(track) - 1
    if (n1 < 1)
        stop("less than two points")
    res <- numeric(n1 - push + 1)
    for (i in seq(along = res))
      res[i] <- spDistsN1(track[i,,drop = FALSE],
                          track[(i + push),,drop = FALSE ], longlat = longlat)
    res
}

