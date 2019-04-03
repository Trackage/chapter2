 ## #source("E:\\build\\spatstatold\\spatstat\\R\\density.psp.R")
 ## load("E:\\mike\\working\\doc\\PhD\\Ch-1-Tag Location Problem\\Ch1-1.RData")
 ##  tr <- tr[tr$seal == "c026", ]
 ##  tr1 <- tr
 ## sub1 <- 17:24
 ## tr <- tr[sub1,]
 ## grd <- makeGridTopology(tr, cells.dim = c(5, 5), xlim = bbox(tr)[1,] + c(-11) * 0.1, ylim = bbox(tr)[2,] + c(-1,1) * 0.1)
 ## image(tripGrid(tr, grid = grd), col = topo.colors(256))
 ## grd <- makeGridTopology(tr, cells.dim = c(15, 15), xlim = bbox(tr)[1,] + c(1,1) * 0.1, ylim = bbox(tr)[2,] + c(-1,1) * 0.1)
 ## image(tripGrid(tr, grid = grd), col = topo.colors(256))

 tripGridEG <- function(duration, exact = FALSE, point.cex = 1, interp.cex = 0.5, cex = 1, main = "", show.text = TRUE, grid.col = "#aaaaaa") {

        if (missing(duration)) duration <- 0.5

        if (exact) {
            trg <- tripGrid(tr, grid = g)
        } else {

        trg <- tripGrid.interp(tr, dur = duration, grid = g)
        interp <- interpequal(tr, dur = duration)

    }

        trg[["z2"]] <- trg[["z"]] ^0.5

         image(as.image.SpatialGridDataFrame(trg["z2"]), col = hsv(0.6, c(0, seq(0.2, 0.8, length = 256)), 1), axes = FALSE)
        plot(pgs, add = TRUE, border = grid.col)

        title(main)
        points(tr, pch = 21, cex = point.cex, fg = "black", bg = "grey")

       if (!exact)        points(interp, pch = 4, interp.cex = 0.5)
        lines(coordinates(tr))
        if (show.text) text(as.data.frame(trg)[trg[["z"]] > 0, c("s1", "s2")], lab = format(trg[["z"]][trg[["z"]] > 0], digits = 3), cex = cex)
        invisible(NULL)

 }


#% par(mfcol = c(2,1))
#% g <- GridTopology(c(-0.1, -0.1), c(0.2, 0.2), c(6, 6))
#% pgs <- as.SpatialPolygons.GridTopology(g)
#% tripGridEG(exact = TRUE, main = "A", cex = 0.8)
