# Let's plot the data in various ways
library(data.table)
data <- as.data.table(read.csv("DeltaR.csv", row.names=1))

datapp <- as.data.table(read.csv("DeltaRProts.csv", row.names=1))

# GetDensityRatio <- function(data, type) {
#     data <- data[data$pairType == type]
#     sameD <- density(data$zDiffs[data$EvType == "Same"], 
#                      from= 0, to = 100, adjust = 2)
#     mixD <- density(data$zDiffs[data$EvType == "Mixed"],
#                      from= 0, to = 100, adjust = 2)
#     ratio <- sameD$y/mixD$y
#     ratioD <- mixD
#     ratioD$y <- ratio
#     
#     list(Same = sameD, Mix = mixD, Ratio = ratioD)
# }
# 
# laD <- GetDensityRatio(data, "LA")
# llD <- GetDensityRatio(data, "LL")
# aaD <- GetDensityRatio(data, "AA")
# 
# plot(aaD$Ratio)

# Now write a function that calculates error bars

GetErrorBars <- function(hists) {
    #seNum <- sqrt(hists$Same$counts)
    #seDen <- sqrt(hists$Mixed$counts)
    numInverseCounts <- 1/hists$Same$counts
    denInverseCounts <- 1/hists$Mix$counts
    numInverseCounts[numInverseCounts == Inf]  <- 0
    fractionalErr <- sqrt(numInverseCounts + denInverseCounts)
    errBar <- fractionalErr * hists$Ratio$density
}


############
# Try passing a column name (deltaR, deltaX, whatever) directly
# Then use data[,column, with = FALSE]
###########


MakeHistograms <- function(data, type, colname, rmin, rmax) {
    data<- data[data$pairType == type]
    column <- unlist(data[,colname, with=FALSE], use.names=FALSE)
    histSame<- hist(column[data$EvType == "Same" &
                               column < rmax & 
                               column > rmin], 
                  xlim = range(c(rmin,rmax)), 
                  breaks = seq(rmin,rmax, by=5))
    histMixed<- hist(column[data$EvType == "Mixed" & 
                                column < rmax &
                                column > rmin], 
                   xlim = range(c(rmin,rmax)), 
                   breaks = seq(rmin,rmax, by=5))

    histRatio<- histSame
    histRatio$counts <- histSame$counts/histMixed$counts
    histRatio$density <- histSame$density/histMixed$density
    hists <- list(Same = histSame, Mix = histMixed, Ratio = histRatio)
    errs <- GetErrorBars(hists)
    hists$Errs <- errs
    hists$Type <- type
    hists
}



MakeHistDataTable <- function(data) {
    dt <- data.table(DeltaR = data$Ratio$mids,
                CF = data$Ratio$density,
                Errs = data$Errs,
                Type = data$Type)
}

MakeAllDT <- function(data, var, rmin, rmax, isLL=TRUE) {   
    types <- character()
    if(isLL) types <- c("LL","AA","LA")
    else types <- c("PP","ApAp","PAp")
    llHists <- MakeHistograms(data, types[1], var, rmin, rmax)
    aaHists <- MakeHistograms(data, types[2], var, rmin, rmax)
    laHists <- MakeHistograms(data, types[3], var, rmin, rmax)
    
    histDT <- rbind(MakeHistDataTable(llHists),
                    MakeHistDataTable(aaHists),
                    MakeHistDataTable(laHists))    
}


deltaRHists <- MakeAllDT(data,"deltaR", 0, 100)
deltaRtHists <- MakeAllDT(data,"deltaRt", 0, 100)
deltaZHists <- MakeAllDT(data,"zDiffs", -50, 50)
deltaTHists <- MakeAllDT(data,"tDiffs", -50, 50)

deltaRHistsP <- MakeAllDT(datapp,"deltaR", 0, 100, FALSE)
deltaRtHistsP <- MakeAllDT(datapp,"deltaRt", 0, 100, FALSE)
deltaZHistsP <- MakeAllDT(datapp,"zDiffs", -50, 50, FALSE)
deltaTHistsP <- MakeAllDT(datapp,"tDiffs", -50, 50, FALSE)
# Then plot

library(Hmisc) # Used for errbar function
library(ggplot2)
MakeErrPlot <- function(hists) {
    x <- hists$Ratio$mids
    y <- hists$Ratio$density
    errs <- hists$Errs
    #errbar(DeltaR, CF, CF + errs, CF - errs  )
    
    dat <- data.table(x,y,errs)
    
    limits <- aes(ymax = y+errs, ymin = y-errs)
    ggplot(dat, aes(x=x, y=y)) +
        geom_errorbar(limits) +
        geom_point()
}

#MakeErrPlot(llHists)
#MakeErrPlot(aaHists)
#MakeErrPlot(laHists)

# Try plotting these side by side

MakeErrPlotFromDT <- function(histDT) {
#    limits <- aes(ymax = y+errs, ymin = y-errs)
    ggplot(histDT, aes(x=DeltaR, y=CF, color = Type)) +
        geom_errorbar(aes(ymin = CF-Errs, ymax = CF+Errs)) +
        geom_point() +
        #geom_line() +
        facet_wrap(~Type, nrow=1) +
        xlab(expression(paste(Delta,"R"))) +
        ylab(expression(paste("CF(",Delta,"R)"))) +
        coord_cartesian(ylim = c(0,2))
}

panelsR <- MakeErrPlotFromDT(deltaRHists)
panelsR

panelsRt <- MakeErrPlotFromDT(deltaRtHists)
panelsRt

panelsZ <- MakeErrPlotFromDT(deltaZHists)
panelsZ

panelsT <- MakeErrPlotFromDT(deltaTHists)
panelsT

panelsRP <- MakeErrPlotFromDT(deltaRHistsP)
panelsRP

panelsRtP <- MakeErrPlotFromDT(deltaRtHistsP)
panelsRtP

panelsZP <- MakeErrPlotFromDT(deltaZHistsP)
panelsZP

panelsTP <- MakeErrPlotFromDT(deltaTHistsP)
panelsTP

