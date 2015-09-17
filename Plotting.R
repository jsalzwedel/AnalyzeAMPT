# Let's plot the data in various ways

data <- read.csv("DeltaR.csv", row.names=1)

library(data.table)
data <- as.data.table(data)

#stripped <- data[,c("deltaR","pairType", "EvType"), with=FALSE]


GetDensityRatio <- function(data, type) {
    data <- data[data$pairType == type]
    sameD <- density(data$deltaR[data$EvType == "Same"], 
                     from= 0, to = 100, adjust = 2)
    mixD <- density(data$deltaR[data$EvType == "Mixed"],
                     from= 0, to = 100, adjust = 2)
    ratio <- sameD$y/mixD$y
    ratioD <- mixD
    ratioD$y <- ratio
    
    list(Same = sameD, Mix = mixD, Ratio = ratioD)
}

laD <- GetDensityRatio(data, "LA")
llD <- GetDensityRatio(data, "LL")
aaD <- GetDensityRatio(data, "AA")

plot(aaD$Ratio)

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


MakeHistograms <- function(data, type) {
    theData <- data[data$pairType == type]
    histSame<- with(theData, hist(deltaR[EvType == "Same" & deltaR < 100], 
                  xlim = range(c(0,100)), 
                  breaks = seq(0,100, by=5)))
    histMixed<- with(theData, hist(deltaR[EvType == "Mixed" & deltaR < 100], 
                   xlim = range(c(0,100)), 
                   breaks = seq(0,100, by=5)))
    histRatio<- histSame
    histRatio$counts <- histSame$counts/histMixed$counts
    histRatio$density <- histSame$density/histMixed$density
    hists <- list(Same = histSame, Mix = histMixed, Ratio = histRatio)
    errs <- GetErrorBars(hists)
    hists$Errs <- errs
    hists$Type <- type
    hists
}

llHists <- MakeHistograms(data, "LL")
aaHists <- MakeHistograms(data, "AA")
laHists <- MakeHistograms(data, "LA")

MakeHistDataTable <- function(data) {
    dt <- data.table(DeltaR = data$Ratio$mids,
                CF = data$Ratio$density,
                Errs = data$Errs,
                Type = data$Type)
}

histDT <- rbind(MakeHistDataTable(llHists),
                MakeHistDataTable(aaHists),
                MakeHistDataTable(laHists))

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

MakeErrPlot(llHists)
MakeErrPlot(aaHists)
MakeErrPlot(laHists)

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

panels <- MakeErrPlotFromDT(histDT)


