data <- read.csv("DeltaR.csv", row.names=1)

library(data.table)
data <- as.data.table(data)

stripped <- data[,c("deltaR","pairType", "EvType"), with=FALSE]


GetDensityRatio <- function(data, type) {
    data <- data[data$pairType == type]
    sameD <- density(data$deltaR[data$EvType == "Same"], 
                     from= 0, to = 100, adjust = 2)
    mixD <- density(data$deltaR[data$EvType == "Mixed"],
                     from= 0, to = 100, adjust = 2)
    ratio <- sameD$y/mixD$y
    #plot(x=sameD$x, y=ratio) 
    ratioD <- mixD
    ratioD$y <- ratio
    
    list(Same = sameD, Mix = mixD, Ratio = ratioD)
}

laD <- GetDensityRatio(data, "LA")
llD <- GetDensityRatio(data, "LL")
aaD <- GetDensityRatio(data, "AA")

plot(laD$Ratio)

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
    hists
}

llHists <- MakeHistograms(data, "LL")
aaHists <- MakeHistograms(data, "AA")
laHists <- MakeHistograms(data, "LA")

# Then plot

MakeErrPlot <- function(hists) {
    DeltaR <- hists$Ratio$mids
    CF <- hists$Ratio$density
    errs <- hists$Errs
    errbar(DeltaR, CF, CF + errs, CF - errs  )
}
# errbar(xpoints, ypoint, yplus, yminus)
library(Hmisc)
MakeErrPlot(llHists)
MakeErrPlot(aaHists)
MakeErrPlot(laHists)
