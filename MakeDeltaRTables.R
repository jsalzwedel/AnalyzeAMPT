ampt <- read.csv("CleanData.csv", row.names=1)

library(data.table)
ampt <- as.data.table(ampt)
# Find particles that pass PDG, eta, and rho cuts
GetParticles <- function(data, thisPDG, etaMin = -0.8, etaMax = 0.8, rho = 10^2) {
    tmp <- data[data$PDG == thisPDG]
    tmp <- tmp[tmp$Eta > etaMin]
    tmp <- tmp[tmp$Eta < etaMax]
    tmp <- tmp[tmp$Rho < rho]
    tmp
}

# Find midrapidity lambdas (throw away obvious secondaries)
lambdas <- GetParticles(ampt, thisPDG = 3122)
antilambdas <- GetParticles(ampt, thisPDG = -3122)

protons <- GetParticles(ampt, thisPDG = 2212)
antiprotons <- GetParticles(ampt, thisPDG = -2212)


# Now lets make Delta R distributions
# For some subset of values, take the difference between each value
DiffIdentical <- function(vals) {
    diffMatrix <- sapply(vals, function(x) x - vals)
    diffs <- diffMatrix[upper.tri(diffMatrix)]
}

# Split the lambda X distribution up by event, then take the differences in X
GetDiffsByEventForColumn <- function(data, column) {
    splitData <- with(data, split(column, Event))
    splitDataDiffs <- sapply(splitData, DiffIdentical)
    diffsCombined <- unlist(splitDataDiffs, use.names=FALSE)
}

GetDeltaRTable <- function(data, pairType) {
    xDiffs <- GetDiffsByEventForColumn(data,data$X)
    yDiffs <- GetDiffsByEventForColumn(data,data$Y)
    zDiffs <- GetDiffsByEventForColumn(data,data$Z)
    tDiffs <- GetDiffsByEventForColumn(data,data$T)
    deltaRt <- sqrt(xDiffs^2 + yDiffs^2)
    deltaR <- sqrt(xDiffs^2 + yDiffs^2 + zDiffs^2)
    dt <- data.table(xDiffs, yDiffs, zDiffs, tDiffs,
                     deltaRt, deltaR,
                     pairType, EvType = "Same")
}


# Try doing this for non-identical particles.  Here, we need to mix between species
DiffNonIdentical <- function(vals1, vals2) {
    diffMatrix <- sapply(vals1, function(x) x - vals2)
}

# This works, but would be nice to just pass a column name in the future
GetDiffsByEventForColumnNotID <- function(data, column, data2, data2column) {
    splitData1 <- with(data, split(column, Event))
    splitData2 <- with(data2, split(data2column, Event))
    splitDataDiffs <- mapply(DiffNonIdentical, splitData1, splitData2)
    diffsCombined <- unlist(splitDataDiffs, use.names=FALSE)
}

GetDeltaRTableNotID <- function(data1, data2, pairType) {
    xDiffs <- GetDiffsByEventForColumnNotID(data1, data1$X, 
                                           data2, data2$X)
    yDiffs <- GetDiffsByEventForColumnNotID(data1, data1$Y, 
                                           data2, data2$Y)
    zDiffs <- GetDiffsByEventForColumnNotID(data1, data1$Z, 
                                           data2, data2$Z)
    tDiffs <- GetDiffsByEventForColumnNotID(data1, data1$T, 
                                            data2, data2$T)
    deltaRt <- sqrt(xDiffs^2 + yDiffs^2)
    deltaR <- sqrt(xDiffs^2 + yDiffs^2 + zDiffs^2)
    dt <- data.table(xDiffs, yDiffs, zDiffs, tDiffs,
                     deltaRt, deltaR,
                     pairType, EvType = "Same")
}



llSame <- GetDeltaRTable(lambdas, "LL")
aaSame <- GetDeltaRTable(antilambdas, "AA")
laSame <- GetDeltaRTableNotID(lambdas, antilambdas, "LA")

ppSame <- GetDeltaRTable(protons, "PP")
apapSame <- GetDeltaRTable(antiprotons, "ApAp")
papSame <- GetDeltaRTableNotID(protons, antiprotons, "PAp")



# Event mixing diffs
GetDiffsMixedEvent <- function(data, column, data2, data2column) {
    splitData1 <- with(data, split(column, Event))
    splitData2 <- with(data2, split(data2column, Event))
    diffsCombined <- numeric()
    for(i in 1:(length(splitData1)-1)) {
        for (j in (i+1):length(splitData2)) {
            splitDataDiffs <- DiffNonIdentical(splitData1[[i]],splitData2[[j]])
            diffVec <- unlist(splitDataDiffs, use.names = FALSE)
            diffsCombined <- c(diffsCombined, diffVec)
        }
    }
    diffsCombined
}


DoMixing <- function(data1, data2, pairType) {
    xDiffs <- GetDiffsMixedEvent(data1, data1$X, data2, data2$X)
    yDiffs <- GetDiffsMixedEvent(data1, data1$Y, data2, data2$Y)
    zDiffs <- GetDiffsMixedEvent(data1, data1$Z, data2, data2$Z)
    tDiffs <- GetDiffsMixedEvent(data1, data1$T, data2, data2$T)
    deltaRt <- sqrt(xDiffs^2 + yDiffs^2)
    deltaR <- sqrt(xDiffs^2 + yDiffs^2 + zDiffs^2)
    dframe <- data.table(xDiffs, yDiffs, zDiffs, tDiffs,
                         deltaRt, deltaR,
                         pairType, EvType = "Mixed")
}

llMix <- DoMixing(lambdas,lambdas, "LL")
aaMix <- DoMixing(antilambdas,antilambdas, "AA")
laMix <- DoMixing(lambdas, antilambdas, "LA")

ppMix <- DoMixing(protons,protons, "PP")
apapMix <- DoMixing(antiprotons, antiprotons, "ApAp")
papMix  <- DoMixing(protons,antiprotons, "PAp")

combined <- rbind(llSame, aaSame, laSame, llMix, aaMix, laMix)
write.csv(combined, "DeltaR.csv")
rm(llSame, aaSame, laSame, llMix, aaMix, laMix)

combinedpp <- rbind(ppSame, apapSame, papSame, ppMix, apapMix, papMix)
write.csv(combinedpp, "DeltaRProts.csv")
rm(ppSame, apapSame, papSame, ppMix, apapMix, papMix)

# Compute mean, standard deviation, and standard error
se <- function(x) sqrt(var(x)/length(x))

MakeResultsTable <- function(data) {
    results <- data.table(MeanDeltaR = mean(data$deltaR),
                          StdErrR = se(data$deltaR),
                          StdR = sd(data$deltaR),
                          MeanDeltaRt = mean(data$deltaRt),
                          StdErrRt = se(data$deltaRt),
                          StdRt = sd(data$deltaRt),
                          PairType = data$pairType[1],
                          EvType = data$EvType[1]
                          )
                          
}

combinedSummary <- rbind(MakeResultsTable(llSame),
                         MakeResultsTable(aaSame),
                         MakeResultsTable(laSame),
                         MakeResultsTable(llMix),
                         MakeResultsTable(aaMix),
                         MakeResultsTable(laMix))


# Exploratory plotting
lambdasXYZ <- lambdas[,intersect(colnames(lambdas), c("X","Y","Z")), with=FALSE]
pairs(lambdasXYZ)
antiXYZ <- antilambdas[,intersect(colnames(lambdas), c("X","Y","Z")), with=FALSE]
pairs(antiXYZ)

