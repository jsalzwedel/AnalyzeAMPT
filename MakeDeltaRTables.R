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
    deltaR <- sqrt(xDiffs^2 + yDiffs^2 + zDiffs^2)
    dframe <- cbind(xDiffs,yDiffs,zDiffs, deltaR, pairType, EvType = "Same")
    as.data.table(dframe)
}

llSame <- GetDeltaRTable(lambdas, "LL")
aaSame <- GetDeltaRTable(antilambdas, "AA")



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

GetDeltaRTableNotID <- function(pairType) {
    xDiffs <- GetDiffsByEventForColumnNotID(lambdas, lambdas$X, 
                                           antilambdas, antilambdas$X)
    yDiffs <- GetDiffsByEventForColumnNotID(lambdas, lambdas$Y, 
                                           antilambdas, antilambdas$Y)
    zDiffs <- GetDiffsByEventForColumnNotID(lambdas, lambdas$Z, 
                                           antilambdas, antilambdas$Z)
    deltaR <- sqrt(xDiffs^2 + yDiffs^2 + zDiffs^2)
    dframe <- cbind(xDiffs,yDiffs,zDiffs, deltaR, pairType, EvType = "Same")
    as.data.table(dframe)
}

laSame <- GetDeltaRTableNotID("LA")





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
    xDiffs <- GetDiffsMixedEvent(data1, data1$Y, data2, data2$Y)
    yDiffs <- GetDiffsMixedEvent(data1, data1$Y, data2, data2$Y)
    zDiffs <- GetDiffsMixedEvent(data1, data1$Y, data2, data2$Y)
    deltaR <- sqrt(xDiffs^2 + yDiffs^2 + zDiffs^2)
    dframe <- cbind(xDiffs,yDiffs,zDiffs, deltaR, pairType, EvType = "Mixed")
    as.data.table(dframe)
}

llMix <- DoMixing(lambdas,lambdas, "LL")
aaMix <- DoMixing(antilambdas,antilambdas, "AA")
laMix <- DoMixing(lambdas, antilambdas, "LA")

combined <- rbind(llSame, aaSame, laSame, llMix, aaMix, laMix)
write.csv(combined, "DeltaR.csv")

# Compute mean, standard deviation, and standard error
se <- function(x) sqrt(var(x)/length(x))

MakeResultsTable <- function(ll, aa, la) {
    results <- data.table(
        PairType = c(ll$pairType, aa$pairType, la$pairType),
        MeanDeltaR = c(mean(ll$deltaR), 
                       mean(aa$deltaR),
                       mean(la$deltaR)),
        Std = c(sd(ll$deltaR), 
                sd(aa$deltaR),
                sd(la$deltaR)),
        StdErr = c(se(ll$deltaR), 
                   se(aa$deltaR),
                   se(la$deltaR))              
    )   
}

resultsSame <- MakeResultsTable(llSame, aaSame, laSame)
resultsMixed <- MakeResultsTable(llMix, aaMix, laMix)


# Exploratory plotting
lambdasXYZ <- lambdas[,intersect(colnames(lambdas), c("X","Y","Z")), with=FALSE]
pairs(lambdasXYZ)
antiXYZ <- antilambdas[,intersect(colnames(lambdas), c("X","Y","Z")), with=FALSE]
pairs(antiXYZ)



# Let's make density functions (histograms) for the deltaR distributions
sameD <- density(laSame$deltaR, from= 0, to = 100, n = 64)
mixD <- density(laMix$deltaR, from= 0, to = 100, n = 64)

plot(sameD)
plot(mixD)
ratio <- sameD$y/mixD$y
plot(x=sameD$x, y=ratio) 
# It would be good to get error bars on this
