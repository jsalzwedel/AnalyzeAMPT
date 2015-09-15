# Read in AMPT data and analyze the relative separation of (anti)lambdas


# Function for reading in data
GetData <- function(file, skip, nrows, nevent) {
    data <- read.table(file=file, skip=skip, nrows=nrows)
    data <- cbind(data, Event = nevent)
    colnames(data) <- c("PDG", "Px", "Py", "Pz", "Mass", "X", "Y", "Z", "T", "Event")
    data
}

# How many rows to skip for each event, and how many to read in
filename <- "ampt.dat"
nSkip <- c(1, 29452, 59407, 87393, 115940)
nRows <- c(29450, 29954, 27985, 28546, 29861)
nEvent <- c(1:5)

# Read in the data
ampt <- GetData(filename, nSkip[1], nRows[1], nEvent[1])
for(i in 2:length(nEvent)) {
    tmp <- GetData(filename, nSkip[i], nRows[i], nEvent[i])
    ampt <- rbind(ampt,tmp)
}


# Convert to data.table for convenience
library(data.table)
ampt <- as.data.table(ampt)

AddColumns <- function(data) {
    # Add columns for magnitude of position and momentum, and for pseudorapidity
    data <- cbind(data, Rho = sqrt(data$X^2 + data$Y^2 + data$Z^2))
    data <- cbind(data, Pmag = sqrt(data$Px^2 + data$Py^2 + data$Pz^2))
    data <- cbind(data, Eta = atanh(data$Pz/data$Pmag))
    data
}

ampt <- AddColumns(ampt)

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
v0s <- rbind(lambdas,antilambdas)
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

GetDeltaRTable <- function(data) {
    xDiffs <- GetDiffsByEventForColumn(data,data$X)
    yDiffs <- GetDiffsByEventForColumn(data,data$Y)
    zDiffs <- GetDiffsByEventForColumn(data,data$Z)
    deltaR <- sqrt(xDiffs^2 + yDiffs^2 + zDiffs^2)
    dframe <- cbind(xDiffs,yDiffs,zDiffs, deltaR)
    as.data.table(dframe)
}

deltaRLambdas <- GetDeltaRTable(lambdas)
deltaRAntilambdas <- GetDeltaRTable(antilambdas)



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

GetDeltaRTableNotID <- function() {
    xDiffs <- GetDiffsByEventForColumnNotID(lambdas, lambdas$X, 
                                           antilambdas, antilambdas$X)
    yDiffs <- GetDiffsByEventForColumnNotID(lambdas, lambdas$Y, 
                                           antilambdas, antilambdas$Y)
    zDiffs <- GetDiffsByEventForColumnNotID(lambdas, lambdas$Z, 
                                           antilambdas, antilambdas$Z)
    deltaR <- sqrt(xDiffs^2 + yDiffs^2 + zDiffs^2)
    dframe <- cbind(xDiffs,yDiffs,zDiffs, deltaR)
    as.data.table(dframe)
}

deltaRNotID <- GetDeltaRTableNotID()

# Compute mean, standard deviation, and standard error
se <- function(x) sqrt(var(x)/length(x))

results <- data.frame(
    PDG = c(3122, -3122, 0),
    MeanDeltaR = c(mean(deltaRLambdas$deltaR), 
                   mean(deltaRAntilambdas$deltaR),
                   mean(deltaRNotID$deltaR)),
    Std = c(sd(deltaRLambdas$deltaR), 
            sd(deltaRAntilambdas$deltaR),
            sd(deltaRNotID$deltaR)),
    StdErr = c(se(deltaRLambdas$deltaR), 
               se(deltaRAntilambdas$deltaR),
               se(deltaRNotID$deltaR))              
)


