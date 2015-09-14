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
nEvent <- c(1, 2, 3, 4, 5)

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

