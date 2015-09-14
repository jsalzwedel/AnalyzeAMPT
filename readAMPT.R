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
nSkip <- c(1, 29452, 59407, 87393)
nRows <- c(29450, 29954, 27985, 28546)
nEvent <- c(1, 2, 3, 4)

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

# Find lambdas

#lambdas <- ampt[ampt$PDG == 3122]
#antilambdas <- ampt[ampt$PDG == -3122]

# Throw away the obvious secondaries
