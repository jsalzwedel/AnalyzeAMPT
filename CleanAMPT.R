# Read in AMPT data and clean it up.


# Function for reading in data
GetData <- function(file, skip, nrows, nevent) {
    data <- read.table(file=file, skip=skip, nrows=nrows)
    data <- cbind(data, Event = nevent)
    colnames(data) <- c("PDG", "Px", "Py", "Pz", "Mass", "X", "Y", "Z", "T", "Event")
    data
}


filename <- "ampt.dat"
#nSkip <- c(1, 29452, 59407, 87393, 115940)
# How many rows does each event occupy
nRows <- c(29450, 29954, 27985, 28546, 29861, 28569, 28516, 28065, 25641, 29403,
           28932, 27734, 25276, 28533, 19728, 30372, 30497, 27872, 25639, 28208)
nEvent <- c(1:length(nRows))

# We need to skip the event information row between each event
nSkip <- 1
for(i in 1:(length(nRows)-1)) {
    nSkip <- c(nSkip, nSkip[i] + 1 + nRows[i])
}

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

# Write out the cleaned data
write.csv(ampt,"CleanData.csv")