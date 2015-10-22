
## PNAS article
## sd01
sd01_outputXCMS <- read.csv("sd01_outputXCMS.csv", header=FALSE, sep=";")
colnames(sd01_outputXCMS) <- as.matrix(sd01_outputXCMS[3,])[1,]
sd01_outputXCMS <- sd01_outputXCMS[-c(1:3),]
sd01_FoldChangePvalue <- read.csv("sd01_FoldchangePvalue.csv", sep=";")
## sd02
sd02_similarityMatrix <- read.csv("sd02_similarityMatrix.csv", sep=";")
rownames(sd02_similarityMatrix) <- sd02_similarityMatrix[,1]
sd02_similarityMatrix <- sd02_similarityMatrix[,-1]
sd02_deconvoluted <- read.csv("sd02_deconvolutedidMSMS.csv", sep=";")


# Camera 
## msp file
# 1_ every one that has one identifier has one spectrum, take average as rt, leave adductionname empty

## identify precursor mz
precursor <- sd02_deconvoluted[,4]
precursor <- as.character(precursor)
splitPrecursor <- strsplit(precursor, split=" _ ")
splitPrecursor <- unlist(splitPrecursor)

lenSplitPrecursor <- length(splitPrecursor)
PrecursorMZ <- splitPrecursor[seq(2, lenSplitPrecursor, 2)]
lenPreMZ <- length(PrecursorMZ)

## add PrecursorMZ to deconvoluted idMSMS
sd02_deconvoluted <- cbind(sd02_deconvoluted, PrecursorMZ)

## change character to numeric
PrecursorMZ <- as.numeric(PrecursorMZ)
## PrecursorMZ <- sort(PrecursorMZ) ## ??
uniquePreMZ <- unique(precursor)
lenUniquePreMZ <- length(uniquePreMZ)
uniquePreMZ_cut <- unique(PrecursorMZ)

## create data frame for MSP file
finalMSP <- matrix(data = NA, nrow = 7 * lenUniquePreMZ + dim(sd02_deconvoluted)[1], ncol = 2) ## 7 new entries + all fragment ion entries
finalMSP <- as.data.frame(finalMSP)

## write to data frame
for (i in 1:lenUniquePreMZ) {
    ind <- which(uniquePreMZ[i] == precursor)    
    entry <- rbind(
        c("NAME: ", "Unknown"),
        c("RETENTIONTIME: ", mean(sd02_deconvoluted[ind,"rt"])),
        c("PRECURSORMZ: ", uniquePreMZ_cut[i]),
        c("METABOLITENAME: ", "Unknown"),
        c("ADDUCTIONNAME: ", "Unknown"),
        c("Num Peaks: ", length(ind)),
        sd02_deconvoluted[ind,c(1,3)],
        c(" ", " ")
    )
    entry <- as.matrix(entry)
    ## determine first empty line
    newstart <- which(is.na(finalMSP[,1]))[1]
    ## determine last line to write to
    newend <- newstart + dim(entry)[1] - 1
    finalMSP[newstart:newend,] <- entry
}

setwd("Documents/University/Master/HiWi/")
## write finalMSP to .msp
write.table(finalMSP, file = "idMSMStoMSP.msp", sep=" ", dec=".",
    row.names=FALSE, col.names=FALSE,quote=FALSE)


