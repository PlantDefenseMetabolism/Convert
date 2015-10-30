## Data from Halle
spectra_0_less <- read.csv("20158251022_spectra_0_less.msp", 
                           comment.char="#")
rawmatrix_0_less <- read.delim("20158251022_rawmatrix_0_less.txt")

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



## 1) convert sd02_deconvoluted to .msp format

## identify precursor mz
precursor <- sd02_deconvoluted[,4]
precursor <- as.character(precursor)
splitPrecursor <- strsplit(precursor, split=" _ ") ## split precursors
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


## write finalMSP to .msp
write.table(finalMSP, file = "idMSMStoMSP.msp", sep=" ", dec=".",
    row.names=FALSE, col.names=FALSE,quote=FALSE)

################################################################################

## 2) convert sd02_deconvoluted precursor ions and sd01_outputXCMS to a format
## as in ..._rawmatrix_0_less.txt by allocating precursor ions to canditate
## m/z values based on minimal distance of m/z and deviance of rt based on 
## an objective function

## cf. group.nearest: Group peaks together across samples by creating a master 
## peak list and assigning corresponding peaks from all samples. 
## Arguments: mzVsRTbalnce Multiplicator for mz value before calculating the (euclidean) distance
## between two peaks (can we apply this here as we use different gradients?, cf. XCMS vignette: 
## "In most cases, LC/MS files that were acquired under different conditions should not
## be compared.  For instance, positive and negative ionization mode files will have no ions
## in common and should thus be preprocessed separately.  Similarly,  data files acquired
## with different elution gradients should not be processed together")
## mzCheck Maximum tolerated distance for mz
## rtCheck Maximum tolerated distance for RT
## kNN Number of nearest Neighbours to check

## isolated mz values from pcgroup_precursorMZ column in sd02_deconvoluted
mz <- uniquePreMZ_cut

## create finalCluster, which is the data.frame to store data
finalCluster <- matrix(nrow = length(mz), ncol = (5 + 183 + 4))
colnames(finalCluster) <- c("Average Rt(min)", "Average mz", "Metabolite Name",
    "Adduct ion name", "Spectrum reference file name", as.character(1:183), 
    "check RT", "dev RT", "check mz", "deviation m/z")
finalCluster <- as.data.frame(finalCluster)

## PARAMETERS
## k-nearest neighbours based on dev from m/z (i.e. the kNN entries with the 
## smallest deviation)
kNN <- 10 
mzCheck <- 1 ## maximum tolerated distance for mz (strong criterion here)
rtCheck <- 30 ## maximum tolerated distance for rt
## Multiplicator for mz value before calculating the (euclidean) distance between two peaks
## high value means that there is a strong weight on the dev m/z value
mzVsRTbalance <- 10000 

## LOOP WHICH WRITES TO finalCluster
for (i in 1:length(mz)) {
    
    sd02mz <- mz[i]
    sd01mz <- sd01_outputXCMS[, "mz"]
    sd01mz <- as.character(sd01mz)
    sd01mz <- as.numeric(sd01mz)
    devmzOld <- devmz <- abs(sd02mz - sd01mz)
    
    ## use only kNN m/z dev
    sortDevMz <- sort(devmz)[1:kNN]
    
    ## get indices in 01_outputXCMS of the smallest deviances
    indSortDevMZOld <- indSortDevMZ <- match(sortDevMz, devmz) 
    ## truncate devmz such that it only includes kNN m/z deviances
    devmz <- devmz[indSortDevMZ]
    
    ## check if devmz is in tolerated distance
    if (any(devmz <= mzCheck)) {
        ## truncate such that indSortDevMZ includes only indices and 
        ## devmz only m/z within the tolerance value
        indSortDevMZ <- indSortDevMZ[devmz <= mzCheck]
        devmz <- devmz[devmz <= mzCheck]
        ToleranceCheckMZ <- TRUE
    } else {
        print (c(i,"Deviation of m/z is greater than tolerance value. I won't truncate the kNN."))
        devmz <- devmz
        ToleranceCheckMZ <- FALSE
    }
    
    ## calculate fake rt from sd02_deconvoluted (from fragment rt values)
    ind <- which(uniquePreMZ[i] == precursor)    
    sd02rt <- sd02_deconvoluted[ind, "rt"]
    sd02rt <- mean(sd02rt) 
    
    ## determine rt values from sd01_outputXCMS
    sd01rt <- sd01_outputXCMS[indSortDevMZ, "rt"]
    sd01rt <- as.character(sd01rt)
    sd01rt <- as.numeric(sd01rt)
    
    devrt <- abs(sd02rt - sd01rt)
    
    if (any(devrt <= rtCheck)) {
        ## truncate devmz and devrt that it is included in the tolerance value
        devmz <- devmz[devrt <= rtCheck] 
        devrt <- devrt[devrt <= rtCheck]
        ToleranceCheckRT <- TRUE
        objective <- mzVsRTbalance * devmz + devrt
    } else {
        print(c(i, "Deviation of rt is greater than tolerance value. I won't use rt as a criterion."))
        ToleranceCheckRT <- FALSE
        objective <- devmz ## use only devmz
    }
    
    ## find smallest value for objective function 
    minInd <- which.min(objective) 
    ## get index in sd01_outputXCMS 
    minInd <- which(devmz[minInd] == devmzOld) 
    
    ## get the entry of sd01_outputXCMS with the smallest value
    XCMS <- sd01_outputXCMS[minInd,]
    
    ## write new entry
    entry <- matrix(0, nrow = 1, ncol = ncol(finalCluster))
    colnames(entry) <- colnames(finalCluster)
    entry[, "Average Rt(min)"] <- sd02rt
    entry[, "Average mz"] <- mz[i]
    entry[, "Metabolite Name"] <- "Unknown"
    entry[, "Adduct ion name"] <- if(nchar(as.character(XCMS[, "adduct"]) == 0)) "Unknown" else XCMS[,"adduct"]

    entry[, "Spectrum reference file name"] <- "Unknown"

    x <- XCMS[,which(colnames(XCMS) == "1"):which(colnames(XCMS) == "183")]
    x <- as.matrix(x)
    x <- as.vector(x)
    entry[, which(colnames(entry)=="1"):which(colnames(entry)=="183")] <- x

    entry[, "check RT"] <- ToleranceCheckRT
    entry[, "dev RT"] <- sd02rt - as.numeric(as.character(XCMS[,"rt"]))
    entry[, "check mz"] <- ToleranceCheckMZ
    entry[, "deviation m/z"] <- sd02mz - as.numeric(as.character(XCMS[,"mz"]))
    
    ## write to finalCluster
    finalCluster[i, ] <- entry
}

## write finalCluster 
write.table(finalCluster, file = "finalCluster.csv", sep=";", dec=".",
    row.names=FALSE, col.names=FALSE,quote=FALSE)
