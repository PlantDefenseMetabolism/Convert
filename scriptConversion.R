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
    row.names=FALSE, col.names=TRUE,quote=FALSE)

################################################################################
## create table with same fragments 
finalMSP
msp <- finalMSP##[1:94,]
IndPrecMZ <- which(msp[,1] == "PRECURSORMZ: ")
precmz <- msp[IndPrecMZ,2]
IndRT <- which(msp[,1] == "RETENTIONTIME: ")
rt <- msp[IndRT,2]
indNumPeaks <- which(msp[,1] == "Num Peaks: ")  ## begining 

l <- list()

indEnd <- as.numeric(msp[indNumPeaks, 2])
inds <- indNumPeaks + 1 
indEnd <- indEnd + inds - 1 

for (i in 1:length(indEnd)) l[[i]] <- c(inds[i]:indEnd[i])
IndFrag <- unlist(l)
frag <- msp[IndFrag,1]
frag <- as.numeric(frag)

frag_s <- sort(frag)
frag_order <- order(frag)
## calculate distance to neighbours  
dist <- list(NA)
for (i in 1:length(frag_s)) {
    if (i != length(frag_s)) {
        dist[[i]] <- c(frag_s[i], frag_s[i+1] - frag_s[i])
    } else dist[[i]] <- c(frag_s[i], Inf)
}
    
distAdapt <- dist ## <- dist[1:30]
##dist[[length(dist)]][2] <- Inf

dist2Adapt <- dist2 <- lapply(dist, "[[", 2)
indConv <- which(unlist(dist2) == 0)

mapping <- list()
for (i in 1:length(dist)) mapping[[i]] <- i

for (i in 1:length(indConv))mapping[[indConv[i]]] <- indConv[i] + 1

inds <- which(unlist(dist2) == 0)
conv <- numeric(length(mapping)) ## conv is the vector which shows convoluted mz
x <- 1
for (i in 1:length(mapping)) {
    if (mapping[i] != i & conv[i] == "0") {
        conv[i] <- paste("M", x, sep="")
        conv[i+1] <- paste("M", x, sep="")
        j <- i
        while (unlist(mapping)[j+1] != (j +1) ) {
            conv[j + 1] <- paste("M", x, sep="")
            conv[j + 2] <- paste("M", x, sep="")
            j <- j +1 
        }
        i <- j
        x <- x +1 
    }
}

tol <- 0.01

indGreater0 <- which(unlist(dist2Adapt) > 0) 
minDist2AdaptGreater0 <- which.min(dist2Adapt[indGreater0])
indGreater0Min <- indGreater0[minDist2AdaptGreater0]

while (distAdapt[indGreater0Min][[1]][2] < tol) {
    ## write all which have Mx value to Mx+1
    if (conv[indGreater0Min + 1 ] == 0) {
        str <- unlist(strsplit(unique(conv),split="M"))
        str <- str[which(str != "")]
        if (0 %in% str) str <- str[which(str != "0")]
        str <- max(as.numeric(str)) + 1
        str <- paste("M", str, sep="")
        if (conv[indGreater0Min] == 0) {
            conv[indGreater0Min] <- str
            
        } else {
            conv[which(conv[indGreater0Min] == conv)] <- str
        }
        conv[indGreater0Min + 1] <- str   
    } else {
        if  (conv[indGreater0Min] == 0) {conv[indGreater0Min] <- conv[indGreater0Min + 1 ]}
        else {
            conv[which(conv[indGreater0Min] == conv)] <- conv[indGreater0Min + 1 ]}
    }
    ## calculate new mean for all instances with Mx+1
    indAdapt <- which(conv[indGreater0Min + 1 ] == conv)
    newMean <- mean(unlist(lapply(distAdapt[indAdapt], "[", 1)))
    ## write new mean
    for (i in indAdapt) distAdapt[[i]][1] <- newMean
    ## write new distances
    distAdaptOld <- distAdapt
    for (i in 1:length(distAdapt)) {
        if (i != length(distAdapt)) {
            distAdapt[[i]] <- c(distAdaptOld[[i]][1], distAdaptOld[[i+1]][1] - distAdaptOld[[i]][1])
        } else distAdapt[[i]] <- c(distAdaptOld[[i]][1], Inf)
    }
    dist2Adapt <- lapply(distAdapt, "[[", 2)
    unlist(dist2Adapt)
    
    indGreater0 <- which(unlist(dist2Adapt) > 0) 
    minDist2AdaptGreater0 <- which.min(dist2Adapt[indGreater0])
    indGreater0Min <- indGreater0[minDist2AdaptGreater0]
}

## write for every conv which has "0" a new Mx 
for (i in which(conv == "0")) {
    str <- unlist(strsplit(unique(conv),split="M"))
    str <- str[which(str != "")]
    if (0 %in% str) str <- str[which(str != "0")]
    str <- max(as.numeric(str)) + 1
    str <- paste("M", str, sep="")
    conv[i] <- str
}
    
uniqueMZ <- unlist(lapply(distAdapt, "[", 1))
uniqueMZ <- unique(uniqueMZ)

mm <- matrix(data = 0, nrow = length(precmz), ncol = length(uniqueMZ))
## convoluted MZ is column names
colnames(mm) <- uniqueMZ
rownames(mm) <- paste(precmz, rt, sep="/")



fragMM <- unlist(l)
## write to mm 
for (i in 1:length(fragMM)) {
    ## which frag is put first? second? ...
    indFragMM <- fragMM[frag_order][i]
    msp[indFragMM,]
    ## get corresponding Precursor MZ
    correspPrecMZ <- precmz[max(which(IndPrecMZ < indFragMM))]
    ## get corresponding rt
    correspPrecRT <- rt[max(which(IndPrecMZ < indFragMM))]
    ## get unique row identifier
    rowInd <- which(paste(correspPrecMZ, correspPrecRT, sep="/") == rownames(mm))
    ## get col index 
    colInd <- which(distAdapt[[i]][1] == colnames(mm))
    ## write 
    mm[rowInd,colInd] <- msp[indFragMM, 2]
}

## NDP, function to calculate normalised dot product
## NDP = sum(Ws1i * Ws2i)^2 / ( sum(Ws1i^2) * sum(Ws2i^2) ) 
## W = [peak intensity] ^ 0.5 [mass] ^ 2
NDP <- function(mat, row1=1, row2=2, m = 0.5, n = 2) {
    len <- dim(mat)[2]
    mass <- colnames(mat)
    mass <- as.numeric(mass)
    WS1 <- numeric(len)
    WS2 <- numeric(len)
    S1 <- as.numeric(mat[row1,])
    S2 <- as.numeric(mat[row2,])
    for (i in 1:len) {
        WS1[i] <- ( S1[i] ) ^ m * ( mass[i] ) ^ n
        WS2[i] <- ( S2[i] ) ^ m * ( mass[i] ) ^ n
    }
    NDP <- ( sum(WS1 * WS2) ) ^ 2 / (sum(WS1 ^ 2 ) * sum(WS2 ^ 2))
    return(NDP)
}

NDP(mm, row1=1, row2=2) 

similarityMatrix <- matrix(0, nrow = length(precmz), ncol = length(precmz))
colnames(similarityMatrix) <- rownames(similarityMatrix) <- paste(precmz, rt, sep="/")

similarityMatrix <- similarityMatrix
for (i in 1:dim(similarityMatrix)[1]) {
    for (j in 1:dim(similarityMatrix)[1]) {
        similarityMatrix[i,j] <- NDP(mm, row1=i, row2=j)
    }
}
##data.frame # rownames unique identifier, colnames = binned m/z
## entries 0 falls keine Ähnlichkeit, oder relative Intensität

