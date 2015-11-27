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
finalMSP <- convert2MSP(sd02_deconvoluted, split = " _ ", splitInd = 2)

## write finalMSP to .msp
write.table(finalMSP, file = "idMSMStoMSP.msp", sep=" ", dec=".",
            row.names=FALSE, col.names=FALSE,quote=FALSE)


############ tissue project
## 1) convert idMSMStissueproject.txt to .msp format
tissue <- read.csv("idMSMStissueproject.txt",sep="\t")
## change tissue to format that it is compatible with writeMSP sd02_deconvoluted
## 1st column: mz, 2nd column: rt, 3rd column: intensity, 
## 4th column: pcgroup_precursorMZ, (5th column: PrecursorMZ)
newtissue <- tissue
newtissue[,2] <- tissue[,3]; colnames(newtissue)[2] <- colnames(tissue)[3]
newtissue[,3] <- tissue[,2]; colnames(newtissue)[3] <- colnames(tissue)[2]
newtissue[,4] <- tissue[,7]; colnames(newtissue)[4] <- colnames(tissue)[7]
newtissue <- newtissue[,1:4]
tissue <- newtissue

finalMSP <- convert2MSP(tissue, split = "_", splitInd=1)


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
uniquePreMZ_cut <- .cutUniquePreMZ(sd02_deconvoluted[,4], splitPattern=" _ ", splitInd=2)

finalCluster <- allocatePrecursor2mz(sd01 = sd01_outputXCMS, sd02 = sd02_deconvoluted, rtCheck = 60) 

## write finalCluster 
write.table(finalCluster, file = "finalCluster.csv", sep=";", dec=".",
            row.names=FALSE, col.names=TRUE,quote=FALSE)

################################################################################
## 3) create table with same fragments (binning)
msp <- finalMSP
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
    
distAdapt <- dist

dist2Adapt <- dist2 <- lapply(dist, "[[", 2)
indConv <- which(unlist(dist2) == 0)

mapping <- list()
for (i in 1:length(dist)) mapping[[i]] <- i

## map to the next element if it has same mz
for (i in 1:length(indConv))mapping[[indConv[i]]] <- indConv[i] + 1

inds <- which(unlist(dist2) == 0)
conv <- numeric(length(mapping)) ## conv is the vector which shows convoluted mz
x <- 1

## for mz values which have distance of 0 create a identifier Mx, where x is an
## increasing number to be able to trace back same mz
for (i in 1:length(mapping)) {
    if (mapping[i] != i & conv[i] == "0") {
        conv[i] <- paste("M", x, sep="")
        conv[i+1] <- paste("M", x, sep="")
        j <- i
        ## check when there is a sequence of mz which have distance 0 and 
        ## allocate then the identical mx
        while (unlist(mapping)[j+1] != (j +1) ) { 
            conv[j + 1] <- paste("M", x, sep="")
            conv[j + 2] <- paste("M", x, sep="")
            j <- j +1 
        }
        i <- j
        x <- x +1 
    }
}

## tolerance value for binning
tol <- 0.01

## actual binning script starts here 
indGreater0 <- which(unlist(dist2Adapt) > 0) ## get distances greater 0
## find smallest distance which is greater than zero
minDist2AdaptGreater0 <- which.min(dist2Adapt[indGreater0])
indGreater0Min <- indGreater0[minDist2AdaptGreater0]

while (distAdapt[indGreater0Min][[1]][2] < tol) {
    ## write all which have Mx value to Mx+1
    if (conv[indGreater0Min + 1 ] == 0) { ## then create new Mx
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
    } else { ## if conv[indGreater0min + 1] != 0, i.e. if it is Mx, then use "old" Mx
        if (conv[indGreater0Min] == 0) {conv[indGreater0Min] <- conv[indGreater0Min + 1 ]}
        else {
            conv[which(conv[indGreater0Min] == conv)] <- conv[indGreater0Min + 1 ]}
    }
    ## calculate new mean for all instances with Mx+1
    indAdapt <- which(conv[indGreater0Min + 1 ] == conv)
    newMean <- mean(unlist(lapply(distAdapt[indAdapt], "[", 1)))
    ## write new mean to all instances with Mx+1
    for (i in indAdapt) distAdapt[[i]][1] <- newMean
    ## calculate new distances and write new distances
    distAdaptOld <- distAdapt
    for (i in 1:length(distAdapt)) { ## calculate for all elements in the list (this can be 
        ## changed, so that we only calculate distance for elements before and after Mx+1)
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

## actual binning script ends here

## write for every conv which has "0" a new Mx 
for (i in which(conv == "0")) {
    str <- unlist(strsplit(unique(conv),split="M"))
    str <- str[which(str != "")]
    if (0 %in% str) str <- str[which(str != "0")] 
    str <- max(as.numeric(str)) + 1 ## find highest x and create new one (+1)
    str <- paste("M", str, sep="")
    conv[i] <- str ## allocate Mx+1 to conv[i]
}

## find all unique bins, these will be the colnames of mm
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
    ##WS1 <- numeric(len)
    ##WS2 <- numeric(len)
    S1 <- as.numeric(mat[row1,])
    S2 <- as.numeric(mat[row2,])
    ##for (i in 1:len) { ## formula according to Li et al 2015, PNAS
    ##    WS1[i] <- ( S1[i] ) ^ m * ( mass[i] ) ^ n ## WS1 <- S1 ^ m * mass ^ n
    ##    WS2[i] <- ( S2[i] ) ^ m * ( mass[i] ) ^ n ## WS2 <- S2 ^ m * mass ^ n
    ##}
    WS1 <- S1 ^ m * mass ^ n
    WS2 <- S2 ^ m * mass ^ n
    NDP <- ( sum(WS1 * WS2) ) ^ 2 / (sum(WS1 ^ 2 ) * sum(WS2 ^ 2))
    return(NDP)
}


## create similarity matrix which contains pairwise similarity measure NDP
similarityMatrix <- matrix(0, nrow = length(precmz), ncol = length(precmz))
colnames(similarityMatrix) <- rownames(similarityMatrix) <- paste(precmz, rt, sep="/")

n <- dim(similarityMatrix)[1]
## write to similarity matrix similarit measure
for (i in 1:n) {
    for (j in 1:n) {
        if (i <= j) {
            similarityMatrix[j,i] <- similarityMatrix[i,j] <- NDP(mm, row1=i, row2=j)
        }
    }
}

## Clustering
library(amap)
hClustMSP <- hcluster(similarityMatrix, method = "spearman")
plot(hClustMSP, labels = FALSE)
colnames(similarityMatrix)[hClustMSP$order]
# hClustMSP2 <- hcluster(similarityMatrix[1:5,1:5], method = "spearman")
# plot(hClustMSP2, labels = FALSE)
# plot(hClustMSP2)

## order newSim according to order of clustering

colnames(similarityMatrix)
hClustMSP$order
newSim <- matrix(data=NA, nrow = dim(similarityMatrix)[1], ncol = dim(similarityMatrix)[1])
newColnames <- colnames(similarityMatrix)[hClustMSP$order]
colnames(newSim) <- rownames(newSim) <- newColnames

for (i in 1:dim(similarityMatrix)[1]) {
    newValues <- similarityMatrix[hClustMSP$order, hClustMSP$order[i]] 
    if (all(names(newSim[,i]) == names(newValues))) {
        newSim[,i] <- newValues
    } else {print("error");break}
}

write.csv(newSim, "orderedSimilarityMatrix.csv", sep=";", dec=".")
##
## circos
NL <- c("CH2", "CH4", "NH3", "H2O", "K+toNH4+", "Na+toH+", "C2H2", 
        "CO", "C2H4", "CH3N", "CH2O", "CH5N", "S", "H2S", 
        "K+toH+", "C2H2O", "C3H6", "CHNO", "CO2", "CH2O2", "C4H8",
        "C3H9N", "C2H4O2", "CH4N2O", "SO2", "C5H8", "C3H6O2", "C6H6",
        "SO3", "C3H2O3", "C4H8O2", "C4H12N2", "H2SO4", "H3PO4", "C5H10O2",
        "C3H4O4", "C6H12O2", "C2H5O4P", "C5H8O4", "C7H19N3", "C6H10O4", "C6H10O5",
        "C6H12O5", "C6H8O6", "C6H12O6", "C6H10O7", "C8H12O6", "C11H10O4", "C10H15N3O6S", 
        "C10H17N3O6S", "C12H20O10", "C12H22O11")

mzNL <- c(14.0157, 16.0313, 17.0265, 18.0106, 20.9293, 21.9819, 26.0157, 
        27.9949, 28.0313, 29.0266, 30.0106, 31.0422, 31.9721, 33.9877,
        37.9559, 42.0106, 42.0470, 43.0058, 43.9898, 46.0055, 56.0626,
        59.0735, 60.0211, 60.0324, 63.9619, 68.0626, 74.0368, 78.0470,
        79.9568, 86.0004, 88.0517, 88.1000, 97.9674, 97.9769, 102.0618,
        104.0110, 116.0861, 123.9926, 132.0423, 145.1579, 146.0579, 162.0528,
        164.0685, 176.0321, 180.0634, 194.0427, 204.0655, 206.0579, 305.0682,
        307.0838, 324.1057, 342.1162)
nl <- matrix(mzNL, nrow=1, ncol = length(NL))
colnames(nl) <- NL

x <- mm[1,]
frag <- names(x)
frag <- as.numeric(frag)
precursor <- rownames(mm)[1]
precursor <- strsplit(precursor, "/")[[1]]
precursor <- as.numeric(precursor[1])
precursorInd <- which.min(abs(frag-precursor))
PRECURSOR <- x[precursorInd]
PRECURSORiMZ <- as.numeric(names(PRECURSOR))
PRECURSORiIntensity <- as.numeric(PRECURSOR)
indsFrag <- which(as.numeric(x) != 0)

FRAGj <- x[indsFrag[2]]
FRAGjMZ <- as.numeric(names(FRAGj))
FRAGjIntensity <- as.numeric(FRAGj)
FRAGjMZ - PRECURSORiMZ
