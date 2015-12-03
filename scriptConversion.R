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
##sd02_deconvoluted[,"rt"] 
##sort(unlist(lapply(strsplit(as.character(sd02_deconvoluted[,4]), split = " _ "), "[", 2)))

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

mm <- binning(msp, tol = 0.01)

################################################################################
## similarity Matrix 
similarityMat <- similarityMatrix(mm)

## Clustering
library(amap)
hClustMSP <- hcluster(similarityMat, method = "spearman")
plot(hClustMSP, labels = FALSE)
colnames(similarityMat)[hClustMSP$order]
# hClustMSP2 <- hcluster(similarityMatrix[1:5,1:5], method = "spearman")
# plot(hClustMSP2, labels = FALSE)
# plot(hClustMSP2)

## order newSim according to order of clustering
colnames(similarityMat)
hClustMSP$order
newSim <- matrix(data=NA, nrow = dim(similarityMat)[1], ncol = dim(similarityMat)[1])
newColnames <- colnames(similarityMat)[hClustMSP$order]
colnames(newSim) <- rownames(newSim) <- newColnames

for (i in 1:dim(similarityMat)[1]) {
    newValues <- similarityMat[hClustMSP$order, hClustMSP$order[i]] 
    if (all(names(newSim[,i]) == names(newValues))) {
        newSim[,i] <- newValues
    } else {print("error");break}
}

write.csv(newSim, "orderedSimilarityMatrix.csv", sep=";", dec=".")
##
## circos

######################### 
## Neutral losses
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
## not finished
msp2FunctionalLossesMSP <- function(msp) {
    precmz <- getPrecursorMZ(msp)
    rt <- getRT(msp)
    indices <- getBegEndIndMSP(msp)
    indBegL <- indices[[1]]
    indEndL <- indices[[2]]
    ## create data frame for MSP file
    finalMSP <- matrix(data = NA, nrow = dim(msp)[1], ncol = 2) 
    finalMSP <- as.data.frame(finalMSP)
    
    ## create MSP from 
    for (i in 1:length(precmz)) {
        print(i)
        indBeg <- indBegL[i]
        indEnd <- indEndL[i]
        
        neutralLoss <- - (as.numeric(precmz[i]) - as.numeric(msp[indBeg:indEnd,1]))
        
        entry <- rbind(
            c("NAME: ", "Unknown"),
            c("RETENTIONTIME: ", rt[i]),
            c("PRECURSORMZ: ", precmz[i]),
            c("METABOLITENAME: ", "Unknown"),
            c("ADDUCTIONNAME: ", "Unknown"),
            c("Num Losses: ", length(indBeg:indEnd)),
            matrix(c(neutralLoss, msp[indBeg:indEnd,2]), ncol = 2),
            c(" ", " ")
        )
        entry <- as.matrix(entry)
        ## determine first empty line
        newstart <- which(is.na(finalMSP[,1]))[1]
        ## determine last line to write to
        newend <- newstart + dim(entry)[1] - 1
        finalMSP[newstart:newend,] <- entry
    }
    
    return(finalMSP)
}

nlMSP <- msp2FunctionalLossesMSP(finalMSP)

## create table with same fragments (binning)
mmNL <- binning(nlMSP, tol = 0.01, is = "Losses")

## similarity Matrix
similarityMat <- similarityMatrix(mmNL)

## Clustering
library(amap)
hClustMSP <- hcluster(similarityMat, method = "spearman")
plot(hClustMSP, labels = FALSE)
colnames(similarityMat)[hClustMSP$order]

## order newSim according to order of clustering
colnames(similarityMat)
hClustMSP$order
newSim <- matrix(data=NA, nrow = dim(similarityMat)[1], ncol = dim(similarityMat)[1])
newColnames <- colnames(similarityMat)[hClustMSP$order]
colnames(newSim) <- rownames(newSim) <- newColnames

for (i in 1:dim(similarityMat)[1]) {
    newValues <- similarityMat[hClustMSP$order, hClustMSP$order[i]] 
    if (all(names(newSim[,i]) == names(newValues))) {
        newSim[,i] <- newValues
    } else {print("error");break}
}

