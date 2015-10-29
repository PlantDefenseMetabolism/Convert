## Data from Halle
spectra_0_less <- read.csv("20158251022_spectra_0_less.msp", 
                           comment.char="#")
rawmatrix_0_less <- read.delim("20158251022_rawmatrix_0_less.txt")

## PNAS article
## sd01
sd01_outputXCMS <- read.csv("~/Documents/University/Master/HiWi/sd01_outputXCMS.csv", header=FALSE, sep=";")
colnames(sd01_outputXCMS) <- as.matrix(sd01_outputXCMS[3,])[1,]
sd01_outputXCMS <- sd01_outputXCMS[-c(1:3),]
sd01_FoldChangePvalue <- read.csv("~/Documents/University/Master/HiWi/sd01_FoldchangePvalue.csv", sep=";")
## sd02
sd02_similarityMatrix <- read.csv("~/Documents/University/Master/HiWi/sd02_similarityMatrix.csv", sep=";")
rownames(sd02_similarityMatrix) <- sd02_similarityMatrix[,1]
sd02_similarityMatrix <- sd02_similarityMatrix[,-1]
sd02_deconvoluted <- read.csv("~/Documents/University/Master/HiWi/sd02_deconvolutedidMSMS.csv", sep=";")



## 1) convert sd02_deconvoluted to .msp format

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


## write finalMSP to .msp
write.table(finalMSP, file = "idMSMStoMSP.msp", sep=" ", dec=".",
    row.names=FALSE, col.names=FALSE,quote=FALSE)

################################################################################

## 2) convert sd02_deconvoluted precursor ions and sd01_outputXCMS to a format
## as in ..._rawmatrix_0_less.txt by allocating precursor ions to canditate
## m/z values based on minimal distance of m/z

mz <- uniquePreMZ_cut
indL <- list(NULL)
offset <- 0.5

minMZ <- as.character(sd01_outputXCMS[,3])
minMZ <- as.numeric(minMZ)
maxMZ <- as.character(sd01_outputXCMS[,4])
maxMZ <- as.numeric(maxMZ)

## preselection of possible entries
while (any(unlist(lapply(indL, length)) == 0)) {
    offset <- offset + 0.01
    for (i in 1:length(uniquePreMZ_cut)) {
        ind <- which( (minMZ-offset) <= mz[i] & mz[i] <= (maxMZ+offset) )
        indL[[i]] <- ind
    }
}

## find from these the preselection the entry with the smallest m/z deviation
## to the precursor mz

## create finalCluster, which is the data.frame to store data
finalCluster <- matrix(nrow = length(mz), ncol = (5 + 183 + 2))
colnames(finalCluster) <- c("Average Rt(min)", "Average mz", "Metabolite Name",
    "Adduct ion name", "Spectrum reference file name", as.character(1:183), 
    "dev RT", "minimum deviation m/z")
finalCluster <- as.data.frame(finalCluster)

for (i in 1:length(mz)) {
    
    ## preselection of m/z values
    mzXCMS <- sd01_outputXCMS[indL[[i]],"mz"]
    mzXCMS <- as.character(mzXCMS)
    mzXCMS <- as.numeric(mzXCMS)
    devMZ <- abs(mzXCMS - mz[i])
    ## find smallest deviation
    indMinmz <- which.min(devMZ)
    XCMS <- sd01_outputXCMS[indL[[i]],][indMinmz,]
    
    ## check if rt from deconvoluted is in range of rtmin to rtmax of sd01_outputXCMS
    ## and assign quality mark
    rtDeconvoluted <- sd02_deconvoluted[i,"rt"]
#     rtminXCMS <- XCMS[,"rtmin"]
#     rtminXCMS <- as.character(rtminXCMS)
#     rtminXCMS <- as.numeric(rtminXCMS)
# 
#     rtmaxXCMS <- XCMS[,"rtmax"]
#     rtmaxXCMS <- as.character(rtmaxXCMS)
#     rtmaxXCMS <- as.numeric(rtmaxXCMS)

    rt_sd01_outputXCMS <- XCMS[, "rt"]
    rt_sd01_outputXCMS <- as.character(rt_sd01_outputXCMS)
    rt_sd01_outputXCMS <- as.numeric(rt_sd01_outputXCMS)
    
    devRT <- rtDeconvoluted - rt_sd01_outputXCMS
    
    ##qualityRT <- rtminXCMS <= rtDeconvoluted & rtDeconvoluted <= rtmaxXCMS
    ##if (qualityRT) {quality <- 1} else {quality <- 0}
    
    indRT <- which(uniquePreMZ[i] == precursor)   
    meanFakeRT <- mean(sd02_deconvoluted[indRT,"rt"])
    
    entry <- matrix(0, nrow = 1, ncol = ncol(finalCluster))
    colnames(entry) <- colnames(finalCluster)
    entry[, "Average Rt(min)"] <- meanFakeRT
    entry[, "Average mz"] <- mz[i]
    entry[, "Metabolite Name"] <- "Unknown"
    entry[, "Adduct ion name"] <- if(nchar(as.character(XCMS[, "adduct"]) == 0)) "Unknown" else XCMS[,"adduct"]
    
    entry[, "Spectrum reference file name"] <- "Unknown"
    
    x <- XCMS[,which(colnames(XCMS) == "1"):which(colnames(XCMS) == "183")]
    x <- as.matrix(x)
    x <- as.vector(x)
    entry[, which(colnames(entry)=="1"):which(colnames(entry)=="183")] <- x
            
    entry[, "dev RT"] <- devRT
    entry[, "minimum deviation m/z"] <- min(devMZ)
    
    finalCluster[i, ] <- entry
}

## 
plot(finalCluster[, "dev RT"], finalCluster[, "minimum deviation m/z"], 
     xlab = "deviance in RT", ylab = "deviation in m/z")
abline(h = 2)
abline(h = 1.0, lty = 2)






##
# 
# rtXCMS <- sd01_outputXCMS[indL[[1]], "rt"]
# rtXCMS <- as.character(rtXCMS)
# rtXCMS <- as.numeric(rtXCMS)
# 
# rtPre <- sd02_deconvoluted[1,"rt"]
# minRTind <- which.min(abs(rtXCMS - rtPre))
# 
# ind <- which(sd01_outputXCMS[,"rt"] == rtXCMS[minRTind])
# sd01_outputXCMS[ind,]
# 
# 
# ## average rt(min), average mz; trio, lvs = replicates, take retention time from average rt from S2, intensities from S1xcmscamera
# ## precursor mz is slightly different from xcmscamera: how to threshold
# 
# range(sd01_outputXCMS[,2], sd01_outputXCMS[,3])
# 
# 
# precursor <- sd02_deconvoluted[,4]
# precursor <- as.character(precursor)
# splitPrecursor <- strsplit(precursor, split=" _ ")
# splitPrecursor <- unlist(splitPrecursor)
# 
# lenSplitPrecursor <- length(splitPrecursor)
# PrecursorMZ <- splitPrecursor[seq(2, lenSplitPrecursor, 2)]
# lenPreMZ <- length(PrecursorMZ)
# 
# ## add PrecursorMZ to deconvoluted idMSMS
# sd02_deconvoluted <- cbind(sd02_deconvoluted, PrecursorMZ)
# 
# ## change character to numeric
# PrecursorMZ <- as.numeric(PrecursorMZ)
# ## PrecursorMZ <- sort(PrecursorMZ) ## ??
# uniquePreMZ <- unique(precursor)
# lenUniquePreMZ <- length(uniquePreMZ)
# uniquePreMZ_cut <- unique(PrecursorMZ)