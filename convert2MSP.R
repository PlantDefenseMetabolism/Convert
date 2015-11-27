## helper function to get unique precursor ions
.cutUniquePreMZ <- function(precursor, splitPattern = splitPattern, splitInd = splitInd) {
    ## split precursors according to split pattern
    precursor <- as.character(precursor)
    splitPrecursor <- strsplit(precursor, split = splitPattern)
    ## extract precursor mz at position splitInd
    splitPrecursor <- lapply(splitPrecursor,"[", splitInd)
    PrecursorMZ <- unlist(splitPrecursor)
    lenPreMZ <- length(PrecursorMZ)
    
    ## change character to numeric
    PrecursorMZ <- as.numeric(PrecursorMZ)
    uniquePreMZ <- unique(precursor)
    lenUniquePreMZ <- length(uniquePreMZ)
    uniquePreMZ_cut <- unique(PrecursorMZ)
    
    return(uniquePreMZ_cut)
}


## function to convert mm into MSP format

## mm needs to be in the format: 
## 1st column: mz, 2nd column: rt, 3rd column: intensity, 
## 4th column: pcgroup_precursorMZ, (5th column: PrecursorMZ)
## concerning 4th column, split is the pattern which separates elements and precursor mz
## splitInd is the position of the precursor mz concerning separatation by split pattern

convert2MSP <- function (mm, splitPattern = "_", splitInd = 1) {
    precursor <- mm[,4]
    precursor <- as.character(precursor)
    
    uniquePreMZ <- unique(precursor)
    uniquePreMZ_cut <- .cutUniquePreMZ(precursor = precursor, split = splitPattern, splitInd = splitInd)
    lenUniquePreMZ <- length(uniquePreMZ_cut)
    
    ## add PrecursorMZ to deconvoluted idMSMS
    ## mm <- cbind(mm, PrecursorMZ)
    
    ## create data frame for MSP file
    finalMSP <- matrix(data = NA, nrow = 7 * lenUniquePreMZ + dim(mm)[1], ncol = 2) ## 7 new entries + all fragment ion entries
    finalMSP <- as.data.frame(finalMSP)
    
    ## write to data frame
    for (i in 1:lenUniquePreMZ) {
        ind <- which(uniquePreMZ[i] == precursor)    
        entry <- rbind(
            c("NAME: ", "Unknown"),
            c("RETENTIONTIME: ", mean(mm[ind,"rt"])),
            c("PRECURSORMZ: ", uniquePreMZ_cut[i]),
            c("METABOLITENAME: ", "Unknown"),
            c("ADDUCTIONNAME: ", "Unknown"),
            c("Num Peaks: ", length(ind)),
            mm[ind,c(1,3)],
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
