
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

## precursor e.g. sd02_deconvoluted
## 1 - 183 (predefined)
allocatePrecursor2mz <- function(sd01, sd02, kNN = 10, mzCheck = 1, rtCheck = 30, mzVsRTbalance = 10000) {
    ## PARAMETERS
    ## k-nearest neighbours based on dev from m/z (i.e. the kNN entries with the 
    ## smallest deviation)
    ## mzCheck: maximum tolerated distance for mz (strong criterion here)
    ## rtCheck maximum tolerated distance for rt
    ## Multiplicator for mz value before calculating the (euclidean) distance between two peaks
    ## high value means that there is a strong weight on the dev m/z value
    ## mzVsRTbalance
    if (kNN < 0) break
    if (mzCheck < 0) break
    if (rtCheck < 0) break
    if (mzVsRTbalance < 0) break
    
    precursor <- sd02[,4]
    
    ## isolated mz values from e.g. pcgroup_precursorMZ column in sd02_deconvoluted
    uniquePrecursor <- .cutUniquePreMZ(precursor, splitPattern=" _ ", splitInd=2)

    ## create finalCluster, which is the data.frame to store data
    finalCluster <- matrix(nrow = length(uniquePrecursor), ncol = (5 + 183 + 4))
    colnames(finalCluster) <- c("Average Rt(min)", "Average mz", "Metabolite Name",
                                "Adduct ion name", "Spectrum reference file name", as.character(1:183), 
                                "check RT", "dev RT", "check mz", "deviation m/z")
    finalCluster <- as.data.frame(finalCluster)
    
    uniquePreMZ <- unique(precursor)
   
    
    ## LOOP WHICH WRITES TO finalCluster
    for (i in 1:length(uniquePrecursor)) {
        
        sd02mz <- uniquePrecursor[i]
        sd01mz <- sd01_outputXCMS[, "mz"]
        sd01mz <- as.character(sd01mz)
        sd01mz <- as.numeric(sd01mz)
        devmzOld <- devmz <- abs(sd02mz - sd01mz)
        
        ## use only kNN m/z dev
        sortDevMz <- sort(devmz)[1:kNN]
        
        ## get indices in sd01 of the smallest deviances
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
        
        ## calculate fake rt from sd02 (from fragment rt values)
        ind <- which(uniquePreMZ[i] == precursor)    
        sd02rt <- sd02[ind, "rt"]
        sd02rt <- mean(sd02rt) 
        
        ## determine rt values from sd01
        sd01rt <- sd01[indSortDevMZ, "rt"]
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
        XCMS <- sd01[minInd,]
        
        ## write new entry
        entry <- matrix(0, nrow = 1, ncol = ncol(finalCluster))
        colnames(entry) <- colnames(finalCluster)
        entry[, "Average Rt(min)"] <- sd02rt
        entry[, "Average mz"] <- uniquePrecursor[i]
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
}

