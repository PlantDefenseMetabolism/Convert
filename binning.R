## function to get precursor mz of a .msp
getPrecursorMZ <- function (msp) {
    ## get indices with precursor mz
    IndPrecMZ <- which(msp[,1] == "PRECURSORMZ: ")
    ## get precursor mz
    precmz <- msp[IndPrecMZ,2]
    
    return(precmz)
}
    
## function to get rt of a .msp
getRT <- function (msp) {
    ## get indices with rt
    IndRT <- which(msp[,1] == "RETENTIONTIME: ")
    ## get rt 
    rt <- msp[IndRT,2]
    
    return(rt)
}

## function to get beginning and end indices of fragments of a .msp
getBegEndIndMSP <- function(msp, is = c("Fragments", "Losses")) {
    
    is <- match.arg(is)
    ## beginning 
    if (is == "Fragments") indNumPeaks <- which(msp[,1] == "Num Peaks: ") 
    if (is == "Losses") indNumPeaks <- which(msp[,1] == "Num Losses: ")
    
    indEnd <- as.numeric(msp[indNumPeaks, 2])
    indBeg <- indNumPeaks + 1 
    indEnd <- indEnd + indBeg - 1 
    
    return(list(indBeg, indEnd))
}

## function to bin mz values 
binning <- function(msp, tol = 0.01, is = c("Fragments", "Losses")) { 
    
    is <- match.arg(is)
    
    ## msp is .msp file
    ## tol is tolerance value for binning
    precmz <- getPrecursorMZ(msp)
    rt <- getRT(msp)
    
    indices <- getBegEndIndMSP(msp, is = is)
    indBeg <- indices[[1]]
    indEnd <- indices[[2]]
    
    ## create list with indices in msp for each precursor
    l <- lapply(1:length(indBeg), function(x) c(indBeg[x]:indEnd[x]))
    
    IndFrag <- unlist(l)
    ## fragments 
    frag <- msp[IndFrag, 1] 
    frag <- as.numeric(frag)
    ## sort frag and get order
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
    
    ## get distance values
    dist2Adapt <- dist2 <- lapply(dist, "[[", 2)
    indConv <- which(unlist(dist2) == 0) ## get indices which have distance of 0
    
    mapping <- lapply(1:length(dist), function(x) x)

    ## map to the next element if it has same mz
    for (i in 1:length(indConv)) mapping[[indConv[i]]] <- indConv[i] + 1
    
    ## conv is the vector which shows convoluted mz
    conv <- numeric(length(mapping)) 
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
    
    IndPrecMZ <- match(precmz, msp[,2])
    
    for (i in 1:length(fragMM)) {
        ## which frag is put first? second? ...
        indFragMM <- fragMM[frag_order][i]
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

    return(mm)
}




