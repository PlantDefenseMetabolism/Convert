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

similarityMatrix <- function(mm) {
    n <- dim(mm)[1]
    similarity <- matrix(0, nrow = n, ncol = n)
    colnames(similarity) <- rownames(similarity) <- rownames(mm)
       ## paste(precmz, rt, sep="/")    
    
    ## write to similarity matrix similarit measure
    for (i in 1:n) {
        for (j in 1:n) {
            if (i <= j) {
                similarity[j,i] <- similarity[i,j] <- NDP(mm, row1=i, row2=j)
            }
        }
    }
    return(similarity)
    
}
    


n <- dim(similarityMatrix)[1]
## write to similarity matrix similarit measure
for (i in 1:n) {
    for (j in 1:n) {
        if (i <= j) {
            similarityMatrix[j,i] <- similarityMatrix[i,j] <- NDP(mm, row1=i, row2=j)
        }
    }
}

