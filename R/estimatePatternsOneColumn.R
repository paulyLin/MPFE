estimatePatternsOneColumn <- function(patternCounts,
                                      epsilon,
                                      eta,
                                      column,
                                      fast,
                                      steps,
                                      reltol)
{
    patternCounts <- patternCounts[, c(1, column + 1)]
    names(patternCounts) <- c('patterns','counts')   
    patternCounts <- patternCounts[patternCounts$counts !=0, ]
    nCpGsites <- nchar(patternCounts$patterns[1])

    # Create the vector of patterns

    if (fast) {
        methData <- patternCounts
    } else {
            binary <- function(x) if (all(x < 2)) x else cbind(binary(x %/% 2), x %% 2)
            cytosineBinary <- binary(0:(2 ^ nCpGsites - 1))
            mPattern <- array(dim=c(2^nCpGsites, 1))
            for(i in 1:(2 ^ nCpGsites)) {
                mPattern[i, 1] <- paste(cytosineBinary[i, ], collapse='')
            }
            counts <- array(0, dim=c(2^nCpGsites,1))
            methData <- data.frame(patterns=mPattern, counts=counts)
            methData$patterns <- as.character(methData$patterns)
            for(i in 1:(2 ^ nCpGsites)) {
                patternMatches <- patternCounts$patterns == mPattern[i, 1]
                if(any(patternMatches)) {
                    methData[i, 2] <- patternCounts[patternMatches, 2]
                }
            }
    }

    yPatterns <- methData$counts
    totalPatterns <- sum(yPatterns)
    readDistribution <- yPatterns / totalPatterns
    patternsMax <- which.max(yPatterns)
    yWithoutMax <- yPatterns[-patternsMax]

    # Create the pattern array

    size <- nrow(methData)

     stringToVector <- function(x) {
         patternVector <- array(dim=nchar(x))
         for (i in 1:nchar(x)) {
             patternVector[i] <- as.numeric(substr(x, i, i))
         }
         return(patternVector)
     }

     patternArray <- array(dim=c(size, nCpGsites))

     for (i in 1:size) {
      patternArray[i, ] <- stringToVector(methData$patterns[i])
     }

     # Construct the conversion matrix

    if (fast) {
        conversionMatrix <- array(0, dim=c(size, size))
        for (i in 1:size) {
            for (j in 1:size) {
                fromIndex <- patternArray[i, ]
                toIndex <- patternArray[j, ]
                conversionMatrix[i, j] <- prod((1 - epsilon - eta + 2 * epsilon * eta) ^ ((fromIndex + toIndex == 0) * 1) *
                                               (epsilon + eta - 2 * epsilon * eta) ^ ((fromIndex - toIndex == -1) * 1) *
                                               eta ^ ((fromIndex - toIndex == 1) * 1) *
                                               (1 - eta) ^ ((fromIndex + toIndex == 2) * 1))
             }
        }
    } else {
        conversionRule <- list(array(0, dim=c(2, 2)))
        for (i in 1:nCpGsites) {
          conversionRule[[i]] <- array(c(1 - epsilon - eta[i] + 2 * epsilon * eta[i],
                                         eta[i],
                                         epsilon + eta[i]- 2 * epsilon * eta[i],
                                         1 - eta[i]),
                                       dim=c(2, 2))
        }
        conversionMatrix <- 1
        for (i in 1:nCpGsites) {
            conversionMatrix <- kronecker(conversionMatrix, conversionRule[[i]])
        }
    }


    # Define the likelihood function
    # Note that likelihoodOpt's argument has one less entry than the function likelihood

    likelihood <- function(theta) {
        phi <- theta %*% conversionMatrix
        if (fast) {
            lhd <- -sum(yPatterns * log(as.vector(phi)))
        } else {
            lhd <- -sum(yPatterns[yPatterns != 0] * log(phi[yPatterns != 0]))
        }
        return(lhd)
    }

    expand <- function(theta, patternsMax) {
        expanded <- c(theta[(1:length(theta)) < patternsMax],
                      1 - sum(theta),
                      theta[(1:length(theta)) >= patternsMax])
        return(expanded)
    }

    likelihoodOpt <- function(theta) {
        likelihoodOpted <- likelihood(expand(theta, patternsMax))
        return(likelihoodOpted)
    }
	
	likelihood_grad <- function(theta){
		theta2 <- expand(theta, patternsMax)
		phi <- as.vector(theta2 %*% conversionMatrix)
		
		if (fast) {
			likelihood_grad2 <- colSums(apply(conversionMatrix, 1, function(x) yPatterns*x/phi))
		} else {
			likelihood_grad2 <- colSums(apply(conversionMatrix, 1, 
											function(x) (yPatterns*x)[yPatterns !=0]/phi[yPatterns !=0]))
		}
		likelihood_grad <- -likelihood_grad2[-patternsMax]+sum(yPatterns*conversionMatrix[patternsMax,]/phi)
		return(likelihood_grad)
	}
	
		
    # Optimisation

    startingVector <- yWithoutMax / totalPatterns
    yZeros <- which(startingVector %in% c(0))
    startingVector[yZeros]<- yPatterns[patternsMax] / (100000 * totalPatterns)

    constraintMatrix <- rbind(diag(size - 1), rep(-1, size - 1))
    constraintVector <- append(rep(0, size - 1), -1)
    opt <- constrOptim(startingVector, likelihoodOpt, grad=likelihood_grad,
                        ui=constraintMatrix,
                        ci=constraintVector,
                        method='BFGS',
                        control=list(maxit=steps, reltol=reltol))

    recovered <- expand(opt$par, patternsMax)

    # Decide the "0's"

    minZeros <- recovered
    minLikelihood <- likelihood(recovered)
    listReduced <- 1:size

    for (i in listReduced[listReduced != patternsMax]) {
        minTemp <- minZeros
        minTemp[patternsMax] <- minTemp[patternsMax] + minTemp[i]
        minTemp[i] <- 0
        if(likelihood(minTemp) <= minLikelihood) {
            minZeros <- minTemp
            minLikelihood <- likelihood(minTemp)
        }
    }

    # Generate output

    compareData <- data.frame(pattern=methData$patterns,
                              coverage=yPatterns,
                              observedDistribution=readDistribution,
                              estimatedDistribution=minZeros)

    if(!fast) {
        compareData <- compareData[compareData$observedDistribution != 0 | compareData$estimatedDistribution != 0, ]
        rownames(compareData) <- 1:nrow(compareData)
    }

    if(opt$convergence!=0) {
        warning('Constrained optimisation did not converge.\n')
    }

    compareData$spurious <- compareData$observedDistribution !=0 & compareData$estimatedDistribution == 0

    return(compareData)
}

