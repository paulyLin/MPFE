estimatePatterns <- function(patternCounts,
                             epsilon=0,
                             eta=0,
                             column=NULL,
                             fast=TRUE,
                             steps=20000)
{
    # Check the patterns
    # First, strip any leading character different from 01
    patterns <- gsub('^[^01]', '', as.character(patternCounts[, 1]))
    # Make sure we only have zeros and ones:
    matches <- grep('[01]+', patterns, invert=TRUE)
    if (length(matches) > 0) {
        stop('Only 0s and 1s allowed in pattern. First (up to ten) offending patterns:', head(patterns[matches], 10), '\n')
    }
    # Make sure that all patterns have the same length
    nCpGsites <- nchar(patterns[1])
    if (any(nchar(patterns) != nCpGsites)) {
        stop('All patterns must be the same length\n')
    }

    patternCounts[, 1] <- patterns
    # Check Epsilon
    if (!is.numeric(epsilon) || epsilon < 0 || epsilon >= 1) {
        stop('epsilon must be a numeric between 0 and 1\n')
    }

    # Check Eta
    nEtas <- length(eta)
    if (nEtas != 1 && nEtas != nCpGsites) {
        stop('Length of eta is not equal to the number of CpG sites\n')
    } else if (!is.numeric(eta) || any(eta < 0) || any(eta >= 1)) {
        stop('eta must be a numeric between 0 and 1\n')
    }

    # Check the column indices
    nColumns <- ncol(patternCounts) - 1
    if (is.null(column)) {
        columns <- 1:nColumns
    } else {
        columns <- column
    }
    if (any(columns < 1 || columns > nColumns)) {
        stop('Column indices must be between 1 and ', nColumns, '\n')
    }

    compareData <- list()
    for (i in 1:length(columns)) {
        compareData[[i]] <- estimatePatternsOneColumn(patternCounts,
                                                      epsilon,
                                                      eta,
                                                      column=columns[i],
                                                      fast,
                                                      steps)
    }

    return(compareData)
}

estimatePatternsOneColumn <- function(patternCounts,
                                      epsilon,
                                      eta,
                                      column,
                                      fast,
                                      steps)
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

    # Optimisation

    startingVector <- yWithoutMax / totalPatterns
    yZeros <- which(startingVector %in% c(0))
    startingVector[yZeros]<- yPatterns[patternsMax] / (100000 * totalPatterns)

    constraintMatrix <- rbind(diag(size - 1), rep(-1, size - 1))
    constraintVector <- append(rep(0, size - 1), -1)
    opt <- constrOptim(startingVector, likelihoodOpt, grad=NULL,
                        ui=constraintMatrix,
                        ci=constraintVector,
                        method='Nelder-Mead',
                        control=list(maxit=steps))

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

# Plot graphs.

plotPatterns <- function(compareData,
                         yLimit1=NULL,
                         yLimit2=NULL)
{
      if(class(compareData) != "data.frame") {
        stop("compareData must be a data.frame")
      }
      if(is.null(yLimit1)) {
        yLimit1 <- ceiling(max(compareData$observedDistribution, compareData$estimatedDistribution) * 10.2) / 10
      }
      if(is.null(yLimit2)) {
        yLimit2 <- quantile(sort(compareData$observedDistribution), 0.9)
      }

      layout(matrix(c(1, 2, 3), 3, 1, byrow=TRUE), heights=c(5, 1, 5))

      par(mar=c(5.4, 5, 2, 5))
      plot(compareData$estimatedDistribution,
           pch=4,
           col='white',
           xlab='pattern',
           ylab='proportion',
           cex.lab=1.5,
           ylim=c(0, yLimit1))
      points(compareData$estimatedDistribution, pch=4, col='red')
      par(new=TRUE)
      plot(compareData$coverage,
           col='blue',
           pch=3,
           xlab='',
           ylab='',
           ylim=c(0, sum(compareData$coverage) * yLimit1),
           axes=FALSE)
      axis(side=4)
      mtext('coverage', side=4, line=2.5)
      abline(0, 0, col='grey')

      par(mar=c(0, 0, 0, 5))
      plot.new()
      legend('right', c('observed distribution','estimated distribution'), pch=c(3,4), col=c('BLUE', 'RED'))

      par(mar=c(5.4, 5, 2, 5))
      plot(compareData$estimatedDistribution,
           pch=4,
           col='white',
           xlab='pattern',
           ylab='proportion',
           cex.lab=1.5,
           ylim=c(0, yLimit2))
      points(compareData$estimatedDistribution, pch=4, col='red')
      par(new=TRUE)
      plot(compareData$coverage,
           col='blue',
           pch=3,
           xlab='',
           ylab='',
           ylim=c(0,sum(compareData$coverage) * yLimit2),
           axes=FALSE)
      axis(side=4, at=c(compareData$coverage, 0, 10))
      mtext('coverage', side=4, line=2.5)
      abline(0, 0, col='grey')
}

# vim:ft=r:ts=4:sw=4:sts=4:expandtab:
