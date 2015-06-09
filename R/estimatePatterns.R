estimatePatterns <- function(patternCounts,
                             epsilon=0,
                             eta=0,
                             column=NULL,
                             fast=TRUE,
                             steps=20000,
                             reltol=1e-12)
{
    # Check the patterns
    # First, strip any leading character different from 01
    patterns <- gsub('[^01]', '', as.character(patternCounts[, 1]))
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
                                                      steps,
                                                      reltol)
    }

    return(compareData)
}
