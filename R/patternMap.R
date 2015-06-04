patternMap <- function(patterns,
                       minFreq=0,
                       maxFreq=1,
                       noSpurious=TRUE,
                       estimatedDistribution=TRUE,
                       topDown=TRUE,
                       allTicks=FALSE,
                       methCol='black',
                       unMethCol='white',
                       ...)
{
    # If we use the estimated distribution, the spurious patterns are at zero
    if (estimatedDistribution)
    {
        noSpurious <- TRUE
    }
    # Remove the spurious patterns if necessary
    if (noSpurious)
    {
        patterns <- patterns[patterns$spurious == FALSE, ]
    }
    distOfInterest <- ifelse(estimatedDistribution,
                             'estimatedDistribution',
                             'observedDistribution')
    # Trim low frequencies
    patterns <- patterns[patterns[, distOfInterest] >= minFreq, ]
    # Trim high frequencies
    patterns <- patterns[patterns[, distOfInterest] <= maxFreq, ]
    # Number of patterns left
    nPatterns <- nrow(patterns)
    # Invert the list if necessary
    patterns <- patterns[order(patterns[, distOfInterest]), ]
    if (!topDown)
    {
        patterns <- patterns[nPatterns:1, ]
    }
    # Size of the patterns
    patternSize <- max(nchar(as.character(patterns$pattern)))
    # Compute the frequencies
    y0 <- sum(patterns[, distOfInterest])
    ys <- y0 - cumsum(patterns[, distOfInterest])

    # Colours

    # Auxilary function to create the colour vector
    getColours <- function(col)
    {
        if (is.function(col))
        {
            cols <- col(nPatterns)
        }
        else
        {
            return(rep(col, length=nPatterns))
        }
    }
    methCol <- getColours(methCol)
    unMethCol <- getColours(unMethCol)

    # Plot the main area
    plot(0,
         type="n",
         xlab="CpG",
         ylab="Proportion",
         xlim=c(0.5, patternSize + 0.5),
         #xaxp=c(1, patternSize, patternSize - 1), 
         ylim=c(0, y0),
         ...)

    # Add ticks for each CpG if required
    if (allTicks)
    {
        axis(1, at=1:patternSize, labels=FALSE)
    }

    # Draw the patterns
    for (i in 1:nPatterns)
    {
        patternString <- as.character(patterns$pattern[i])
        pattern       <- as.integer(strsplit(patternString, '')[[1]])
        freq          <- patterns[i, distOfInterest]
        x0            <- 0.5
        y1            <- ys[i]
        for (c in pattern)
        {
            colour <- ifelse(c == 1, methCol[i], unMethCol[i])
            x1     <- x0 + 1
            rect(x0, y1, x1, y0, col=colour, border='NA')
            x0     <- x1
        }
        y0 <- y1
    }
}

# vim:ft=r:ts=4:sw=4:sts=4:expandtab:
