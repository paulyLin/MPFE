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
