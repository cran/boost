summarize <- function(boost.out, resp, mout=ncol(boost.out), grafik=TRUE)
  {
    mcra <- apply(((boost.out>0.5)*1)!=resp, 2, mean)
    mini <- which.min(mcra)
    mcrs <- round(min(mcra), 4)
    mcrf <- round(mean(((boost.out[,mout]>0.5)*1)!=resp),4)
    cat("\n")
    cat("Minimal mcr:  ",mcrs,"achieved after",mini,"boosting step(s)\n")
    cat("Fixed mcr:    ",mcrf,"achieved after",mout,"boosting step(s)\n")
    if (grafik)
      {
        xax <- "Boosting steps"
        yax <- "Error rate"
        ttl <- "LogitBoost"
        plot(mcra, xlab=xax, ylab=yax, main=ttl, type="l")
      }
  }
    
