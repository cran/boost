simulator <- function(x, y, respmod = c("none", "resp1", "resp2", "resp3"),
                      nos = 1200, gene = NULL, signs = NULL)
  {
    ## Argument matching
    respmod <- match.arg(respmod)

    ## Determine the correlation structure and the means
    x.0               <- x[which(y==0),]
    sigma.0           <- var(x.0)
    zerleg.0          <- eigen(sigma.0, symmetric = TRUE)
    ew.0              <- zerleg.0$values
    ew.0[ew.0<10^-10] <- 0
    ev.0              <- zerleg.0$vectors
    wurzel.0          <- ev.0%*%diag(sqrt(ew.0))
    mean.0            <- apply(x.0,2,mean)*1
    x.1               <- x[which(y==1),]
    sigma.1           <- var(x.1)
    zerleg.1          <- eigen(sigma.1, symmetric = TRUE)
    ew.1              <- zerleg.1$values
    ew.1[ew.1<10^-10] <- 0
    ev.1              <- zerleg.1$vectors
    wurzel.1          <- ev.1%*%diag(sqrt(ew.1))
    mean.1            <- apply(x.1,2,mean)*1

    ## Simulation of gene expression profiles
    nvars         <- ncol(x)
    n1            <- round(nos/2)
    n2            <- nos-n1
    u             <- matrix(rnorm(nvars*n1), nvars)
    mittel        <- matrix(rep(mean.0,n1),n1,nvars,byrow=TRUE)
    zufall.0      <- t(wurzel.0%*%u) + mittel
    u             <- matrix(rnorm(nvars*n2), nvars)
    mittel        <- matrix(rep(mean.1,n2),n2,nvars,byrow=TRUE)
    zufall.1      <- t(wurzel.1%*%u) + mittel
    zufall.x      <- rbind(zufall.0, zufall.1)
    zufall.y      <- c(rep(0,n1),rep(1,n2))

    ## Calling the response model
    response <- switch(respmod,
                       none  = list(y=zufall.y),
                       resp1 = response1(zufall.x, zufall.y, gene, signs),
                       resp2 = response2(zufall.x, zufall.y, gene, signs),
                       resp3 = response3(zufall.x, zufall.y, gene, signs))

    ## That's it!
    dimnames(zufall.x) <- NULL
    c(list(x=zufall.x), response)
  }

