response3 <- function(zufall.x, zufall.y, gene = NULL, signs = NULL)
  {
    ## Searching discriminative genes
    if(is.null(gene) && is.null(signs))
      {
        nog           <- 25
        scores        <- apply(zufall.x,2,score,zufall.y)
        signs         <- scores<((table(zufall.y)[1]*table(zufall.y)[2])/2)
        signs         <- (signs-0.5)*2
        scores        <- abs(scores-(table(zufall.y)[1]*table(zufall.y)[2])/2)
        gene          <- rev(order(scores))[1:nog]
        signs         <- signs[gene]
      }

    ## Matrix of (flipped) genes for the response model
    nos                   <- dim(zufall.x)[1]
    zufall.klein          <- zufall.x[,gene]
    change                <- which(signs[gene]<0)
    zufall.klein[,change] <- -zufall.klein[,change]
    
    ## Determine probabilities, y-labels and bayes error
    betas         <- runif(25, 0, 2)
    gammas        <- runif(25, 0, 0.2)
    deltas        <- runif(25, 0, 0.1)
    term3         <- (zufall.klein%*%deltas)
    term2         <- (zufall.klein%*%gammas)
    term1         <- (zufall.klein%*%betas)
    term1         <- term1-mean(term1)
    term2         <- term2-mean(term2)
    term3         <- term3-mean(term3)
    wert          <- term1*(1+term2)*(1+term3)
    wert          <- 10*((wert-mean(wert))/sd(wert))
    wert          <- pmin(200, wert)
    probab        <- exp(wert)/(1+exp(wert))
    simu.y        <- numeric(nos)
    for(i in 1:nos)    simu.y[i] <- rbinom(1, 1, probab[i])
    bayes         <- sum((probab*(probab<=0.5))+((1-probab)*(probab>0.5)))/nos

    ## Output
    attributes(probab) <- NULL
    list(probab=probab, y=simu.y, bayes=bayes, gene=gene, signs=signs)
  }

