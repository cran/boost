response1 <- function(zufall.x, zufall.y, gene = NULL, signs = NULL)
  {
    ## Searching discriminative genes
    if(is.null(gene) && is.null(signs))
      {
        nog           <- 10
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

    ## Determine probabilities, y-labels and bayes-error
    wert          <- rowSums(zufall.klein)
    wert          <- 5*((wert-mean(wert))/sd(wert))      
    probab        <- exp(wert)/(1+exp(wert))
    simu.y        <- numeric(nos)
    for(i in 1:nos)  simu.y[i] <- rbinom(1, 1, probab[i])
    bayes         <- sum((probab*(probab<=0.5))+((1-probab)*(probab>0.5)))/nos

    ## Output
    attributes(probab) <- NULL
    list(probab=probab, y=simu.y, bayes=bayes, gene=gene, signs=signs)
  }

