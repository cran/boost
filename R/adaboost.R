adaboost <- function(xlearn, ylearn, xtest, presel = 200, mfinal = 100)
  {
    ## Feature Preselection
    if (presel > 0)
      {
        s       <- apply(xlearn, 2, score, ylearn)
        quality <- apply(rbind(s,-s+(sum(ylearn==0)*sum(ylearn==1))),2,max)
        genes   <- rev(order(quality))[1:presel]
        xlearn  <- xlearn[, genes]
        xtest   <- xtest[ , genes, drop = FALSE]
      }
    
    ## Initialization
    learn    <- dim(xlearn)[1]         
    test     <- dim(xtest)[1]
    Flearn   <- numeric(learn)
    Ftest    <- numeric(test)
    ptest    <- matrix(0, test, mfinal)
    w        <- rep(1/learn, learn)

    ## Shifting the labels
    ylearn[ylearn==0] <- -1

    ## Boosting Iterations
    for (m in 1:mfinal)
      {
        ## Fitting the tree
        update <- learner(ylearn, w, xlearn, xtest, method, args, bag = 1)
        flearn <- ((update$learn>0.5)*2)-1
        ftest  <- ((update$test> 0.5)*2)-1
         
        ## Updating and probabilities
        fehler   <- sum(w*(flearn!=ylearn))
        if (fehler>0)
          {
            consti   <- log((1-fehler)/fehler)
            w        <- w*exp(consti*(flearn!=ylearn))     
            w        <- pmax(w/sum(w), 1e-24)
            Flearn   <- Flearn + (consti*flearn)
            Ftest    <- Ftest  + (consti*ftest)
          }
        if (fehler==0)
          {
            w        <- w
            Flearn   <- Flearn + (1*flearn)
            Ftest    <- Ftest  + (1*ftest)
          }
        ptest[,m] <- 1/(1+exp((-2)*Ftest))
      }
   
    ## Output
    ptest
  }

