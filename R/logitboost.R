logitboost <- function(xlearn, ylearn, xtest, presel = 200, mfinal = 100)
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

    ## Length of training and test data
    learn    <- dim(xlearn)[1]         
    test     <- dim(xtest)[1]
    
    ## Initialization
    Flearn   <- numeric(learn)             
    Ftest    <- numeric(test)        
    plearn   <- rep(1/2, learn)
    ptest    <- matrix(0, test, mfinal)
    
    ## Boosting Iterations
    for (m in 1:mfinal)
      {
        ## Computation of working response and weights
        w    <- pmax(plearn*(1-plearn), 1e-24)
        z    <- (ylearn-plearn)/w               

        ## Fitting the weak learner
        update <- learner(z, w, xlearn, xtest, method, args, bag = 1)
        flearn <- update$learn
        ftest  <- update$test
         
        ## Updating and probabilities
        Flearn    <- Flearn + (1/2)*flearn
        Ftest     <- Ftest + (1/2)*ftest
        plearn    <- 1/(1+exp((-2)*Flearn))
        ptest[,m] <- 1/(1+exp((-2)*Ftest))
      }
   
    ## Output
    ptest
  }

