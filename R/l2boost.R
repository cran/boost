l2boost <- function(xlearn, ylearn, xtest, presel = 200, mfinal = 100)
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
    Flearn   <- rep(0, learn)
    Ftest    <- rep(0, test)
    ptest    <- matrix(0, test, mfinal)

    ## Shifting the labels
    ylearn[ylearn == 0] <- -1

    ## Boosting Iterations
    for (m in 1:mfinal)
      {
        ## Fitting the tree
        ulearn <- ylearn-Flearn
        update <- learner(ulearn,rep(1,learn),xlearn,xtest,method,args,bag=1)
        flearn <- update$learn
        ftest  <- update$test
         
        ## Updating and probabilities
        Flearn    <- Flearn + flearn
        Ftest     <- Ftest  + ftest
        Flearn    <- pmin(Flearn, 1)
        Flearn    <- pmax(Flearn,-1)
        Ftest     <- pmin(Ftest,  1)
        Ftest     <- pmax(Ftest, -1)  
        ptest[,m] <- 1/(1+exp((-2)*Ftest))
      }
   
    ## Output
    ptest
  }

