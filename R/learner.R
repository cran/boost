learner <- function(y, w, xlearn, xtest, method, args, bag)
  {
    ## Definitions
    learn  <- dim(xlearn)[1]
    test   <- dim(xtest)[1]
    blearn <- matrix(0, bag, learn)
    btest  <- matrix(0, bag, test)
    
    ## Currently only stumps as learners are supported, no choice of args!!!
    cntrl <- rpart.control(maxdepth = 1, minsplit = learn-1, #minbucket = 1,
                           maxsurrogate = 0, usesurrogate=0, maxcompete = 1,
                           cp = 0, xval = 0)
    
    ## Bagging stumps/trees
    if (bag==1)
      {
        bx         <- xlearn
        fit        <- rpart(y~bx, weights = w/mean(w), control = cntrl)
        bx         <- xtest
        blearn[1,] <- predict(fit)
        btest[1,]  <- predict(fit, newdata = data.frame(bx))
      }
    if (bag>1)
      {
        for (b in 1:bag)
          {
            indices    <- sample(1:learn, learn, replace = TRUE)
            by         <- y[indices]
            bw         <- w[indices]
            bx         <- xlearn[indices,]
            fit        <- rpart(by~bx, weights=bw/mean(bw), control=cntrl)
            bx         <- xlearn
            blearn[b,] <- predict(fit, newdata = data.frame(bx))
            bx         <- xtest
            btest[b,]  <- predict(fit, newdata = data.frame(bx))
          }
      }

    ## Output
    list(learn = apply(blearn, 2, mean), test = apply(btest, 2, mean))
  }

