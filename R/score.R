score <- function(x,y)
  {
    wilcox.test(x[which(y==0)], x[which(y==1)])$statistic
  }
