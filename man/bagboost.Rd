\name{bagboost}
\alias{bagboost}
\title{bagboost}
\description{{An implementation of the BagBoost algorithm for binary
  classification}
}
\usage{bagboost(xlearn, ylearn, xtest, presel = 200, mfinal = 100, bag = 50)}
\arguments{
  \item{xlearn}{A (n x p)-matrix, where rows correspond to training
    instances and columns contain the predictor variables.}
  \item{ylearn}{A vector of length n containing the class labels, which
    need to be coded by 0 and 1.}
  \item{xtest}{A (m x p)-matrix, where rows correspond to test
      instances and columns contain the predictor variables.}
  \item{presel}{An integer, giving the number of features to be
    pre-selected according to the Wilcoxon test statistic. Default is
    presel=200 features. If presel=0, no feature preselection is carried
    out.}
  \item{mfinal}{An integer, the number of iterations for which boosting
    is run. Defaults to mfinal=100 iterations}
  \item{bag}{An integer, the number of bagging steps that shall be done
    to obtain the weak learner. Defaults to bag=50 bagging iterations.}
}

\value{The function outputs an array, whose rows contain out-of-sample
  probabilities that the class labels are predicted as being of class 1,
  for every boosting iteration.}
\references{
  \item{o}{"Boosting for Tumor Classification with Gene Expression
    Data", Marcel Dettling and Peter B\"uhlmann. Bioinformatics (2003),
    Vol. 19, p. 1061--1069.}
  \item{o}{"BagBoosting for Tumor Classification with Gene Expression
    Data", Marcel Dettling. To appear in Bioinformatics (2005).}
  \item{o}{Further information is available from the webpage
    \url{http://stat.ethz.ch/~dettling}}
  }
\author{Marcel Dettling}
\examples{
data(leukemia, package = "boost")

## Dividing the leukemia dataset into training and test data
xlearn <- leukemia.x[c(1:20, 34:38),]
ylearn <- leukemia.y[c(1:20, 34:38)]
xtest  <- leukemia.x[21:33,]
ytest  <- leukemia.y[21:33]

## Classification with bagboost
fit <- bagboost(xlearn, ylearn, xtest, presel=50, mfinal=20, bag=5)
summarize(fit, ytest)
}
\keyword{classif}

