\name{summarize}
\alias{summarize}
\title{Summarize the output of classification with boosting functions}
\description{Yields text and graphical output that summarizes the
  misclassifcation error rates that have been achieved with boosting methods}
\usage{summarize(boost.out, resp, mout = ncol(boost.out), grafik = TRUE)}
\arguments{
  \item{boost.out}{An R-object, as obtained from one of the functions
    'bagboost', 'logitboost', 'adaboost' or 'l2boost'.}
  \item{resp}{A vector containing the class labels of the test
    instances. Needs to be coded by 0 and 1.}
  \item{mout}{The number of boosting iterations for which the error rate
  shall be printed. Defaults to the number of iterations boosting has
  been run for.}
  \item{grafik}{Logical, indicating whether a plot of the error rates is
  desired or not.}
}
\value{Just verbatim and graphical output.}
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
\seealso{bagboost, logitboost, adaboost, l2boost}
\examples{
data(leukemia, package = "boost")

## Dividing the leukemia dataset into training and test data
xlearn <- leukemia.x[c(1:20, 34:38),]
ylearn <- leukemia.y[c(1:20, 34:38)]
xtest  <- leukemia.x[21:33,]
ytest  <- leukemia.y[21:33]

## Classification with logitboost
fit <- logitboost(xlearn, ylearn, xtest, presel=50, mfinal=20)
summarize(fit, ytest)
}
\keyword{classif}
