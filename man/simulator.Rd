\name{simulator}
\alias{simulator}
\title{simulator}
\description{Simulation of (microarray) data according to correlation
  and mean structures from real datasets.}
\usage{simulator(x, y, respmod = c("none", "resp1", "resp2", "resp3"),
nos = 1200, gene = NULL, signs = NULL)}
\arguments{
  \item{x}{A (n x p)-matrix, whose correlation and mean structure is to
    be used for simulating data. Its rows correspond to training
    instances and columns contain the predictor variables.}
  \item{y}{A vector of length n containing the class labels, which
    need to be coded by 0 and 1.}
  \item{respmod}{A character string. Either "none" where the simulated
    gene expression labels are determined model-free depending which
    class mean and correlation structure had been used for their
    determination. The choice of "resp1", "resp2" and "resp3"
    means that a response model is applied. For "resp1", 10 genes are
    selected and determine conditional proabilities via a logistic model
    with equal weights. The class labels are then regarded as having a
    Bernoulli distribution with probability p. For "resp2", 25 genes
    are plugged into the logistic model with non-equal weights. With
    "resp3", 25 genes are chosen for a logistic model with second
    and third order interactions.}
  \item{nos}{An integer, giving the number of instances which are
    simulated.}
  \item{gene}{A vector giving the index of the genes which shall be used
    for model based class label simulation. Defaults to NULL. This
    argument should only be used for specially designed simulation
    studies, where it is important that the same predictor variables are
    repeatedly used for simulating class label.}
  \item{signs}{A vector containing entries of +1 and -1. Defaults to
    NULL and is only of importance in specially designed simulation
    studies, where it is important that the same predictor variables are
    repeatedly used for simulating class label. }
}
\details{The new instances are simulated according to a multivariate
  normal distribution with means and correlation structure taken from a
  real (gene expression) dataset. This structure is obtained by
  transforming a standard multivariate normal distribution, which
  requires a eigenvalue decomposition of the provided real dataset. For
  datasets with many predictors (>500), this can be fairly time
  consuming. Simulating data without applying a response model is fine
  for most purposes, only special analysis tasks usually require it.}
\value{Returns a list containing
  \item{x}{An (nos x p)-matrix, containing the simulated data}
  \item{y}{A vector of length nos, containing the class labels of the
    simulated data.}
  \item{probab}{A vector of length nos, containing the conditional
    probabilities of the simulated data. Is empty if respmod="none".}
  \item{bayes}{An integer, giving the Bayes error (theoretically minimal
    misclassification risk) for the simulated data. Is empty if 
    respmod="none".}
  \item{gene}{A vector, containing the indices of the variables which
    had been used in the logistic model for either "resp1", "resp2" or
    "resp3". Is empty if respmod="none".}
  \item{signs}{A vector, containing -1 and +1. Indicates with what
    polarization a predictor variable had been used in the logistic
    model. Is empty if respmod="none".}b
}
\references{
  \item{o}{"BagBoosting for Tumor Classification with Gene Expression
    Data", Marcel Dettling. To appear in Bioinformatics (2005).}
  \item{o}{Further information is available from the webpage
    \url{http://stat.ethz.ch/~dettling}}
}
\examples{
set.seed(21)
data(leukemia)

## Simulation of gene expression data
simu <- simulator(leukemia.x, leukemia.y, nos=200)

## Defining training and test data
xlearn <- simu$x[1:150,]
ylearn <- simu$y[1:150]
xtest  <- simu$x[151:200,]
ytest  <- simu$y[151:200]

## Classification with logitboost
fit <- logitboost(xlearn, ylearn, xtest, mfinal=20, presel=50)
summarize(fit, ytest)
}
\keyword{classif}
