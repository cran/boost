\name{leukemia}
\alias{leukemia}
\alias{leukemia.x}
\alias{leukemia.y}
\alias{leukemia.z}
\non_function{}
\title{A part of the famous AML/ALL-leukemia dataset}
\description{This is the training set of the famous AML/ALL-leukemia
  dataset from the Whitehead Institute. It has been reduced to 250
  genes, about the half of which are very informative for
  classification, whereas the other half was chosen randomly.}
\usage{data(leukemia)}
\format{Contains three R-objects: The expression matrix leukemia.x, the
  associated binary response variable leukemia.y, and the associated
  3-class response variable leukemia.z}
\source{http://www.genome.wi.mit.edu/MPR}
\references{First published in Golub et al: Molecular Classification of
  Cancer: Class Discovery and Class Prediction by Gene Expression
  Monitoring. Science 1999, 286: 531-538.}
\examples{
data(leukemia)
str(leukemia.x)
str(leukemia.y)
str(leukemia.z)
par(mfrow=c(1,2))
plot(leukemia.x[,56], leukemia.y)
plot(leukemia.x[,174],leukemia.z)}
\keyword{datasets}



