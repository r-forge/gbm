\name{print.gbm}
\alias{print.gbm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  Display basic information about a \code{gbm} object.
}
\usage{
print.gbm(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class \code{gbm}. }
  \item{\dots}{ arguments passed to \code{print.default}. }
}
\details{
  Prints some information about the model object. In particular, the call
  to the model fitting function is given, and the type of loss function
  that was used is given, as is the total number of iterations.
  
  If cross-validation was performed, the 'best' number of trees as
  estimated
   by cross-validation error is
  dispalyed. If a test set was used, the 'best' number
  of trees as estimated by the test set error is displayed.
  
  The number of available predictors, and the number of those having
  non-zero influence on predictions is given (which might be interesting
  in data mining applications).
  
  If K-class, Bernoulli or adaboost classification was performed,
  the confusion matrix and prediction accuracy are printed (objects
  being allocated to the class with highest probability for K-class
  and Bernoulli). These classifications are performed on the entire
  training
  data using the model with the 'best' number of trees as described
  above, or the maximum number of trees if the 'best' can't be
  computed.
  
  If the 'distribution' was specified as Gaussian, Laplace, quantile,
  bisquare or t-distribution, a summary of the residuals is displayed.
  The residuals are for the training data with the model at the 'best' number of trees, as
  described above, or the maximum number of trees if the 'best' can't
  be computed.
  
}


\author{ Harry Southworth, Daniel Edwards }

\seealso{ \code{\link{gbm}} }
\examples{
library( gbm )
data( iris )
iris.mod <- gbm( Species ~ ., distribution="kclass", data=iris,
                 n.trees=2000, shrinkage=.01, cv.folds=5 )
iris.mod
data( lung )
lung.mod <- gbm( Surv(time, status) ~ ., distribution="coxph", data=lung,
                 n.trees=2000, shrinkage=.01, cv.folds=5 )
lung.mod
}
\keyword{models}
\keyword{nonlinear}
\keyword{survival}
\keyword{nonparametric}
\keyword{tree}
