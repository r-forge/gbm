\name{estBisqParam}
\Rdversion{1.1}
\alias{estBisqParam}
\title{
Estimate the Bisquare Parameter
}
\description{
Converts an efficiency percentage into the associated bisquare parameter.
}
\usage{
estBisqParam(x = 0.85)
}
\arguments{
  \item{x}{
bisquare efficiency. Only bisquare efficiency between 0.8 and 0.95 is currently supported
}
}
\details{
Borrowed from \code{lmRob.fit.compute()} in the robust package
}
\value{
Returns the parameter for the bisquare last function for the given efficiency
}
\references{
See the robust package
}
\author{
Harry Southworth
}

\seealso{
\code{\link{lmRob.fit.compute}}
}
\examples{
estBisqParam(x = 0.85) 
}

\keyword{robust }
