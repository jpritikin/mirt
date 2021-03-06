\name{personfit}
\alias{personfit}
\title{Person fit statistics}
\usage{
  personfit(x, method = "EAP", ...)
}
\arguments{
  \item{x}{a computed model object of class
  \code{ExploratoryClass}, \code{ConfirmatoryClass}, or
  \code{MultipleGroupClass}}

  \item{method}{type of factor score estimation method. Can
  be expected a-posteriori (\code{"EAP"}), Bayes modal
  (\code{"MAP"}), weighted likelihood estimation
  (\code{"WLE"}), or maximum likelihood (\code{"ML"})}

  \item{...}{additional arguments to be passed to
  \code{fscores()}}
}
\description{
  \code{personfit} calculates the Zh values from Drasgow,
  Levine and Williams (1985) for unidimensional and
  multidimensional models. For Rasch models infit and
  outfit statistics are also produced. The returned object
  is a \code{data.frame} consisting either of the tabulated
  data or full data with the statistics appended to the
  rightmost columns.
}
\examples{
\dontrun{

#make some data
set.seed(1234)
a <- matrix(rlnorm(20),ncol=1)
d <- matrix(rnorm(20),ncol=1)
items <- rep('dich', 20)
data <- simdata(a,d, 2000, items)

x <- mirt(data, 1)
fit <- personfit(x)
head(fit)

  }
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\references{
  Drasgow, F., Levine, M. V., & Williams, E. A. (1985).
  Appropriateness measurement with polychotomous item
  response models and standardized indices. \emph{Journal
  of Mathematical and Statistical Psychology, 38}, 67-86.

  Reise, S. P. (1990). A comparison of item- and person-fit
  methods of assessing model-data fit in IRT. \emph{Applied
  Psychological Measurement, 14}, 127-137.
}
\seealso{
  \code{\link{itemfit}}
}
\keyword{fit}
\keyword{person}

