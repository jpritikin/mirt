#' Full information maximum likelihood estimation of multidimensional IRT models
#'
#' Analysis of dichotomous and polytomous response data using latent
#' trait models under the Item Response Theory paradigm. Exploratory models can be
#' estimated via quadrature or stochastic methods, a generalized confirmatory bi-factor
#' analysis is included, and confirmatory models can be fit with a
#' Metropolis-Hastings Robbins-Monro algorithm which can include polynomial or
#' product constructed latent traits. Multiple group analysis
#' and mixed effects designs may be performed for unidimensional or
#' multidimensional item response models for detecting differential item
#' functioning and modelling item and person covariates.
#'
#' Users interested in the most recent version of this package can visit
#' \code{https://github.com/philchalmers/mirt} and follow the instructions
#' for installing the package from source.
#'
#'
#'
#' @name mirt-package
#' @docType package
#' @title Full information maximum likelihood estimation of IRT models.
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @useDynLib mirt
#' @importFrom stats anova fitted residuals
#' @import lattice GPArotation mvtnorm
#' @exportMethod anova
#' @exportMethod fitted
#' @exportMethod residuals
#' @exportMethod summary
#' @keywords package
NULL

#' Description of Science data
#'
#' A 4 item data set borrowed from \code{\link[ltm]{ltm}} package, first example
#' of the \code{\link[ltm]{grm}} function. See more complete documentation therein.
#'
#'
#' @name Science
#' @docType data
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords data
NULL

#' Description of SAT12 data
#'
#' Data obtained from the TESTFACT (Woods et al., 2003) manual, with 32 response pattern
#' scored items for a grade 12 science assessment test (SAT) measuring topics of chemistry,
#' biology, and physics. The scoring key for these data is
#' [1, 4, 5, 2, 3, 1, 2, 1, 3, 1, 2, 4, 2, 1, 5, 3, 4, 4, 1, 4, 3, 3, 4, 1, 3, 5, 1, 3, 1, 5, 4, 5],
#' respectively.
#'
#'
#' @name SAT12
#' @docType data
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Wood, R., Wilson, D. T., Gibbons, R. D., Schilling, S. G., Muraki, E., & Bock, R. D. (2003).
#' TESTFACT 4 for Windows: Test Scoring, Item Statistics, and Full-information Item Factor Analysis
#' [Computer software]. Lincolnwood, IL: Scientific Software International.
#'
#' @keywords data
NULL

#' Description of LSAT7 data
#'
#' Data from Bock & Lieberman (1970); contains 5 dichotomously scored
#' items obtained from the Law School Admissions Test, section 7.
#'
#'
#' @name LSAT7
#' @docType data
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Bock, R. D., & Lieberman, M. (1970). Fitting a response model for \emph{n}
#' dichotomously scored items. \emph{Psychometrika, 35}(2), 179-197.
#'
#' @keywords data
NULL

#' Description of LSAT6 data
#'
#' Data from Thissen (1982); contains 5 dichotomously scored
#' items obtained from the Law School Admissions Test, section 6.
#'
#'
#' @name LSAT6
#' @docType data
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Thissen, D. (1982). Marginal maximum likelihood estimation for the one-parameter logistic model.
#' \emph{Psychometrika, 47}, 175-186.
#'
#' @keywords data
NULL

#' Description of deAyala data
#'
#' Mathematics data from de Ayala (2009; pg. 14); 5 item dataset in table format.
#'
#'
#' @name deAyala
#' @docType data
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' de Ayala, R. J. (2009). \emph{The theory and practice of item response theory}. Guilford Press.
#'
#' @keywords data
NULL
