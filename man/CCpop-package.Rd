\name{CCpop-package}
\alias{CCpop-package}
\alias{CCpop}

\docType{package}

\title{One and two locus GWAS of binary phenotype with case-control-population design}

\description{Tests for marginal and pairwise SNP associations with binary phenotypes in case-control, case-population, and case-control-population studies.}

\details{
  \tabular{ll}{
    Package: \tab CCpop\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2014-03-18\cr
    License: \tab GPL-2\cr
  }

  Inputs to test functions are vectors (for marginal tests) and matrices (for pairwise tests) of genotypic counts and are given separately for cases, controls, and population cohorts (where relevant). The constrained maximum likelihood estimation based tests (Kaufman and Rosset, 2014) also require a value for the known phenotypic prevalence in the population (a crude estimate is usually sufficient). When applicable, tests assuming Hardy-Weinberg equilibrium (HWE) and linkage equilibrium (LE) will be considerably more powerful.
  
  See \link{marginal.assoc.tests} for single SNP tests, and \link{pairwise.assoc.tests} for pairwise joint, pure-interaction, and conditional tests.
}

\author{Shachar Kaufman <shachark@post.tau.ac.il>}

\references{
Kaufman, S., & Rosset, S. (2014). Exploiting Population Samples To Enhance Genome-Wide Association Studies of Disease. Genetics, genetics-114.
}

\keyword{
  GWAS
  case-control
  case-population
  case-control-population
  Hardy-Weinberg equilibrium
  Linkage equilibrium
}

%\seealso{
% ~~ Optional links to other man pages, e.g. ~~
% ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}

\examples{
  ## An example marginal/pairwise association
  
  # Controls
  t0 = matrix(c(375,  240,   46,
                640,  405,   62,
                300,  169,   19), nrow = 3, byrow = TRUE)
                
  # Cases
  t1 = matrix(c(317,  162,   15,
                459,  209,   22,
                120,   76,   13), nrow = 3, byrow = TRUE)
  
  # Independent population sample, marginalized for SNP1 and SNP2
  tp1 = c(2410, 4253, 1945)
  tp2 = c(4972, 3140,  496)
  
  ## The prevalence of the studied disease in the population
  prevalence = 0.001

  marginal.assoc.test.pop.hwe.kpy(t0 = rowSums(t0), t1 = rowSums(t1), tp = tp1, prevalence)
  marginal.assoc.test.pop.hwe.kpy(t0 = colSums(t0), t1 = colSums(t1), tp = tp2, prevalence)
  pairwise.assoc.test.pop.hwe.le.kpy(t0, t1, tp1, tp2, prevalence)  
  conditional.assoc.test.pure.pop.hwe.le.kpy(t0, t1, tp1, tp2, prevalence)
}
