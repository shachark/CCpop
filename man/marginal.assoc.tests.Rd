\name{marginal.assoc.test.unconstrained.chisq}
\alias{marginal.assoc.tests}
\alias{marginal.assoc.test.unconstrained.chisq}
\alias{marginal.assoc.test.unconstrained.gsq}
\alias{marginal.assoc.test.hwe.in.controls}
\alias{marginal.assoc.test.pop.hwe.kpy}
\alias{marginal.assoc.test.kpx.kpy}

\title{Marginal association tests}
\description{Tests for association of single SNP locus and a binary phenotype in case-control, case-population, and case-control-population designs, under various assumptions.}

\usage{
  marginal.assoc.test.unconstrained.chisq(t0, t1)
  marginal.assoc.test.unconstrained.gsq(t0, t1)
  marginal.assoc.test.hwe.in.controls(t0, t1)
  marginal.assoc.test.pop.hwe.kpy(t0, t1, tp, prevalence, 
    pen.initial = NULL, f.initial = NULL)
  marginal.assoc.test.kpx.kpy(t0, t1, prevalence, px, pen.initial = NULL)
}

\arguments{
  \item{t0}{A 3-component vector of genotype counts in the controls cohort}
  \item{t1}{A 3-component vector of genotype counts in the cases cohort}
  \item{tp}{A 3-component vector of genotype counts in the population cohort}
  \item{prevalence}{The population prevalence of the phenotype \code{Pr(y=1) < 0.5}}
  \item{px}{A 3-component vector of known genotypic distribution in the population \code{Pr(x)}}
  \item{pen.initial}{A 3-component vector of initial guess for phenotypic penetrance \code{Pr(y=1|x)}}
  \item{f.initial}{Initial guess for the SNP's minor allele frequency (MAF)}
}

\details{
 marginal.assoc.test.unconstrained.chisq -- Pearson's chi-squared test for the 2x3 contingency table of a standard case-control design.
 
 marginal.assoc.test.unconstrained.gsq -- Likelihood ratio ("G") test for the 2x3 contingency table of a standard case-control design.
 
 marginal.assoc.test.hwe.in.controls -- A test for association in a standard case-control design, assuming HWE among the controls (Chen and Chatterjee, 2007).
 
 marginal.assoc.test.pop.hwe.kpy -- GLRT based on constrained maximum likelihood estimation assuming HWE in the general population, and that disease prevalence is known.
 
 marginal.assoc.test.kpx.kpy -- GLRT based on constrained maximum likelihood estimation assuming the marginal genotypic and phenotypic distributions are known.
}

\value{
  \item{p.value}{P-value indicating significance of the association}
  \item{statistic}{Value of the relevant test statistic}
  \item{pen}{The estimated case-control peneterance \code{Pr(y=1|x,ascertained)}}
}

\references{
  Chen J., & Chatterjee, N. (2007). Exploiting hardy-weinberg equilibrium for effcient screening of single SNP associations from case-control studies. Human heredity, 63(3-4), 196-204.
  
  Kaufman, S., & Rosset, S. (2014). Exploiting Population Samples To Enhance Genome-Wide Association Studies of Disease. Genetics, 114 (early electronic).
}

\author{
  Shachar Kaufman
}

\note{
  The implementation assumes that contingency tables are provided with homozygote major allele as the first element.
  
  In order to perform case-population testing, the user may set \code{t0} to a vector of zeros. 
  
  An inferior alternative to case-control-population analysis which may be desirable in the interest of analysis speed and can be acceptable for low prevalence phenotypes is to extend the controls cohort with any population samples. This can be achieved by providing \code{t0.new = t0.old + tp.old} and using one of the case-control tests.
  
  Tests based on the 2x9 contingency table treat the two SNPs as one 9-category variable.
}

\seealso{
  \code{\link{pairwise.assoc.tests}},
}

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
