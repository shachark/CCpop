\name{pairwise.assoc.test.unconstrained.chisq}
\alias{pairwise.assoc.tests}
\alias{pairwise.assoc.test.unconstrained.chisq}
\alias{pairwise.assoc.test.unconstrained.gsq}
\alias{pairwise.assoc.test.case.only}
\alias{pairwise.assoc.test.ind.3d}
\alias{pairwise.assoc.test.pure.unconstrained}
\alias{pairwise.assoc.test.kpy}
\alias{pairwise.assoc.test.hwe.le.kpy}
\alias{pairwise.assoc.test.kpx.kpy}
\alias{pairwise.assoc.test.pop.kpy}
\alias{pairwise.assoc.test.pop.hwe.le.kpy}
\alias{pairwise.assoc.test.pure.pop.kpy}
\alias{pairwise.assoc.test.pure.pop.hwe.le.kpy}
\alias{conditional.assoc.test.pure.pop.hwe.le.kpy}

\title{Pairwise association tests}
\description{Tests for joint association of a pair of SNP loci and a binary phenotype in case-only, case-control, case-population, and case-control-population designs, under various assumptions.}

\usage{
  pairwise.assoc.test.unconstrained.chisq(t0, t1)
  pairwise.assoc.test.unconstrained.gsq(t0, t1)
  pairwise.assoc.test.case.only(t1)
  pairwise.assoc.test.ind.3d(t0, t1)
  pairwise.assoc.test.pure.unconstrained(t0, t1)
  pairwise.assoc.test.kpy(t0, t1, prevalence, pen.initial = NULL, pxx.initial = NULL)
  pairwise.assoc.test.hwe.le.kpy(t0, t1, prevalence, 
    pen.initial = NULL, f1.initial = NULL, f2.initial = NULL)
  pairwise.assoc.test.kpx.kpy(t0, t1, prevalence, pxx, pen.initial = NULL)
  pairwise.assoc.test.pop.kpy(t0, t1, tp, prevalence, 
    pen.initial = NULL, pxx.initial = NULL)
  pairwise.assoc.test.pop.hwe.le.kpy(t0, t1, tp1, tp2, prevalence, 
    pen.initial = NULL, f1.initial = NULL, f2.initial = NULL)
  pairwise.assoc.test.pure.pop.kpy(t0, t1, tp, prevalence, 
    pen.initial = NULL, pxx.initial = NULL)
  pairwise.assoc.test.pure.pop.hwe.le.kpy(t0, t1, tp1, tp2, prevalence, 
    pen.initial = NULL, f1.initial = NULL, f2.initial = NULL)
  conditional.assoc.test.pure.pop.hwe.le.kpy(t0, t1, tp1, tp2, prevalence, 
    pen.initial = NULL, f1.initial = NULL, f2.initial = NULL)
}

\arguments{
  \item{t0}{A 3x3 matrix of genotype counts in the controls cohort}
  \item{t1}{A 3x3 matrix of genotype counts in the cases cohort}
  \item{tp}{A 3x3 matrix of genotype counts in the population cohort}
  \item{tp1}{A 3-component vector of genotype counts for the first (\code{t0/t1} rows) locus in the population cohort}
  \item{tp2}{A 3-component vector of genotype counts for the second (\code{t0/t1} columns) locus in the population cohort}
  \item{prevalence}{The population prevalence of the phenotype \code{Pr(y=1) < 0.5}}
  \item{pxx}{A 3x3 matrix of known genotypic distribution in the population \code{Pr(x1,x2)}}
  \item{pen.initial}{A 3x3 matrix of initial guess for the phenotypic penetrance \code{Pr(y=1|x1,x2)}}
  \item{pxx.initial}{A 3x3 matrix of initial guess for the genotypic distribution in the population \code{Pr(x1,x2)}}
  \item{f1.initial}{Initial guess for the MAF at the first (rows) locus}
  \item{f2.initial}{Initial guess for the MAF at the second (columns) locus}
}

\details{
  pairwise.assoc.test.unconstrained.chisq -- Pearson's chi-squared test for the 2x9 contingency table of a standard case-control design.
  
  pairwise.assoc.test.unconstrained.gsq -- Likelihood ratio ("G") test for the 2x9 contingency table of a standard case-control design.
  
  pairwise.assoc.test.case.only -- The likelihood ratio ("G") test for the 3x3 contingency table of the two SNPs in cases only.
  
  pairwise.assoc.test.ind.3d -- The likelihood ratio ("G") test for the 2x3x3 contingency table of the two SNPs and phenotype value (i.e., a case-control design).
  
  pairwise.assoc.test.pure.unconstrained -- The logistic regression based GLRT for full model vs. main effects model (within a standard case-control design).
  
  pairwise.assoc.test.kpy -- The GLRT based on constrained maximum likelihood estimation assuming that phenotype prevalence in the population is known (within a standard case-control design).
  
  pairwise.assoc.test.hwe.le.kpy -- The GLRT based on constrained maximum likelihood estimation assuming Hardy-Weinberg equilibrium (HWE) and linkage equilibrium (LE) in the general population, and that phenotype prevalence in the population is known (within a standard case-control design).
  
  pairwise.assoc.test.kpx.kpy -- The GLRT based on constrained maximum likelihood estimation assuming that marginal genotypic and phenotypic distributions in the population are both known (within a standard case-control design).
  
  pairwise.assoc.test.pop.kpy -- The GLRT based on constrained maximum likelihood estimation assuming that phenotype prevalence in the population is known (within a case-control-population design).
  
  pairwise.assoc.test.pop.hwe.le.kpy -- The GLRT based on constrained maximum likelihood estimation assuming Hardy-Weinberg equilibrium (HWE) and linkage equilibrium (LE) in the general population, and that phenotype prevalence in the population is known (within a case-control-population design).
  
  pairwise.assoc.test.pure.pop.kpy -- A GLRT based on constrained maximum likelihood estimation assuming that phenotype prevalence is known (within a case-control-population design). The test contrasts a general pairwise association versus a logit-additive "main-effects" model.
  
  pairwise.assoc.test.pure.pop.hwe.le.kpy -- Like pairwise.assoc.test.pure.pop.kpy, but done under the assumption of HWE and LE.
  
  conditional.assoc.test.pure.pop.hwe.le.kpy -- A GLRT based on constrained maximum likelihood estimation assuming that phenotype prevalence is known and assuming HWE and LE (within a case-control-population design). This tests for a pairwise association above and beyond a marginal association of x1 alone, or of x2 alone (and the *maximum* p-value of these two options is returned, i.e., the test is conservative)
}

\value{
  \item{p.value}{P-value indicating significance of the association}
  \item{statistic}{Value of the relevant test statistic}
  \item{pen}{The estimated case-control peneterance \code{Pr(y=1|x,ascertained)}}
}

\references{
  Kaufman, S., & Rosset, S. (2014). Exploiting Population Samples To Enhance Genome-Wide Association Studies of Disease. Genetics, 114 (early electronic).
}

\author{
  Shachar Kaufman
}

\note{
  The implementation assumes that contingency tables are provided with homozygote major allele as the first row and column.
  
  In order to perform case-population testing, the user may set \code{t0} to a matrix of zeros. 
  
  An inferior alternative to case-control-population analysis which may be desirable in the interest of analysis speed and can be acceptable for low prevalence phenotypes is to extend the controls cohort with any population samples. This can be achieved by providing \code{t0.new = t0.old + tp.old} and using one of the case-control tests.
  
  \code{pen} is not returned by the case-only test.
}

\seealso{
  \code{\link{marginal.assoc.tests}},
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
