# Tests for association (or independence) between one categorical variable (like
# a SNP) and a binary phenotype, conditional on a second categorical variable 
# (like sex or ancestry, or another SNP).

# NOTE: In all matrix inputs and outputs (e.g., contingency tables, known 
# pairwise distributions, initial and fitted parameter values) the variable
# being conditioned on is the row (first coordinate) variable. It will be named
# x1 here. The covariate being conditioned on is referred to as x2.

# Standard saturated logistic regression based test
standard.conditional.test = function(t0, t1) {
  yg = cbind(c(t1), c(t0))
  xg = expand.grid(x1 = factor(0:(nrow(t0) - 1)), x2 = factor(0:(ncol(t0) - 1)))
  
  # FIXME this could be done faster since both are saturated models (of their contingency respective tabels)
  null.model   = glm(yg ~ xg$x1        , family = binomial)
  altern.model = glm(yg ~ xg$x1 * xg$x2, family = binomial)
  
  anv = anova(null.model, altern.model, test = 'Chisq') # Why not F?
  
  pen = matrix(predict(altern.model, type = 'response'), nrow(t0))
  statistic = anv$Deviance[2]
  p.value = anv$'Pr(>Chi)'[2]
  
  return (list(pen = pen, statistic = statistic, p.value = p.value))
}

cmle.conditional.test = function(t0, t1, tp, prevalence, pen.initial = NULL, pxx.initial = NULL) {
  # 1. Fit the full model assuming x2 has no effect. This is equivalent to fitting a marginal model
  # with x1 only and forgetting entirely about x2 (because the constraint actually does not involve 
  # x2, and its effect on the likelihood cancels from the likelihood ratio statistic)
  if (is.null(pen.initial)) {
    pen.initial.x1 = NULL
    pxx.initial.x1 = NULL
  } else {
    pxx.initial.x1 = rowSums(pxx.initial)
    pen.initial.x1 = rowSums(pen.initial * pxx.initial) / pxx.initial.x1
  }
  
  fit.x1 = marginal.assoc.test.pop.kpy(rowSums(t0), rowSums(t1), rowSums(tp),
    prevalence, pen.initial.x1, pxx.initial.x1)
  
  # 2. Fit the full model using x1 and x2
  if (is.null(pen.initial)) {
    pen.initial = cbind(fit.x1$pen, fit.x1$pen, fit.x1$pen)
    px1.initial = fit.x1$px
    tx = t0 + t1 + tp
    px2.initial = colSums(tx) / sum(tx)
    pxx.initial = px1.initial %o% px2.initial
    pxx.initial[pxx.initial < .Machine$double.eps] = .Machine$double.eps
  }
  
  fit.x1x2 = pairwise.assoc.test.pop.kpy(t0, t1, tp, prevalence, pen.initial, pxx.initial)
  
  # 3. Perform the GLRT for 1. vs 2.
  statistic = fit.x1x2$statistic - fit.x1$statistic
  p.value = pchisq(statistic, df = length(t0) - nrow(t0), lower.tail = F) # FIXME what happens with empty cells?
  
  return (list(statistic = statistic, p.value = p.value))
}

# Assuming the tested SNP is in HWE
cmle.conditional.test.hwe = function(t0, t1, tp, prevalence, pen.initial = NULL, p.x1.initial = NULL, f2.given.x1.initial = NULL) {
  # 1. Fit the full model using x1 alone
  if (is.null(pen.initial)) {
    pen.initial.x1 = NULL
  } else {
    pxx.initial = p.x1.initial * cbind((1 - f2.given.x1.initial) ^ 2, 2 * f2.given.x1.initial * (1 - f2.given.x1.initial), f2.given.x1.initial ^ 2)
    pen.initial.x1 = rowSums(pen.initial * pxx.initial) / p.x1.initial
  }
  
  fit.x1 = marginal.assoc.test.pop.kpy(rowSums(t0), rowSums(t1), rowSums(tp),
    prevalence, pen.initial.x1, p.x1.initial)
  
  # 2. Fit the full model using x1 and x2
  if (is.null(pen.initial)) {
    pen.initial = cbind(fit.x1$pen, fit.x1$pen, fit.x1$pen)
    p.x1.initial = fit.x1$px
    tx = t0 + t1 + tp
    tx1 = rowSums(tx)
    f2.given.x1.initial = (0.5 * tx[, 2] + tx[, 3]) / tx1
    f2.given.x1.initial[f2.given.x1.initial > 0.5] = 0.5
    f2.given.x1.initial[f2.given.x1.initial == 0] = .Machine$double.eps
  }
  
  fit.x1x2 = pairwise.assoc.test.pop.hwe2.kpy(t0, t1, tp, prevalence, pen.initial, p.x1.initial, f2.given.x1.initial)
  
  # 3. Perform the GLRT for 1. vs 2.
  statistic = fit.x1x2$statistic - fit.x1$statistic
  p.value = pchisq(statistic, df = length(t0) - nrow(t0), lower.tail = F) # FIXME what happens with empty cells?
  
  return (list(statistic = statistic, p.value = p.value))
}

# Assuming independence of the test locus and the covariate being conditioned on.
cmle.conditional.test.ind = function(t0, t1, tp, prevalence, pen.initial = NULL, px1.initial = NULL, px2.initial = NULL) {
  # 1. Fit the full model using x1 alone (whatever x1, the row variable, stands for)
  if (is.null(pen.initial)) {
    pen.initial.x1 = NULL
  } else {
    pen.initial.x1 = rowSums(pen.initial * (px1.initial %o% px2.initial)) / px1.initial
  }
  
  fit.x1 = marginal.assoc.test.pop.kpy(rowSums(t0), rowSums(t1), rowSums(tp),
    prevalence, pen.initial.x1, px1.initial)
  
  # 2. Fit the full model using x1 and x2
  if (is.null(pen.initial)) {
    pen.initial = cbind(fit.x1$pen, fit.x1$pen, fit.x1$pen)
    px1.initial = fit.x1$px
    tx = t0 + t1 + tp
    px2.initial = colSums(tx) / sum(tx)
  }
  
  fit.x1x2 = pairwise.assoc.test.pop.le.kpy(t0, t1, tp, prevalence, pen.initial, px1.initial, px2.initial)
  
  # 3. Perform the GLRT for 1. vs 2.
  statistic = fit.x1x2$statistic - fit.x1$statistic
  p.value = pchisq(statistic, df = length(t0) - nrow(t0), lower.tail = F)
  
  return (list(statistic = statistic, p.value = p.value))
}

# Test SNP (x2) effect conditional on the effect of the covariate (x1), assuming x1 and x2 are 
# independent in the population, and that x1 is in HWE.
cmle.conditional.test.hwe.ind = function(t0, t1, tp, prevalence, pen.initial = NULL, px1.initial = NULL, f2.initial = NULL) {
  # 1. Fit the full model using x1 alone
  if (is.null(pen.initial)) {
    pen.initial.x1 = NULL
  } else {
    px2.initial = c((1 - f2.initial) ^ 2, 2 * f2.initial * (1 - f2.initial), f2.initial ^ 2)
    pen.initial.x1 = rowSums(pen.initial * (px1.initial %o% px2.initial)) / px1.initial
  }
  
  fit.x1 = marginal.assoc.test.pop.kpy(rowSums(t0), rowSums(t1), rowSums(tp),
    prevalence, pen.initial.x1, px1.initial)
  
  # 2. Fit the full model using x1 and x2
  if (is.null(pen.initial)) {
    pen.initial = cbind(fit.x1$pen, fit.x1$pen, fit.x1$pen)
    px1.initial = fit.x1$px
    tx2 = colSums(t0 + t1 + tp)
    f2.initial = min((0.5 * tx2[2] + tx2[3]) / sum(tx2), 0.5)
  }
  
  fit.x1x2 = pairwise.assoc.test.pop.hwe2.ind.kpy(t0, t1, rowSums(tp), colSums(tp), 
    prevalence, pen.initial, px1.initial, f2.initial)
  
  # 3. Perform the GLRT for 1. vs 2.
  statistic = fit.x1x2$statistic - fit.x1$statistic
  p.value = pchisq(statistic, df = length(t0) - nrow(t0), lower.tail = F)
  
  return (list(statistic = statistic, p.value = p.value))
}

# Test one SNP (x2) conditional on the effect of another SNP (x1), assuming they are in LE and that 
# both are in HWE.
cmle.conditional.test.hwe.le = function(t0, t1, tp, prevalence, pen.initial = NULL, f1.initial = NULL, f2.initial = NULL) {
  # 1. Fit the full model using x1 alone (x1, the row variable, is assumed to be a ternari SNP in HWE)
  if (is.null(pen.initial)) {
    pen.initial.x1 = NULL
  } else {
    px1.initial = c((1 - f1.initial) ^ 2, 2 * f1.initial * (1 - f1.initial), f1.initial ^ 2)
    px2.initial = c((1 - f2.initial) ^ 2, 2 * f2.initial * (1 - f2.initial), f2.initial ^ 2)
    pen.initial.x1 = rowSums(pen.initial * (px1.initial %o% px2.initial)) / px1.initial
  }
  
  fit.x1 = marginal.assoc.test.pop.hwe.kpy(rowSums(t0), rowSums(t1), rowSums(tp),
    prevalence, pen.initial.x1, f1.initial)
  
  # 2. Fit the full model using x1 and x2
  if (is.null(pen.initial)) {
    pen.initial = cbind(fit.x1$pen, fit.x1$pen, fit.x1$pen)
    f1.initial = fit.x1$f
    tx2 = colSums(t0 + t1 + tp)
    f2.initial = min((0.5 * tx2[2] + tx2[3]) / sum(tx2), 0.5)
  }
  
  fit.x1x2 = pairwise.assoc.test.pop.hwe.le.kpy(t0, t1, rowSums(tp), colSums(tp), 
    prevalence, pen.initial, f1.initial, f2.initial)
  
  # 3. Perform the GLRT for 1. vs 2.
  statistic = fit.x1x2$statistic - fit.x1$statistic
  p.value = pchisq(statistic, df = length(t0) - nrow(t0), lower.tail = F)
  
  return (list(statistic = statistic, p.value = p.value))
}

cmle.conditional.known.pxx = function(t0, t1, prevalence, p.xx, pen.initial = NULL) {
  p.x1 = rowSums(p.xx)
  
  # 1. Fit the full model using x1 alone (whatever x1, the row variable stands for)
  if (is.null(pen.initial)) {
    pen.initial.x1 = NULL
  } else {
    pen.initial.x1 = rowSums(pen.initial * p.xx) / p.x1
  }
  
  fit.x1 = marginal.assoc.test.kpx.kpy(rowSums(t0), rowSums(t1), prevalence, p.x1, pen.initial.x1)
  
  # 2. Fit the full model using x1 and x2
  fit.x1x2 = pairwise.assoc.test.kpx.kpy(t0, t1, prevalence, p.xx, pen.initial)
  
  # 3. Perform the GLRT for 1. vs 2.
  statistic = fit.x1x2$statistic - fit.x1$statistic
  p.value = pchisq(statistic, df = length(t0) - nrow(t0), lower.tail = F) # FIXME what happens with empty cells?
  
  return (list(statistic = statistic, p.value = p.value))
}
