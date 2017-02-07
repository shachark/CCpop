# Tests for association of one categorical covariate with the phenotype, above and beyond other
# covariates, where the latters can be continuous.

# TODO: give these better names within the package?

# The standard logistic regression GLRT. According to Prentice and Pyke (1979), if there is an
# intercept and if we assume nothing about the marginal probability of X, then this test is
# consistent and efficient under both the prospective and the retrospective designs.

prospective.logistic.glrt = function(formula.null, formula.altern, dat) {
  null.model = glm(formula.null, dat, family = binomial())
  altern.model = glm(formula.altern, dat, family = binomial())
  
  null.ll = logLik(null.model)
  altern.ll = logLik(altern.model)
  null.df = attr(null.ll, 'df')
  altern.df = attr(altern.ll, 'df')
  
  statistic = as.vector(2 * (altern.ll - null.ll))
  p.value = pchisq(statistic, df = altern.df - null.df, lower.tail = F)
  
  return (list(statistic = statistic, p.value = p.value))
}

# The (prospective) probit regression GLRT

prospective.probit.glrt = function(formula.null, formula.altern, dat) {
  null.model = glm(formula.null, dat, family = binomial(link = 'probit'))
  altern.model = glm(formula.altern, dat, family = binomial(link = 'probit'))
  
  null.ll = logLik(null.model)
  altern.ll = logLik(altern.model)
  null.df = attr(null.ll, 'df')
  altern.df = attr(altern.ll, 'df')
  
  statistic = as.vector(2 * (altern.ll - null.ll))
  p.value = pchisq(statistic, df = altern.df - null.df, lower.tail = F)
  
  return (list(statistic = statistic, p.value = p.value))
}

# The retrospective logistic model that exploits population samples and known prevalence.
# The distribution of explanatory variables is a nuisance parameter, and modeled semiparametrically.
# The computational approach of Scott & Wild (2001) is used.

# Helper: fit the model
fit.scott.wild = function(formula, data, K) {
  mf = model.frame(formula, data, na.action = NULL)
  y = model.extract(mf, 'response')
  current.na.action = options('na.action')
  options(na.action = 'na.pass')
  ax = model.matrix(mf, data)
  options(na.action = current.na.action$na.action)

  y0.idx = which(y == 0)
  y1.idx = which(y == 1)
  ym.idx = which(is.na(y))
  
  n0 = length(y0.idx)
  n1 = length(y1.idx)
  nm = length(ym.idx)
  n = n0 + n1 + nm
  
  w0 = nm + n0 / (1 - K)
  w1 = nm + n1 / K
  
  obj = function (pars) {
    eta = c(ax %*% pars)
    
    # Logit
    P1 = exp(eta) / (1 + exp(eta))
    P0 = 1 - P1
    p1 = exp(eta) / (1 + exp(eta)) ^ 2
    p0 = -p1
    
    # Probit
    #P0 = pnorm(eta, lower.tail = T)
    #P1 = pnorm(eta, lower.tail = F)
    #p0 = dnorm(eta)
    #p1 = -p0
    
    P = rep(1, n)
    P[y0.idx] = P0[y0.idx]
    P[y1.idx] = P1[y1.idx]
    Pw = w0 * P0 + w1 * P1
    
    p = rep(0, n)
    p[y0.idx] = p0[y0.idx]
    p[y1.idx] = p1[y1.idx]
    pw = w0 * p0 + w1 * p1
    
    obj = sum(log(Pw)) - sum(log(P))
    grad = c((pw / Pw - p / P) %*% ax)
    
    return (list(objective = obj, gradient = grad))
  }

  # Use impute-as-controls logistic regression to obtain initial param values
  y.imp = y
  y.imp[ym.idx] = 0
  pars.init = coef(glm.fit(ax, y.imp, family = binomial()))

  sol = nloptr(pars.init, eval_f = obj, opts = list(algorithm = 'NLOPT_LD_MMA', xtol_rel = 1e-3, maxeval = 1000))
  # FIXME check convergence
  
  ret = -sol$objective
  attr(ret, 'df') = length(pars.init)
  
  return (ret)
}

retrospective.logit.glrt = function(formula.null, formula.altern, dat, K) {
  null.ll = fit.scott.wild(formula.null, dat, K)
  altern.ll = fit.scott.wild(formula.altern, dat, K)
  null.df = attr(null.ll, 'df')
  altern.df = attr(altern.ll, 'df')
  
  statistic = 2 * (altern.ll - null.ll); attr(statistic, 'df') = NULL
  p.value = pchisq(statistic, df = altern.df - null.df, lower.tail = F)
  
  return (list(statistic = statistic, p.value = p.value))
}

# Helper: solve the extended C&C/S&W MLE under null and alternative
snw.cnc.core = function (ax, y, ax.expand, g.counts, null.idxs.dropped, K) {
  y0.idx = which(y == 0)
  y1.idx = which(y == 1)
  ym.idx = which(is.na(y))
  
  n = length(y)
  n0 = length(y0.idx)
  n1 = length(y1.idx)
  c0 = length(ym.idx)
  
  n.coeffs.altern = ncol(ax)
  n.coeffs.null = n.coeffs.altern - length(null.idxs.dropped)
  nr.g.levels = length(g.counts)
  
  g.idxs = rep(1:nr.g.levels, each = n)
  
  w0 = n0 / (1 - K)
  w1 = n1 / K
  
  beta.idx.null = 1:n.coeffs.null
  q.idx.null = n.coeffs.null + (1:nr.g.levels)
  beta.idx.altern = 1:n.coeffs.altern
  q.idx.altern = n.coeffs.altern + (1:nr.g.levels)
  
  obj.core = function (beta, q) {
    # The CC-only part    
    eta = c(ax %*% beta)
    exp.eta = exp(eta)
    P1 = exp.eta / (1 + exp.eta)
    P0 = 1 - P1
    p1 = exp.eta / (1 + exp.eta) ^ 2
    p0 = -p1
    
    P = rep(1, n)
    P[y0.idx] = P0[y0.idx]
    P[y1.idx] = P1[y1.idx]
    
    p = rep(0, n)
    p[y0.idx] = p0[y0.idx]
    p[y1.idx] = p1[y1.idx]
    
    ll.cc = sum(log(P))
    dll.cc = c((p / P) %*% ax, rep(0, nr.g.levels))
    
    # The CCpop G part
    ll.ccpop.g = sum(g.counts * log(q))
    dll.ccpop.g = c(rep(0, n.coeffs.altern), g.counts / q)
    
    # The CCpop E part
    eta = c(ax.expand %*% beta)
    exp.eta = exp(eta)
    P1 = exp.eta / (1 + exp.eta)
    P0 = 1 - P1
    p1 = exp.eta / (1 + exp.eta) ^ 2
    p0 = -p1
    Pw = w0 * P0 + w1 * P1
    pw = w0 * p0 + w1 * p1
    
    Q = q[g.idxs]
    
    c0mQPw = c0 + rowSums(matrix(Q * Pw, ncol = nr.g.levels))
    Qpw = (Q * pw) * ax.expand
    if (nr.g.levels == 2) {
      qPw = cbind(Pw[1:n], Pw[n + (1:n)])
      Qpw = Qpw[1:n, ] + Qpw[n + (1:n), ]
    } else {
      qPw = cbind(Pw[1:n], Pw[n + (1:n)], Pw[2*n + (1:n)])
      Qpw = Qpw[1:n, ] + Qpw[n + (1:n), ] + Qpw[2*n + (1:n), ]
    }
    
    ll.ccpop.e = sum(log(c0mQPw))
    dll.ccpop.e = colSums(cbind(Qpw, qPw) / c0mQPw)
    
    # Sum up
    obj = ll.ccpop.e - ll.cc - ll.ccpop.g
    grad = dll.ccpop.e - dll.cc - dll.ccpop.g
    
    return (list(objective = obj, gradient = grad))
  }
  
  obj.null = function (pars) {
    beta = rep(0, n.coeffs.altern)
    beta[-null.idxs.dropped] = pars[beta.idx.null]
    q = pars[q.idx.null]
    
    #res = obj.core(beta, q)
    res = snw_cnc_core_obj_core(beta, q, ax, ax.expand, y, g.counts, K)
    res$gradient = c(res$gradient[-null.idxs.dropped])

    return (res)
  }
  
  obj.altern = function (pars) {
    beta = pars[beta.idx.altern]
    q = pars[q.idx.altern]

    #return (obj.core(beta, q))
    return (snw_cnc_core_obj_core(beta, q, ax, ax.expand, y, g.counts, K))
  }
  
  con.null = function (pars) {
    return (list(constraints = sum(pars[q.idx.null]) - 1,
                 jacobian = rep(as.double(0:1), c(n.coeffs.null, nr.g.levels))))
  }

  con.altern = function (pars) {
    return (list(constraints = sum(pars[q.idx.altern]) - 1,
                 jacobian = rep(as.double(0:1), c(n.coeffs.altern, nr.g.levels))))
  }
  
  lb.null   = c(rep(-Inf, n.coeffs.null  ), rep(0, nr.g.levels))
  ub.null   = c(rep( Inf, n.coeffs.null  ), rep(1, nr.g.levels))
  lb.altern = c(rep(-Inf, n.coeffs.altern), rep(0, nr.g.levels))
  ub.altern = c(rep( Inf, n.coeffs.altern), rep(1, nr.g.levels))
  
  # Use impute-as-controls logistic regression to obtain initial param values
  y.imp = y
  y.imp[ym.idx] = 0
  pars.init.null = c(coef(glm.fit(ax[, -null.idxs.dropped], y.imp, family = binomial())), g.counts / sum(g.counts))
  pars.init.altern = c(coef(glm.fit(ax, y.imp, family = binomial())), g.counts / sum(g.counts))

  # Fit under null and alternative
  sol.null   = nloptr(pars.init.null  , eval_f = obj.null  , eval_g_ineq = con.null  , lb = lb.null  , ub = ub.null  , opts = list(algorithm = 'NLOPT_LD_MMA', xtol_rel = 1e-3, maxeval = 1000))
  sol.altern = nloptr(pars.init.altern, eval_f = obj.altern, eval_g_ineq = con.altern, lb = lb.altern, ub = ub.altern, opts = list(algorithm = 'NLOPT_LD_MMA', xtol_rel = 1e-3, maxeval = 1000))

  # FIXME check convergence
  
  null.ll = -sol.null$objective
  altern.ll = -sol.altern$objective
  attr(null.ll, 'df') = length(pars.init.null)
  attr(altern.ll, 'df') = length(pars.init.altern)
  
  return (list(null.ll = null.ll, altern.ll = altern.ll))
}

# The retrospective logistic model that exploits known prevalence and covariate independence.
# It's essentially Chatterjee and Carroll (2005), but I've implemented it differently.

cnc.glrt = function(formula.null, formula.altern, g.name, dat, K) {
  mf = model.frame(formula.altern, dat)
  y = model.extract(mf, 'response')
  ax = model.matrix(mf, dat)
  
  tg = as.data.frame(table(dat[, g.name]))
  nr.g.levels = nrow(tg)
  dat.expand = NULL
  for (i in 1:nr.g.levels) {
    dat.e = dat
    dat.e[, g.name] = i - 1
    dat.expand = rbind(dat.expand, dat.e)
  }
  rm(dat.e)
  mf.expand = model.frame(formula.altern, dat.expand)
  ax.expand = model.matrix(mf.expand, dat.expand)
    
  # Figure out the additional covariates added when moving from null to altern formula
  # FIXME there must be a better way to do this
  mf.null = model.frame(formula.null, dat[1, ])
  ax.null = model.matrix(mf.null, dat[1, ])
  null.idxs.dropped = which(!is.element(colnames(ax), colnames(ax.null)))

  res = snw.cnc.core(ax, y, ax.expand, tg$Freq, null.idxs.dropped, K)

  statistic = 2 * (res$altern.ll - res$null.ll)
  p.value = pchisq(statistic, df = length(null.idxs.dropped), lower.tail = F)
  
  return (list(statistic = statistic, p.value = p.value))
}

# The retrospective logistic model that exploits known prevalence, covariate independence and 
# population samples. Follows both Scott and Whild (2001) and Chatterjee and Carroll (2005)

snw.cnc.glrt = function(formula.null, formula.altern, g.name, dat, K) {
  mf = model.frame(formula.altern, dat, na.action = NULL)
  y = model.extract(mf, 'response')

  tg = as.data.frame(table(dat[, g.name]))
  nr.g.levels = nrow(tg)
  dat.expand = NULL
  for (i in 1:nr.g.levels) {
    dat.e = dat
    dat.e[, g.name] = i - 1
    dat.expand = rbind(dat.expand, dat.e)
  }
  rm(dat.e)
  mf.expand = model.frame(formula.altern, dat.expand, na.action = NULL)
  
  current.na.action = options('na.action')
  options(na.action = 'na.pass')
  ax = model.matrix(mf, dat)
  ax.expand = model.matrix(mf.expand, dat.expand)
  options(na.action = current.na.action$na.action)
  
  # Figure out the additional covariates added when moving from null to altern formula
  # FIXME there must be a better way to do this
  mf.null = model.frame(formula.null, dat[1, ])
  ax.null = model.matrix(mf.null, dat[1, ])
  null.idxs.dropped = which(!is.element(colnames(ax), colnames(ax.null)))
  
  res = snw.cnc.core(ax, y, ax.expand, tg$Freq, null.idxs.dropped, K)
  
  statistic = 2 * (res$altern.ll - res$null.ll)
  p.value = pchisq(statistic, df = length(null.idxs.dropped), lower.tail = F)
  
  return (list(statistic = statistic, p.value = p.value))
}

# A retrospective logistic model with known prevalence and population samples, where the conditional
# distribution of G|E is Balding-Nichols with known Fst.

bn.logistic.glrt = function(formula.null, formula.altern, g.name, a.name, dat, K, Fst) {
  # Extract the models for y (under the alternative) and g
  mf = model.frame(formula.altern, dat, na.action = NULL)
  y = model.extract(mf, 'response')

  # Construct an expanded data matrix of all combinations of values a of with those of g
  tg = as.data.frame(table(dat[, g.name]))
  nr.g.levels = nrow(tg)
  dat.expand = NULL
  for (i in 1:nr.g.levels) {
    dat.e = dat
    dat.e[, g.name] = i - 1
    dat.expand = rbind(dat.expand, dat.e)
  }
  rm(dat.e)
  mf.expand = model.frame(formula.altern, dat.expand, na.action = NULL)

  # Extract the actual data matrices
  ax = model.matrix(mf, dat) # this is only for the CC part
  current.na.action = options('na.action')
  options(na.action = 'na.pass')
  ax.expand = model.matrix(mf.expand, dat.expand)
  options(na.action = current.na.action$na.action)
  
  # Figure out the additional covariates added when moving from null to altern formula
  # FIXME there must be a better way to do this
  mf.null = model.frame(formula.null, dat[1, ])
  ax.null = model.matrix(mf.null, dat[1, ])
  null.idxs.dropped = which(!is.element(colnames(ax), colnames(ax.null)))
  
  # Prepare for optimization
  ym.idx = which(is.na(y))
  if (length(ym.idx) > 0) {
    y = y[-ym.idx]
  }
  
  g = dat[, g.name]
  a = dat[, a.name]
  g.expand = dat.expand[, g.name]
  a.expand = dat.expand[, a.name]
  
  y0.idx = which(y == 0)
  y1.idx = which(y == 1)
  
  g0.idx = which(g == 0)
  g1.idx = which(g == 1)
  g2.idx = which(g == 2)
  
  g.expand0.idx = which(g.expand == 0)
  g.expand1.idx = which(g.expand == 1)
  g.expand2.idx = which(g.expand == 2)
  
  n = length(y)
  n0 = length(y0.idx)
  n1 = length(y1.idx)
  c0 = length(ym.idx)
  N = n + c0
  N.expand = nrow(ax.expand)
  
  n.coeffs.altern = ncol(ax)
  n.coeffs.null = n.coeffs.altern - length(null.idxs.dropped)
  
  w0 = n0 / (1 - K)
  w1 = n1 / K
  
  beta.idx.null = 1:n.coeffs.null
  f0.idx.null = n.coeffs.null + 1
  beta.idx.altern = 1:n.coeffs.altern
  f0.idx.altern = n.coeffs.altern + 1

  obj.core = function (beta, f0) {
    # The CC-only part    
    eta = c(ax %*% beta)
    exp.eta = exp(eta)
    P = exp.eta / (1 + exp.eta)
    P[y0.idx] = 1 - P[y0.idx]
    p = exp.eta / (1 + exp.eta) ^ 2
    p[y0.idx] = -p[y0.idx]
    
    ll.cc = sum(log(P))
    dll.cc = c((p / P) %*% ax, 0)
    
    # The CCpop G|E part (FIXME computation here can probably be sped up considerably!)
    tmp.V = Fst * f0 * (1 - f0)
    d.tmp.V = Fst * (1 - 2 * f0)
    q = Q = rep(0, N)
    Q[g0.idx] = 1 + f0 * (f0 - 2) + tmp.V + 2 * a[g0.idx] * tmp.V * (a[g0.idx] - 1)
    Q[g1.idx] = 2 * (f0 - tmp.V - f0 ^ 2 + 2 * a[g1.idx] * tmp.V * (1 - a[g1.idx]))
    Q[g2.idx] = tmp.V + f0 ^ 2 - 2 * a[g2.idx] * tmp.V * (1 - a[g2.idx])
    q[g0.idx] = 2 * f0 - 2 + (1 - 2 * a[g0.idx] + 2 * a[g0.idx] ^ 2) * d.tmp.V
    q[g1.idx] = 2 * (1 - 2 * f0 + (-1 + 2 * a[g1.idx] - 2 * a[g1.idx] ^ 2) * d.tmp.V)
    q[g2.idx] = 2 * f0 + (1 - 2 * a[g2.idx] + 2 * a[g2.idx] ^ 2) * d.tmp.V
    
    ll.ccpop.ge = sum(log(Q))
    dll.ccpop.ge = c(rep(0, n.coeffs.altern), sum(q / Q))
    
    # The CCpop E part
    eta = c(ax.expand %*% beta)
    exp.eta = exp(eta)
    P1 = exp.eta / (1 + exp.eta)
    P0 = 1 - P1
    p1 = exp.eta / (1 + exp.eta) ^ 2
    p0 = -p1
    Pw = w0 * P0 + w1 * P1
    pw = w0 * p0 + w1 * p1
    
    q = Q = rep(0, N.expand)
    Q[g.expand0.idx] = 1 + f0 * (f0 - 2) + tmp.V + 2 * a.expand[g.expand0.idx] * tmp.V * (a.expand[g.expand0.idx] - 1)
    Q[g.expand1.idx] = 2 * (f0 - tmp.V - f0 ^ 2 + 2 * a.expand[g.expand1.idx] * tmp.V * (1 - a.expand[g.expand1.idx]))
    Q[g.expand2.idx] = tmp.V + f0 ^ 2 - 2 * a.expand[g.expand2.idx] * tmp.V * (1 - a.expand[g.expand2.idx])
    q[g.expand0.idx] = 2 * f0 - 2 + (1 - 2 * a.expand[g.expand0.idx] + 2 * a.expand[g.expand0.idx] ^ 2) * d.tmp.V
    q[g.expand1.idx] = 2 * (1 - 2 * f0 + (-1 + 2 * a.expand[g.expand1.idx] - 2 * a.expand[g.expand1.idx] ^ 2) * d.tmp.V)
    q[g.expand2.idx] = 2 * f0 + (1 - 2 * a.expand[g.expand2.idx] + 2 * a.expand[g.expand2.idx] ^ 2) * d.tmp.V
    
    QPw = rowSums(matrix(Q * Pw, nrow = N))
    qPw = rowSums(matrix((q * Pw), nrow = N))
    Qpw = (Q * pw) * ax.expand
    Qpw = Qpw[g.expand0.idx, ] + Qpw[g.expand1.idx, ] + Qpw[g.expand2.idx, ]
    c0mQPw = c0 + QPw
    
    ll.ccpop.e = sum(log(c0mQPw))
    dll.ccpop.e = colSums(cbind(Qpw, qPw) / c0mQPw)
    
    # Sum up
    obj = ll.ccpop.e - ll.cc - ll.ccpop.ge
    grad = dll.ccpop.e - dll.cc - dll.ccpop.ge
    
    return (list(objective = obj, gradient = grad))
  }
  
  obj.null = function (pars) {
    beta = rep(0, n.coeffs.altern)
    beta[-null.idxs.dropped] = pars[beta.idx.null]
    f0 = pars[f0.idx.null]
    
    #res = obj.core(beta, f0)
    res = bn_logistic_glrt_obj_core(beta, f0, y, g, a, ax, ax.expand, g.expand, a.expand, K, Fst, N, n0, n1, c0)
    res$gradient = c(res$gradient)[-null.idxs.dropped]
    
    return (res)
  }
  
  obj.altern = function (pars) {
    beta = pars[beta.idx.altern]
    f0 = pars[f0.idx.altern]
    
    #res = obj.core(beta, f0)
    res = bn_logistic_glrt_obj_core(beta, f0, y, g, a, ax, ax.expand, g.expand, a.expand, K, Fst, N, n0, n1, c0)
    res$gradient = c(res$gradient)
    
    return (res)
  }
  
  # Use impute-as-controls estimators to obtain initial param values
  y.imp = y
  y.imp[ym.idx] = 0
  f0.init = min((length(g2.idx) * 2 + length(g1.idx)) / (2 * N), 0.45)
  pars.init.null = c(coef(glm.fit(ax[, -null.idxs.dropped], y.imp, family = binomial())), f0.init)
  pars.init.altern = c(coef(glm.fit(ax, y.imp, family = binomial())), f0.init)
  
  lb.null   = c(rep(-Inf, n.coeffs.null  ), 0  )
  lb.altern = c(rep(-Inf, n.coeffs.altern), 0  )
  ub.null   = c(rep( Inf, n.coeffs.null  ), 0.5)
  ub.altern = c(rep( Inf, n.coeffs.altern), 0.5)
  
  # Fit under null and alternative
  sol.null   = nloptr(pars.init.null  , eval_f = obj.null  , lb = lb.null  , ub = ub.null  , opts = list(algorithm = 'NLOPT_LD_MMA', xtol_rel = 1e-3, maxeval = 1000))
  sol.altern = nloptr(pars.init.altern, eval_f = obj.altern, lb = lb.altern, ub = ub.altern, opts = list(algorithm = 'NLOPT_LD_MMA', xtol_rel = 1e-3, maxeval = 1000))
  
  # FIXME check convergence
  
  null.ll = -sol.null$objective
  altern.ll = -sol.altern$objective
  
  statistic = 2 * (altern.ll - null.ll)
  p.value = pchisq(statistic, df = length(null.idxs.dropped), lower.tail = F)
  
  return (list(statistic = statistic, p.value = p.value))
}

# A retrospective logistic model with known prevalence and population samples, where the conditional
# distribution of G|E is also logistic.

double.logistic.glrt = function(formula.null, formula.altern, formula.g, dat, K) {
  # FIXME there is a lot of data replication here. Optimize?
  
  # Extract the models for y (under the alternative) and g
  mf = model.frame(formula.altern, dat, na.action = NULL)
  mf.g = model.frame(formula.g, dat, na.action = NULL)
  y = model.extract(mf, 'response')
  g = model.extract(mf.g, 'response')
  
  # Construct an expanded data matrix of all combinations of values of e (maybe vector valued) with
  # those of g
  g.name = all.vars(formula.g)[1]
  tg = as.data.frame(table(g))
  nr.g.levels = nrow(tg)
  nr.g.alleles = ifelse(nr.g.levels == 3, 2, 1) # FIXME do we want to just always assume 2 alleles in HWE?
  dat.expand = NULL
  for (i in 1:nr.g.levels) {
    dat.e = dat
    dat.e[, g.name] = as.numeric(tg$g[i]) - 1
    dat.expand = rbind(dat.expand, dat.e)
  }
  rm(dat.e)
  mf.expand = model.frame(formula.altern, dat.expand, na.action = NULL)
  mf.g.expand = model.frame(formula.g, dat.expand, na.action = NULL)
  g.expand = model.extract(mf.g.expand, 'response')
  
  # Extract the actual data matrices
  current.na.action = options('na.action')
  options(na.action = 'na.omit')
  ax = model.matrix(mf, dat) # this is only for the CC part
  options(na.action = 'na.pass')
  ax.g = model.matrix(mf.g, dat)
  ax.expand = model.matrix(mf.expand, dat.expand)
  ax.g.expand = model.matrix(mf.g.expand, dat.expand)
  options(na.action = current.na.action$na.action)
  
  # Figure out the additional covariates added when moving from null to altern formula
  # FIXME there must be a better way to do this
  mf.null = model.frame(formula.null, dat[1, ])
  ax.null = model.matrix(mf.null, dat[1, ])
  null.idxs.dropped = which(!is.element(colnames(ax), colnames(ax.null)))

  # Prepare for optimization
  ym.idx = which(is.na(y))
  if (length(ym.idx) > 0) {
    y = y[-ym.idx]
  }
  
  y0.idx = which(y == 0)
  y1.idx = which(y == 1)

  g0.idx = which(g == 0)
  g1.idx = which(g == 1)
  g2.idx = which(g == 2)

  g.expand0.idx = which(g.expand == 0)
  g.expand1.idx = which(g.expand == 1)
  g.expand2.idx = which(g.expand == 2)
  
  n = length(y)
  n0 = length(y0.idx)
  n1 = length(y1.idx)
  c0 = length(ym.idx)
  N = n + c0
  
  n.coeffs.altern = ncol(ax)
  n.coeffs.null = n.coeffs.altern - length(null.idxs.dropped)
  n.coeffs.g = ncol(ax.g)
  
  w0 = n0 / (1 - K)
  w1 = n1 / K
  
  beta.idx.null = 1:n.coeffs.null
  theta.idx.null = n.coeffs.null + (1:n.coeffs.g)
  beta.idx.altern = 1:n.coeffs.altern
  theta.idx.altern = n.coeffs.altern + (1:n.coeffs.g)
  
  obj.core = function (beta, theta) {
    # The CC-only part    
    eta = c(ax %*% beta)
    exp.eta = exp(eta)
    P = exp.eta / (1 + exp.eta)
    P[y0.idx] = 1 - P[y0.idx]
    p = exp.eta / (1 + exp.eta) ^ 2
    p[y0.idx] = -p[y0.idx]
        
    ll.cc = sum(log(P))
    dll.cc = c((p / P) %*% ax, rep(0, n.coeffs.g))
    
    # The CCpop G|E part
    eta.g = c(ax.g %*% theta)
    exp.eta.g = exp(eta.g)
    maf.g = exp.eta.g / (1 + exp.eta.g)
    Q = dbinom(g, nr.g.alleles, maf.g)
    if (nr.g.alleles == 2) {
      q = 2 * exp.eta.g / (1 + exp.eta.g) ^ 2
      q[g0.idx] = q[g0.idx] * (maf.g[g0.idx] - 1)
      q[g1.idx] = q[g1.idx] * (1 - 2 * maf.g[g1.idx])
      q[g2.idx] = q[g2.idx] * maf.g[g2.idx]
    } else {
      q = exp.eta.g / (1 + exp.eta.g) ^ 2
      q[g0.idx] = -q[g0.idx]
    }
    
    ll.ccpop.ge = sum(log(Q))
    dll.ccpop.ge = c(rep(0, n.coeffs.altern), (q / Q) %*% ax.g)
    
    # The CCpop E part
    eta = c(ax.expand %*% beta)
    exp.eta = exp(eta)
    P1 = exp.eta / (1 + exp.eta)
    P0 = 1 - P1
    p1 = exp.eta / (1 + exp.eta) ^ 2
    p0 = -p1
    Pw = w0 * P0 + w1 * P1
    pw = w0 * p0 + w1 * p1
    
    eta.g = c(ax.g.expand %*% theta)
    exp.eta.g = exp(eta.g)
    maf.g = exp.eta.g / (1 + exp.eta.g)
    Q = dbinom(g.expand, nr.g.alleles, maf.g)
    if (nr.g.alleles == 2) {
      q = 2 * exp.eta.g / (1 + exp.eta.g) ^ 2
      q[g.expand0.idx] = q[g.expand0.idx] * (maf.g[g.expand0.idx] - 1)
      q[g.expand1.idx] = q[g.expand1.idx] * (1 - 2 * maf.g[g.expand1.idx])
      q[g.expand2.idx] = q[g.expand2.idx] * maf.g[g.expand2.idx]
    } else {
      q = exp.eta.g / (1 + exp.eta.g) ^ 2
      q[g.expand0.idx] = -q[g.expand0.idx]
    }
    
    QPw = rowSums(matrix(Q * Pw, nrow = N))
    qPw = apply(array((q * Pw) * ax.g.expand, dim = c(N, nr.g.levels, n.coeffs.g)), c(1, 3), sum)
    Qpw = apply(array((Q * pw) * ax.expand, dim = c(N, nr.g.levels, n.coeffs.altern)), c(1, 3), sum)
    c0mQPw = c0 + QPw
  
    ll.ccpop.e = sum(log(c0mQPw))
    dll.ccpop.e = colSums(cbind(Qpw, qPw) / c0mQPw)
    
    # Sum up
    obj = ll.ccpop.e - ll.cc - ll.ccpop.ge
    grad = dll.ccpop.e - dll.cc - dll.ccpop.ge
    
    return (list(objective = obj, gradient = grad))
  }
  
  obj.null = function (pars) {
    beta = rep(0, n.coeffs.altern)
    beta[-null.idxs.dropped] = pars[beta.idx.null]
    theta = pars[theta.idx.null]
    
    #res = obj.core(beta, theta)
    res = double_logistic_glrt_obj_core(beta, theta, y, g, ax, ax.g, ax.expand, ax.g.expand, g.expand, w0, w1)
    res$gradient = c(res$gradient)[-null.idxs.dropped]

    #cat('NULL @')
    #cat(sprintf('% .02f ', beta))
    #cat(sprintf('% .02f ', theta))
    #cat('=> value =', res$objective, ' grad =', res$gradient, '\n')
    
    return (res)
  }
  
  obj.altern = function (pars) {
    beta = pars[beta.idx.altern]
    theta = pars[theta.idx.altern]

    #res = obj.core(beta, theta)
    res = double_logistic_glrt_obj_core(beta, theta, y, g, ax, ax.g, ax.expand, ax.g.expand, g.expand, w0, w1)
    res$gradient = c(res$gradient)
    
    #cat('ALTN @')
    #cat(sprintf('% .02f ', beta))
    #cat(sprintf('% .02f ', theta))
    #cat('=> value =', res$objective, ' grad =', res$gradient, '\n')
    
    return (res)
  }
  
  # Use impute-as-controls logistic regression to obtain initial param values
  dat[, g.name] = factor(g, levels = 0:2)
  theta.init = coef(glm(formula.g, dat, family = binomial()))
  dat[, g.name] = g
  pars.init.null = c(coef.imputed.logistic(formula.null, dat), theta.init)
  #pars.init.altern = c(coef.imputed.logistic(formula.altern, dat), theta.init)
  
  # Fit under null and alternative
  sol.null   = nloptr(pars.init.null  , eval_f = obj.null  , opts = list(algorithm = 'NLOPT_LD_MMA', xtol_rel = 1e-3, maxeval = 1000))
  
  # Use the null solution as the starting point for finding the MLE under the alternative
  beta.init = rep(0, n.coeffs.altern)
  beta.init[-null.idxs.dropped] = sol.null$solution[beta.idx.null]
  pars.init.altern = c(beta.init, sol.null$solution[theta.idx.null])
  
  sol.altern = nloptr(pars.init.altern, eval_f = obj.altern, opts = list(algorithm = 'NLOPT_LD_MMA', xtol_rel = 1e-3, maxeval = 1000))
    
  # FIXME check convergence
  
  null.ll = -sol.null$objective
  altern.ll = -sol.altern$objective
    
  statistic = 2 * (altern.ll - null.ll)
  p.value = pchisq(statistic, df = length(null.idxs.dropped), lower.tail = F)
  
  return (list(statistic = statistic, p.value = p.value))
}

coef.imputed.logistic = function(formula, dat) {
  mf = model.frame(formula, dat, na.action = NULL)
  y = model.extract(mf, 'response')
  y[is.na(y)] = 0
  current.na.action = options('na.action')
  options(na.action = 'na.pass')
  ax = model.matrix(mf, dat)
  options(na.action = current.na.action$na.action)
  model = glm.fit(ax, y, family = binomial())
  return (coef(model))
}
