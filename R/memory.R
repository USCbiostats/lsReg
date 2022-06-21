# Base allocation
allocate.lsregmem <- function(mdl, mdltype, q, testtype) {
  data <- mdl$model
  attr(data, "terms") <- NULL
  n <- nrow(data)
  p <- ncol(data)
  q <- q
  if (mdltype == 1) {
    return (allocate.lslinregmem(data = data,
                                 beta0 = mdl$coefficients,
                                 n = n,
                                 p = p,
                                 q = q,
                                 testtype = testtype))
  } else if (mdltype == 2) {
    return (allocate.lslogregmem(data = data,
                                 beta0 = mdl$coefficients,
                                 n = n,
                                 p = p,
                                 q = q,
                                 testtype = testtype))
  }
  return (allocate.lspoisregmem(data = data,
                                beta0 = mdl$coefficients,
                                n = n,
                                p = p,
                                q = q,
                                testtype = testtype))
}
