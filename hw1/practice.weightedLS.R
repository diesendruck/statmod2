ols.heterosked.example = function(n) {
  y = 3-2*x + rnorm(n,0,sapply(x,function(x){1+0.5*x^2}))
  fit.ols = lm(y~x)
  # Return the errors
  return(fit.ols$coefficients - c(3,-2))
}

ols.heterosked.error.stats = function(n,m=10000) {
  ols.errors.raw = t(replicate(m,ols.heterosked.example(n)))
  # transpose gives us a matrix with named columns
  intercept.sd = sd(ols.errors.raw[,"(Intercept)"])
  slope.sd = sd(ols.errors.raw[,"x"])
  return(list(intercept.sd=intercept.sd,slope.sd=slope.sd))
}

wls.heterosked.example = function(n) {
  y = 3-2*x + rnorm(n,0,sapply(x,function(x){1+0.5*x^2}))
  fit.wls = lm(y~x,weights=1/(1+0.5*x^2))
  # Return the errors
  return(fit.wls$coefficients - c(3,-2))
}

wls.heterosked.error.stats = function(n,m=10000) {
  wls.errors.raw = t(replicate(m,wls.heterosked.example(n)))
  # transpose gives us a matrix with named columns
  intercept.sd = sd(wls.errors.raw[,"(Intercept)"])
  slope.sd = sd(wls.errors.raw[,"x"])
  return(list(intercept.sd=intercept.sd,slope.sd=slope.sd))
}

n=100
x = rnorm(n,0,3)

ols.heterosked.example(n)
ols.heterosked.error.stats(n, m=10000)

wls.heterosked.example(n)
wls.heterosked.error.stats(n, m=10000)
