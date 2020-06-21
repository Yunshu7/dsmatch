# method of bisection to find a x such that f[x] = y for an increasing f
bisect <- function (f, lo, hi, y = 0, ytol = 1e-12, itmax = 100) {
  mi <- (lo + hi) / 2
  f_mi <- f(mi)
  i <- 0
  
  while ((i < itmax) & (abs(f_mi - y) > ytol)) {
    if (f_mi > 0)
      hi <- mi
    else
      lo <- mi
    mi <- (lo + hi) / 2
    f_mi <- f(mi)
    i <- i + 1
  }
  
  return(mi)
}