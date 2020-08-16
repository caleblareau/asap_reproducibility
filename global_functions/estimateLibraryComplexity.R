
# Function translated from java version: https://github.com/broadinstitute/picard/blob/master/src/main/java/picard/sam/DuplicationMetrics.java
# Not vectorized!!! 
estimateLibrarySize <- function(nTotal, nUnique){
  
  f <- function(x, c, n) {
    return(c / x - 1 + exp(-n / x))
  }
  
  m = 1
  M = 100
  
  nDuplicates <- (nTotal - nUnique) + 1 # charity to handle only unique reads observed
  
  if (nUnique > nTotal | (f(m * nUnique, nUnique, nTotal) < 0) | nUnique < 0 | nTotal < 0 | nDuplicates < 0) {
    message("Library size returns 0 -- invalid inputs; check this cell more closely")
    return(0)
  }
  
  while (f(M * nUnique, nUnique, nTotal) > 0) {
    M <- M*10.0
  }
  
  for(i in seq(0, 40)) {
    r <-  (m + M) / 2.0
    u <- f(r * nUnique, nUnique, nTotal);
    if (u == 0) {
      break
    } else if (u > 0) {
      m = r
    } else if (u < 0) {
      M = r
    }
  }
  
  return(round(nUnique * (m + M) / 2.0))
}