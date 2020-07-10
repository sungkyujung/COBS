 bisect.dec <- function(h, l, t1, t2, Tol = .Machine$double.eps) {
  # Rootfinding for strictly decreasing function h
  # with h(l) = infty, h(t1) > 0 > h(t2)
  # 
  # Search range (t1,t2)
  #  
  # If sign of h(t1) is negative, then extend the search range
  while(h(t1) < 0) {
    t1 <- t1 - (t1 - l)/2
  }
  # If sign of h(t2) is positive, then extend the search range
  while(h(t2) > 0) {
    t2 <- t2 + (t2 - t1)/2
  } 
  
  for (i in 1:100) {
    
    tmid <- (t1 + t2) / 2 # Calculate midpoint
    h.tmid <- h(tmid)
    
    
    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the 
    # function and return the root.
    if ( h.tmid == 0 || (t2 - tmid) < Tol ){
      return(tmid)
    }
    
    # Update the search range 
    ifelse( h.tmid > 0,  t1 <- tmid, t2 <- tmid) 
  } 
  return(tmid)
}
