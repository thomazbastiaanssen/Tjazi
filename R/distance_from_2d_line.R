distance_from_2d_line <- function(a,b,c) {
  #a = point, b and c are two points on your line of interest.
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
}