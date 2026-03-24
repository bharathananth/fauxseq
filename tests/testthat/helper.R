
is.wholenumber <- function(x, tol = 2 * .Machine$double.eps) {
  abs(x - round(x)) < tol
}
