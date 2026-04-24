# ============================================================
# TEST DE EXPONENCIALIDAD DE CABANA & CABANA
# ============================================================

# --- Funciones auxiliares ---------------------------------------------------

# Polinomios de Laguerre L_0..L_l evaluados en x: devuelve matriz (l+1) x n
laguerre_matrix <- function(x, l) {
  n <- length(x)
  L <- matrix(0, nrow = l + 1, ncol = n)
  L[1, ] <- 1
  if (l == 0) return(L)
  L[2, ] <- 1 - x
  for (k in seq_len(l - 1)) {
    L[k + 2, ] <- ((2*k + 1 - x) * L[k + 1, ] - k * L[k, ]) / (k + 1)
  }
  L
}


TEPest <- function(X){
	  n <- length(X)
	Y <- X / mean(X)

  # Estadistico
  Phi <- laguerre_matrix(Y, 12)
  lb  <- rowSums(Phi[3:13, , drop=FALSE]) / sqrt(n)
  Qn  <- as.numeric(t(lb) %*% Cdat %*% lb)	
  return(Qn)
}

#'@title Exponentiality test based on TEP's
#'
#' @description Performs a test of exponentiality consistent
#' and focused to detect Weibul alternatives
#'
#' @details The tests statistics are quadratic forms in the vector of evaluations
#' of the normalized Laguerre polynomials at the scaled sample points,
#' introduced by the authors in
#' \emph{Goodness-of-fit to the Exponential Distribution, focused on Weibull alternatives}, 
#'  Communications in Statistics, Simulation and Computation \bold{34}, 711-723, 2005.
#' @param X Vector containing the elements of the sample
#' @param l Integer, Number of Laguerre polynomials evaluated at the sample points -1 (defalts to 10)
#' @param nrep Integer, Number of Monte Carlo evaluations to estimate the p-value
#' (defalts to 1000)
#'
#' @return A list with
#'   \itemize{
#'     \item `"pval"` — The p-value of the sample estimated by MC
#'     \item `"est"` — The value of the test statistic 
#'   }
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' X <- rexp(100)
#' TEPexp(X)
#' }
#'
#' @export
TEPexp <- function(X,nrep=1000){
	Qdat <- TEPest(X)
	Qrep <- sapply(1:nrep,function(i){Xsim=rexp(length(X))
		TEPest(Xsim)})
		pval <- mean(Qrep>rep(Qdat,length(Qrep)))
		return(list(pval=pval,est=Qdat))
}
