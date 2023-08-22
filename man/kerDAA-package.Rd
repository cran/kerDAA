\name{kerDAA}
\alias{kerDAA-package}
\title{New kernel-based test for differential association analysis}
\description{This package can be used to determine whether two high-dimensional samples have similar dependence relationships across two conditions.}
\author{
Hoseung Song and Michael C. Wu

Maintainer: Hoseung Song (hsong3@fredhutch.org)
}
\references{
Song, H. and Wu, M.C. (2023). Multivariate differential association analysis. arXiv:2307.15268
}

\seealso{
  \code{\link{kerdaa}}
}
\examples{

# Dimension of variables.
d = 100

# The first covariance matrix
SIG = matrix(0, d, d)
for (i in 1:d) {
  for (j in 1:d) {
    SIG[i,j] = 0.4^(abs(i-j))
  }
}

# The second covariance matrix
SIG1 = matrix(0, d, d)
for (i in 1:d) {
  for (j in 1:d) {
    SIG1[i,j] = (0.4+0.5)^(abs(i-j))
  }
}

set.seed(500)

# We use 'rmvnorm' in 'mvtnorm' package to generate multivariate normally distributed samples
require(mvtnorm)
Z = rmvnorm(100, mean = rep(0,100), sigma = SIG)
X1 = Z[,1:50]
Y1 = Z[,51:100]

Z = rmvnorm(100, mean = rep(0,100), sigma = SIG1)
X2 = Z[,1:50]
Y2 = Z[,51:100]

a = kerdaa(X1, Y1, X2, Y2, perm=1000)
# output results based on the permutation and the asymptotic results
# the test statistic values can be found in a$stat_g and a$stat_l
# p-values can be found in a$pval and a$pval_perm

}
