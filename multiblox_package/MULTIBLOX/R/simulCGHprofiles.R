simulCGHprofiles <-
function (M, n, k.mean, SNR, lambda, seed) 
{
  ### M number of profiles
  ### n number of probes (variables) for each profile
  ### k.mean average number of aberrant segment for each profile
  ### SNR Signal/Noise Ratio
  ### lambda Variance ratio between measurement noise and bakground noise
  ### seed
  
  set.seed(seed)
  P = c(0.5, 0, 0.5)
  Pcum = cumsum(P)
  K = 1 + rpois(M, k.mean - 1) # number of segments for each profile
  k.max = max(K)
  X = matrix(nrow = M, ncol = n)
  E = matrix(nrow = M, ncol = n)
  mu = matrix(nrow = M, ncol = n)
  tau = matrix(NaN, nrow = M, ncol = k.max) # levels of each segment over the M profiles
  for (m in (1:M)) { # for each profile
    tau.tmp = ceiling((n - 1) * runif(K[m] - 1))
    tau.tmp = sort(tau.tmp)
    tau[m, (1:K[m])] = c(tau.tmp, n)
    level = vector(length = K[m])
    level[1] = 0
    mu.tmp = rep(level[1], tau[m, 1])
    if (K[m] > 1) {
      for (k in (2:K[m])) {
        if (mu.tmp[tau[m, k - 1]] != 0) 
          level[k] = 0
        else level[k] = -1 + sum(runif(1) > Pcum)
        mu.tmp = c(mu.tmp, rep(level[k], (tau[m, k] - 
                                            tau[m, k - 1])))
      }
    }
    mu[m, ] = mu.tmp
  }
  mu[mu == 1] = 0.58
  sigma0 = sqrt(sum(mu^2)/(n * M * SNR))
  TrueRef = t(matrix(rnorm(n, 0, sigma0/lambda) + (sigma0/lambda) * 
                       sin(2 * pi * (1:n)/n * 10), ncol = M, nrow = n))
  for (m in (1:M)) {
    E[m, ] = rnorm(n, 0, sigma0)
  }
  X = mu + TrueRef + E
  Y = data.frame(t(X))
  colnames(Y) = paste(rep("Ind", M), c(1:M), sep = "")
  invisible(list(Y = Y, Km = K, mu0 = mu, sigma0 = sigma0, 
                 theta0 = TrueRef[1, ]))
}
