calc_nbinom_weighted <- Vectorize(function(q0, q1, new_mu, new_theta,
                                           p0, p0_new, p1, p1_new) {
  qvec <- q0:q1
  dvec <- c(0, dnbinom((q0+1):(q1-1), mu = new_mu, size = new_theta), 0)
  dvec[1] <- p0_new - p0
  dvec[length(dvec)] <- p1_new - p1
  dvec <- abs(dvec)
  sum(qvec*dvec)/sum(dvec)
})

get_nbinomx <- function(x, mu, theta, new_mu, new_theta) {
  out <- rep(0, length(x))
  p0 <- pnbinom(x-1, mu = mu, size = theta, log.p = TRUE)
  p1 <- pnbinom(x, mu = mu, size = theta, log.p = TRUE)
  q0 <- qnbinom(p0, mu = new_mu, size = new_theta, log.p = TRUE)
  q1 <- qnbinom(p1, mu = new_mu, size = new_theta, log.p = TRUE)

  ind1 <- which(q0 == q1 & q1 != Inf)
  ind2 <- which(q0 != q1 & q1 != Inf)
  ind3 <- which(q1 == Inf)

  if (length(ind1) > 0) {
    out[ind1] <- q1[ind1]
  }

  if (length(ind2) > 0) {
    p0_new <- pnbinom(q0[ind2], mu = new_mu[ind2], size = new_theta[ind2])

    q1away <- abs(q0[ind2] - q1[ind2]) == 1

    if (sum(q1away) > 0) {
      d1 <- abs(p0_new[q1away] - exp(p0[ind2][q1away]))
      d2 <- abs(exp(p1[ind2][q1away]) - p0_new[q1away])
      out[ind2[q1away]] <- (q0[ind2[q1away]]*d1 + q1[ind2[q1away]]*d2)/(d1+d2)
    }

    if (sum(!q1away) > 0) {
      p1_new <- pnbinom(q1[ind2[!q1away]]-1, mu = new_mu[ind2[!q1away]],
                        size = new_theta[ind2[!q1away]])
      out[ind2[!q1away]] <- calc_nbinom_weighted(
        q0[ind2[!q1away]], q1[ind2[!q1away]], new_mu[ind2[!q1away]],
        new_theta[ind2[!q1away]], exp(p0[ind2[!q1away]]), p0_new[!q1away],
        exp(p1[ind2[!q1away]]), p1_new)

    }
  }

  if (length(ind3) > 0) {
    x1 <- x[ind3]
    mu1 <- mu[ind3]
    theta1 <- theta[ind3]
    new_mu1 <- new_mu[ind3]
    new_theta1 <- new_theta[ind3]
    sdx <- sqrt(mu1 + mu1^2/theta1)
    sdnewx <- sqrt(new_mu1 + new_mu1^2/new_theta1)
    diff <- (x1 - mu1)/sdx*sdnewx
    lower <- (new_mu1 + diff)*1/2
    upper <- (new_mu1 + diff)*2
    for (i in seq_along(ind3)) {
      k <- ind3[i]
      cand <- floor(seq(lower[i], upper[i], length = 10))
      d1 <- dnbinom(x[k], mu = mu[k], size = theta[k], log = TRUE)
      d2 <- dnbinom(cand, mu = new_mu[k], size = new_theta[k], log = TRUE)
      out[k] <- cand[which.min(abs(d1 - d2))]
    }
  }
  round(out, 4)
}
