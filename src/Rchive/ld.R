LD <- function (nAB, nAb, naB, nab) {
  n <- nAB + nAb + naB + nab
  nA <- nAB + nAb
  na <- naB + nab
  nB <- nAB + naB
  nb <- nAb + nab

  print(
    str_c(
      "nA = ", round(nA, 2),
      ", na = ", round(na, 2),
      ", nB = ", round(nB, 2),
      ", nb = ", round(nb, 2)
    )
  )
  
  pA <- nA / n
  pa <- 1 - pA
  pB <- nB / n
  pb <- 1 - pB

  print(
    str_c(
      "pA = ", round(pA, 2),
      ", pa = ", round(pa, 2),
      ", nB = ", round(pB, 2),
      ", pb = ", round(pb, 2)
    )
  )


  DA <- pA - pA^2
  DB <- pB - pB^2

  delta <- (((nAB + nab) - (naB + nAb)) / (2 * n)) -
    ((na - nA) * (nb - nB) / (2 * n^2))

  t <- (pA * pa + DA) * (pB * pb + DB)

  delta / sqrt(t)
}

LD(300, 0, 10, 50)