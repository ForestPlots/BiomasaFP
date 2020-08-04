# Helper functions to calculate instantanious productivity, mortality and stem turnover

# Functions adapted from Kohyama et al J Ecology (https://github.com/kohyamat/p-B)

#' @author T.S. Kohyama et al.

productivity <- function(dbh1, dbh2, w1, w2, t, area, Alive,Recruit) {
  si <- ifelse(Alive==1 & Recruit==0,1,0) # survival
  di <- ifelse(Alive==0,1,0) # death
  ri <- Recruit # recruitment
  Ns0 <- sum(si, na.rm = TRUE)
  N0 <- Ns0 + sum(di, na.rm = TRUE)
  NT <- Ns0 + sum(ri, na.rm = TRUE)
  Bs0 <- sum(si * w1, na.rm = TRUE)
  BsT <- sum(si * w2, na.rm = TRUE)
  B0 <- Bs0 + sum(di * w1, na.rm = TRUE)
  BT <- BsT + sum(ri * w2, na.rm = TRUE)
  # period-mean biomass and abundance
  Nw <- ifelse(NT != N0, (NT - N0) / log(NT / N0), N0)
  N <- Nw / area # (per ha)
  Bw <- ifelse(BT != B0, (BT - B0) / log(BT / B0), B0)
  B <- Bw / area

  # Standardized maximum tree mass for initial population
  W_max <- as.numeric(quantile(w1[ri != 1], 0.99)) # Mg

  # turnover rates
  r <- try(turnover_est(si + ri, si, t),silent=TRUE)
  if(inherits(r,"try-error")){
    r<-NA
  }
  m <- try(turnover_est(si + di, si, t),silent=TRUE)
  if(inherits(m,"try-error")){
    m<-NA
  }
  p <- try(turnover_est(w2, si * w1, t),silent=TRUE)
  if(inherits(p,"try-error")){
    p<-NA
  }
  l <- try(turnover_est(w1, si * w1, t),silent=TRUE)
  if(inherits(l,"try-error")){
    l<-NA
  }

  # absolute productivity (Mg per ha per year)
  P <- p * B
  P_simple <- sum(((si + ri) * w2 - si * w1) / t)
  P_simple <- P_simple / area

  return(c(
    "B" = B,
    "N" = N,
    "W_max" = W_max,
    "p" = p,
    "l" = l,
    "r" = r,
    "m" = m,
    "P" = P,
    "P_simple" = P_simple
  ))
}

turnover_est <- function(y, z, t) {
  f <- function(rho) {
    sum(y * exp(-rho * t) - z)
  }
  df <- function(rho) {
    sum(-t * y * exp(-rho * t))
  }
  # Newton-Rapton iteration
  rho <- 0.02
  precision <- 1.0e-12 # to stop iteration
  change <- precision + 1.0

  while (change > precision) {
    rho2 <- rho - f(rho) / df(rho)
    change <- abs(rho2 - rho)
    rho <- rho2
  }
  return(rho)
}
