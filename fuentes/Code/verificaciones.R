# ==============================================================================
# VERIFICACIONES BAYESIANAS (3.2, 3.3 y 3.4)
# ------------------------------------------------------------------------------
# Sobre las cadenas guardadas en Code/resultados/:
#   1. Log-verosimilitud MARGINAL por iteracion (via espectral) y DIC
#      (3.3 y 3.4; el Ising se excluye: verosimilitud doblemente intractable)
#   2. Probabilidades posteriores de contrastes (3.2, 3.3 y 3.4)
#   3. Prior predictive check: simulacion desde las previas (3.3 y 3.4)
#   4. Posterior predictive checks (ppp) con imagenes replicadas (3.3 y 3.4)
# La log-verosimilitud marginal de y ~ N(mu*1, Sigma_y) se evalua con la base
# de autovalores de la Laplaciana (validada contra la densa en una rejilla
# pequena al inicio; ver validacion interna).
# Uso: Rscript verificaciones.R
# ==============================================================================
suppressPackageStartupMessages({
  library(Matrix); library(mvtnorm); library(coda)
})

# ------------------------------------------------------------------------------
# 0. CONFIGURACION Y RUTAS
# ------------------------------------------------------------------------------
if (!exists("dir_out")) {
  archivo_config <- if (identical(.Platform$OS.type, "windows")) "config_windows.R" else "config_linux.R"
  dirs <- c(".", "Code", "..", file.path("..", "Code"))
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
  if (length(f) == 1L && !is.na(f)) dirs <- c(dirname(gsub("~+~", " ", f, fixed = TRUE)), dirs)
  for (i in seq_len(sys.nframe())) {
    g <- sys.frame(i)$ofile
    if (!is.null(g) && nzchar(g)) { dirs <- c(dirname(gsub("~+~", " ", g, fixed = TRUE)), dirs); break }
  }
  cargado <- FALSE
  for (d in unique(dirs)) {
    if (file.exists(file.path(d, archivo_config))) {
      source(file.path(d, archivo_config)); cargado <- TRUE; break
    }
  }
  if (!cargado) stop("No se encontro '", archivo_config, "'.")
}
dir_ver <- file.path(dir_out, "verificacion")
dir.create(dir_ver, showWarnings = FALSE, recursive = TRUE)
cat("Resultados de verificacion en:", dir_ver, "\n")

# ------------------------------------------------------------------------------
# 1. ESTRUCTURA ESPACIAL, BASE DE AUTOVECTORES Y LOG-VEROSIMILITUD MARGINAL
# ------------------------------------------------------------------------------
construir_reja <- function(M, N) {
  P <- M * N
  idx <- matrix(seq_len(P), M, N)
  i <- c(as.vector(idx[, -N]), as.vector(idx[-M, ]))
  j <- c(as.vector(idx[, -1]), as.vector(idx[-1, ]))
  A <- sparseMatrix(i = c(i, j), j = c(j, i), x = 1, dims = c(P, P))
  L <- Diagonal(x = rowSums(A)) - A
  w <- outer(4 * sin(pi * (0:(M - 1)) / (2 * M))^2,
             4 * sin(pi * (0:(N - 1)) / (2 * N))^2, "+")
  list(L = L, w = as.vector(w), aristas = cbind(i, j), P = P)
}

Cmat <- function(n) {
  k <- 0:(n - 1); a <- c(1 / sqrt(n), rep(sqrt(2 / n), n - 1))
  R <- matrix(0, n, n)
  for (ki in seq_along(k)) R[ki, ] <- a[ki] * cos(pi * k[ki] * (1:n - 0.5) / n)
  R
}

# y (vector P o matriz P x 3) -> coeficientes en la base de autovalores
transformar <- function(Y, M, N) {
  CM <- Cmat(M); CN <- Cmat(N)
  if (is.matrix(Y)) {
    out <- matrix(0, nrow(Y), ncol(Y))
    for (cc in 1:ncol(Y)) out[, cc] <- as.vector(CM %*% matrix(Y[, cc], M, N) %*% t(CN))
    out
  } else as.vector(CM %*% matrix(Y, M, N) %*% t(CN))
}

# ll marginal de 3.3: y ~ N(mu 1, s2 I + (tau I + lam L)^{-1})
ll_marg_33 <- function(mu, s2, tau, lam, Yt, w, P) {
  sk <- s2 + 1 / (tau + lam * w)
  yt <- Yt; yt[1] <- yt[1] - mu * sqrt(P)
  -P / 2 * log(2 * pi) - 0.5 * sum(log(sk)) - 0.5 * sum(yt * yt / sk)
}

# ll marginal de 3.4 (eigen-trick simetrico):
# y ~ N(1 (x) mu, I (x) Sigma + (I + lam L)^{-1} (x) Omega^{-1})
ll_marg_34 <- function(mu3, Sg3, Om, lam, Yt, w, P) {
  ew <- eigen(Om, symmetric = TRUE)
  Omh <- ew$vectors %*% diag(sqrt(ew$values)) %*% t(ew$vectors)
  H <- Omh %*% Sg3 %*% Omh
  eh <- eigen(H, symmetric = TRUE)
  rho <- eh$values; U <- eh$vectors
  yt <- Yt; yt[1, ] <- yt[1, ] - mu3 * sqrt(P)
  Acoef <- yt %*% Omh
  Al <- t(U) %*% t(Acoef)
  ck <- 1 / (1 + lam * w)
  q <- numeric(P)
  for (j in 1:3) q <- q + Al[j, ]^2 / (ck + rho[j])
  ld <- -P * log(det(Om)) + sum(sapply(1:3, function(j) sum(log(ck + rho[j]))))
  -3 * P / 2 * log(2 * pi) - 0.5 * ld - 0.5 * sum(q)
}

vech2mat <- function(v) {
  m <- matrix(0, 3, 3)
  m[1, 1] <- v[1]; m[2, 1] <- m[1, 2] <- v[2]; m[3, 1] <- m[1, 3] <- v[3]
  m[2, 2] <- v[4]; m[3, 2] <- m[2, 3] <- v[5]; m[3, 3] <- v[6]
  m
}

# Muestreador exacto del GMRF previo x | mu, tau, lam (3.3) y X | mu, Om, lam (3.4)
factorizar_A <- function(tau, lam, reja, perm = NULL) {
  Q <- as(tau * Diagonal(reja$P) + lam * reja$L, "dgCMatrix")
  if (is.null(perm)) perm <- Cholesky(Q, LDL = FALSE, perm = TRUE, Imult = 0)@perm
  ch <- Cholesky(Q, perm = perm, LDL = FALSE, Imult = 0)
  list(Lm = expand(ch)$L, perm = perm)
}
rmv_gmrf33 <- function(mu, tau, lam, reja, fac) {
  z <- rnorm(reja$P)
  u <- as.matrix(solve(t(fac$Lm), z))
  x <- numeric(reja$P); x[fac$perm + 1] <- u
  x + mu
}
rmv_gmrf34 <- function(mu3, Om, lam, reja, fac) {
  Z <- matrix(rnorm(reja$P * 3), reja$P, 3)
  M1 <- as.matrix(solve(t(fac$Lm), Z))
  ew <- eigen(Om, symmetric = TRUE)
  Gm <- ew$vectors %*% diag(1 / sqrt(ew$values)) %*% t(ew$vectors)  # Gm Gm' = Om^{-1}
  Xt <- as.matrix(M1 %*% t(Gm))
  Xp <- matrix(0, reja$P, 3); Xp[fac$perm + 1, ] <- Xt
  Xp + rep(mu3, each = reja$P)
}

# ------------------------------------------------------------------------------
# 2. VALIDACION INTERNA (rejilla pequena contra la densa)
# ------------------------------------------------------------------------------
cat("\nValidacion interna de la log-verosimilitud marginal y del muestreador...\n")
set.seed(11)
Ms <- 6; Ns <- 5; Ps <- Ms * Ns
rj <- construir_reja(Ms, Ns)
ys <- rnorm(Ps, 0.5, 0.2)
Yts <- transformar(ys, Ms, Ns)
Qs <- 2 * diag(Ps) + 7 * as.matrix(rj$L)
ll_d <- dmvnorm(ys, rep(0.4, Ps), 0.01 * diag(Ps) + solve(Qs), log = TRUE)
ll_s <- ll_marg_33(0.4, 0.01, 2, 7, Yts, rj$w, Ps)
stopifnot(abs(ll_d - ll_s) < 1e-8)

Sg3s <- matrix(c(0.01, 0.001, 0, 0.001, 0.02, 0, 0, 0, 0.015), 3, 3)
Oms <- matrix(c(2, 0.3, 0.1, 0.3, 1.5, 0.2, 0.1, 0.2, 1), 3, 3)
Ys <- matrix(rnorm(Ps * 3, 0.5, 0.2), Ps, 3)
Yt4s <- transformar(Ys, Ms, Ns)
A34 <- diag(Ps) + 7 * as.matrix(rj$L)   # I + lam L (factor espacial del 3.4)
Sg_big <- kronecker(Sg3s, diag(Ps)) + kronecker(solve(Oms), solve(A34))
ll_d4 <- dmvnorm(as.vector(Ys), rep(c(0.5, 0.4, 0.3), each = Ps), Sg_big, log = TRUE)
ll_s4 <- ll_marg_34(c(0.5, 0.4, 0.3), Sg3s, Oms, 7, Yt4s, rj$w, Ps)
stopifnot(abs(ll_d4 - ll_s4) < 1e-8)

facS <- factorizar_A(2, 7, rj)
Xs <- replicate(2000, rmv_gmrf33(0.4, 2, 7, rj, facS))
err_cov <- max(abs(cov(t(Xs)) - solve(Qs)) / sqrt(outer(diag(solve(Qs)), diag(solve(Qs)))))
stopifnot(err_cov < 0.25)   # ruido de MC con 2000 replicas
cat("  OK (ll 3.3 y 3.4 vs densa: error < 1e-8; covarianza del muestreador dentro del ruido de MC)\n")

# ------------------------------------------------------------------------------
# 3. LOG-VEROSIMILITUD Y DIC (3.3 y 3.4)
#    (con CDE_SOLO_FUNCIONES=1 se cargan las funciones sin ejecutar nada)
# ------------------------------------------------------------------------------
if (!identical(Sys.getenv("CDE_SOLO_FUNCIONES"), "1")) {
M <- 128; N <- 128; P <- M * N
reja <- construir_reja(M, N)

procesar_modelo <- function(nombre, ruta, es34) {
  obj <- readRDS(ruta)
  filas <- list(); trazas <- list()
  for (nm in names(obj$resultados)) {
    esc <- obj$resultados[[nm]]
    res <- esc$res; post <- res$post
    y <- if (es34) esc$Y else esc$y
    Yt <- transformar(y, M, N)
    Tn <- length(post)
    if (es34) {
      mu_c <- res$mu[post, , drop = FALSE]
      lam_c <- res$lambda[post]
      ll <- numeric(Tn)
      for (t in seq_len(Tn)) {
        Sg3 <- vech2mat(res$Sigma[post[t], ])
        Om <- vech2mat(res$Omega[post[t], ])
        ll[t] <- ll_marg_34(mu_c[t, ], Sg3, Om, lam_c[t], Yt, reja$w, P)
      }
      Sg_bar <- vech2mat(colMeans(res$Sigma[post, , drop = FALSE]))
      Om_bar <- vech2mat(colMeans(res$Omega[post, , drop = FALSE]))
      ll_bar <- ll_marg_34(colMeans(mu_c), Sg_bar, Om_bar, mean(lam_c), Yt, reja$w, P)
    } else {
      ll <- numeric(Tn)
      for (t in seq_len(Tn)) {
        ll[t] <- ll_marg_33(res$mu[post[t]], res$sig2[post[t]],
                            res$tau[post[t]], res$lambda[post[t]], Yt, reja$w, P)
      }
      ll_bar <- ll_marg_33(mean(res$mu[post]), mean(res$sig2[post]),
                           mean(res$tau[post]), mean(res$lambda[post]), Yt, reja$w, P)
    }
    p_D <- 2 * (ll_bar - mean(ll))
    dic <- -2 * ll_bar + 2 * p_D
    filas[[nm]] <- data.frame(
      modelo = nombre, escenario = nm, nivel_ruido = esc$sigma0,
      ll_media = mean(ll), ll_min = min(ll), ll_max = max(ll),
      ll_en_theta_barra = ll_bar, p_D = p_D, DIC = dic,
      stringsAsFactors = FALSE
    )
    trazas[[nm]] <- ll
    png(file.path(dir_ver, sprintf("ll_traza_%s_%.2f.png", nombre, esc$sigma0)),
        width = 1300, height = 380, res = 120)
    plot(ll, type = "l", col = "grey30", xlab = "Iteracion (tras calentamiento)",
         ylab = "log p(y | theta)", main = sprintf("%s: log-verosimilitud marginal (s0 = %.2f)",
                                                   nombre, esc$sigma0))
    abline(h = mean(ll), col = "#1b9e77", lwd = 2)
    abline(h = ll_bar, col = "#d95f02", lwd = 2, lty = 2)
    legend("topright", legend = c("media de la cadena", "en theta barra"),
           col = c("#1b9e77", "#d95f02"), lwd = 2, bty = "n", cex = 0.8)
    dev.off()
    cat(sprintf("  %s | %s: ll media = %.1f | p_D = %.1f | DIC = %.1f\n",
                nombre, nm, mean(ll), p_D, dic))
  }
  do.call(rbind, filas)
}

cat("\n----- Log-verosimilitud marginal y DIC -----\n")
t33 <- procesar_modelo("3.3", file.path(dir_out, "3.3", "3.3_resultados.rds"), FALSE)
t34 <- procesar_modelo("3.4", file.path(dir_out, "3.4", "3.4_resultados.rds"), TRUE)
dic_tab <- rbind(t33, t34)
write.csv(dic_tab, file.path(dir_ver, "dic.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. PROBABILIDADES POSTERIORES DE CONTRASTES
# ------------------------------------------------------------------------------
cat("\n----- Probabilidades posteriores de contrastes -----\n")
obj32 <- readRDS(file.path(dir_out, "3.2", "3.2_resultados.rds"))
obj33 <- readRDS(file.path(dir_out, "3.3", "3.3_resultados.rds"))
obj34 <- readRDS(file.path(dir_out, "3.4", "3.4_resultados.rds"))

contr <- list()
beta_c <- log(1 + sqrt(2)) / 2   # ~0.4407: temperatura critica del Ising 2D
for (nm in names(obj32$resultados)) {
  esc <- obj32$resultados[[nm]]; post <- esc$post
  b <- esc$beta[post]; r <- esc$rho[post]
  contr[[length(contr) + 1]] <- data.frame(
    modelo = "Ising", escenario = nm, contraste = "P(beta > beta_c | y)",
    prob = mean(b > beta_c), stringsAsFactors = FALSE)
  contr[[length(contr) + 1]] <- data.frame(
    modelo = "Ising", escenario = nm, contraste = "P(rho < rho0 | y)",
    prob = mean(r < esc$rho0), stringsAsFactors = FALSE)
}
for (nm in names(obj33$resultados)) {
  esc <- obj33$resultados[[nm]]; post <- esc$res$post
  s2 <- esc$res$sig2[post]; tau <- esc$res$tau[post]; lam <- esc$res$lambda[post]
  s0 <- esc$sigma0
  contr[[length(contr) + 1]] <- data.frame(
    modelo = "GMRF grises", escenario = nm,
    contraste = sprintf("P(sigma2 < %.4f | y)", s0^2),
    prob = mean(s2 < s0^2), stringsAsFactors = FALSE)
  contr[[length(contr) + 1]] <- data.frame(
    modelo = "GMRF grises", escenario = nm, contraste = "P(tau < 0.01 | y)",
    prob = mean(tau < 0.01), stringsAsFactors = FALSE)
  contr[[length(contr) + 1]] <- data.frame(
    modelo = "GMRF grises", escenario = nm, contraste = "P(lambda > 20 | y)",
    prob = mean(lam > 20), stringsAsFactors = FALSE)
}
for (nm in names(obj34$resultados)) {
  esc <- obj34$resultados[[nm]]; post <- esc$res$post
  s0 <- esc$sigma0
  s2d <- apply(esc$res$Sigma[post, , drop = FALSE], 1, function(v) diag(vech2mat(v)))
  for (cc in 1:3) {
    contr[[length(contr) + 1]] <- data.frame(
      modelo = "GMRF RGB", escenario = nm,
      contraste = sprintf("P(Sigma_%d%d < %.4f | y)", cc, cc, s0^2),
      prob = mean(s2d[cc, ] < s0^2), stringsAsFactors = FALSE)
  }
  contr[[length(contr) + 1]] <- data.frame(
    modelo = "GMRF RGB", escenario = nm, contraste = "P(lambda > 400 | y)",
    prob = mean(esc$res$lambda[post] > 400), stringsAsFactors = FALSE)
}
contr <- do.call(rbind, contr)
write.csv(contr, file.path(dir_ver, "contrastes.csv"), row.names = FALSE)
print(contr, row.names = FALSE, digits = 3)

# ------------------------------------------------------------------------------
# 5. PRIOR PREDICTIVE CHECK (3.3 y 3.4)
# ------------------------------------------------------------------------------
cat("\n----- Prior predictive check -----\n")
set.seed(20261125)
Nrep <- 30
est_img <- function(y, aristas) {
  difs <- abs(y[aristas[, 1]] - y[aristas[, 2]])
  c(media = mean(y), DE = sd(y), contraste = mean(difs),
    min = min(y), max = max(y))
}

facP <- NULL
prior33 <- matrix(0, Nrep, 5); colnames(prior33) <- c("media", "DE", "contraste", "min", "max")
for (r in 1:Nrep) {
  mu <- rnorm(1, 0.5, sqrt(0.25))
  s2 <- 1 / rgamma(1, shape = 2, rate = 0.01)
  tau <- rgamma(1, shape = 1, rate = 0.05)
  lam <- rgamma(1, shape = 2, rate = 0.1)
  facP <- factorizar_A(tau, lam, reja)
  x <- rmv_gmrf33(mu, tau, lam, reja, facP)
  y <- x + rnorm(P, 0, sqrt(s2))
  prior33[r, ] <- est_img(y, reja$aristas)
}

prior34 <- matrix(0, Nrep, 5); colnames(prior34) <- colnames(prior33)
for (r in 1:Nrep) {
  mu3 <- rnorm(3, 0.5, sqrt(0.25))
  Sg3 <- solve(rWishart(1, df = 6, Sigma = solve(0.01 * diag(3)))[, , 1])
  Om <- rWishart(1, df = 5, Sigma = 0.1 * diag(3))[, , 1]
  lam <- rgamma(1, shape = 2, rate = 0.1)
  facP <- factorizar_A(1, lam, reja)   # factor espacial de I + lam L (tau = 1)
  X <- rmv_gmrf34(mu3, Om, lam, reja, facP)
  Y <- X + matrix(rnorm(P * 3), P, 3) %*% chol(Sg3)
  prior34[r, ] <- est_img(as.vector(rowMeans(Y)), reja$aristas)
}

write.csv(rbind(data.frame(modelo = "GMRF grises", prior33),
                data.frame(modelo = "GMRF RGB", prior34)),
          file.path(dir_ver, "prior_predictivo.csv"), row.names = FALSE)

obs33 <- sapply(obj33$resultados, function(esc) est_img(esc$y, reja$aristas))
obs34 <- sapply(obj34$resultados, function(esc) est_img(rowMeans(esc$Y), reja$aristas))
png(file.path(dir_ver, "prior_predictivo.png"), width = 1300, height = 800, res = 120)
op <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (est in c("media", "DE", "contraste")) {
  hist(prior33[, est], breaks = 15, col = "grey80", border = "white", main = paste("3.3 previa:", est),
       xlab = est, xlim = range(c(prior33[, est], obs33[est, ])))
  abline(v = obs33[est, ], col = c("#1b9e77", "#d95f02", "#7570b3"), lwd = 2)
  legend("topright", legend = paste0("y s0=", c(0.10, 0.05, 0.20)),
         col = c("#1b9e77", "#d95f02", "#7570b3"), lwd = 2, bty = "n", cex = 0.7)
}
for (est in c("media", "DE", "contraste")) {
  hist(prior34[, est], breaks = 15, col = "grey80", border = "white", main = paste("3.4 previa:", est),
       xlab = est, xlim = range(c(prior34[, est], obs34[est, ])))
  abline(v = obs34[est, ], col = c("#1b9e77", "#d95f02", "#7570b3"), lwd = 2)
}
par(op); dev.off()
cat("  Simuladas", Nrep, "imagenes previas por modelo; cobertura de escalas en prior_predictivo.png\n")
for (est in c("media", "DE", "contraste")) {
  cat(sprintf("  %s | 3.3 previa [%.3f, %.3f] vs observado (%.3f-%.3f) | 3.4 previa [%.3f, %.3f] vs obs (%.3f-%.3f)\n",
              est, min(prior33[, est]), max(prior33[, est]), min(obs33[est, ]), max(obs33[est, ]),
              min(prior34[, est]), max(prior34[, est]), min(obs34[est, ]), max(obs34[est, ])))
}

# ------------------------------------------------------------------------------
# 6. POSTERIOR PREDICTIVE CHECKS (3.3 y 3.4)
# ------------------------------------------------------------------------------
cat("\n----- Posterior predictive checks -----\n")
set.seed(20261126)
Nrep <- 30
ppp_modelo <- function(nombre, obj, es34) {
  filas <- list()
  for (nm in names(obj$resultados)) {
    esc <- obj$resultados[[nm]]
    res <- esc$res; post <- res$post
    y <- if (es34) esc$Y else esc$y
    y_obs <- if (es34) as.vector(rowMeans(y)) else y
    t_obs <- est_img(y_obs, reja$aristas)
    idx <- post[round(seq(1, length(post), length.out = Nrep))]
    t_rep <- matrix(0, Nrep, 5); colnames(t_rep) <- names(t_obs)
    for (r in seq_along(idx)) {
      t <- idx[r]
      if (es34) {
        Sg3 <- vech2mat(res$Sigma[t, ])
        Om <- vech2mat(res$Omega[t, ])
        fac <- factorizar_A(1, res$lambda[t], reja)
        X <- rmv_gmrf34(res$mu[t, ], Om, res$lambda[t], reja, fac)
        Y <- X + matrix(rnorm(P * 3), P, 3) %*% chol(Sg3)
        y_star <- as.vector(rowMeans(Y))
      } else {
        fac <- factorizar_A(res$tau[t], res$lambda[t], reja)
        x <- rmv_gmrf33(res$mu[t], res$tau[t], res$lambda[t], reja, fac)
        y_star <- x + rnorm(P, 0, sqrt(res$sig2[t]))
      }
      t_rep[r, ] <- est_img(y_star, reja$aristas)
    }
    ppp <- colMeans(sweep(t_rep, 2, t_obs, "<="))
    filas[[nm]] <- data.frame(
      modelo = nombre, escenario = nm, nivel_ruido = esc$sigma0,
      estadistico = names(t_obs), t_observado = as.numeric(t_obs),
      t_rep_media = colMeans(t_rep),
      t_rep_p2.5 = apply(t_rep, 2, quantile, 0.025),
      t_rep_p97.5 = apply(t_rep, 2, quantile, 0.975),
      ppp = as.numeric(ppp), stringsAsFactors = FALSE
    )
    png(file.path(dir_ver, sprintf("ppp_%s_%.2f.png", nombre, esc$sigma0)),
        width = 1300, height = 800, res = 120)
    op <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
    for (est in c("media", "DE", "contraste")) {
      hist(t_rep[, est], breaks = 15, col = "grey80", border = "white",
           main = sprintf("%s s0=%.2f: t = %s (ppp = %.2f)", nombre, esc$sigma0, est,
                          ppp[est]),
           xlab = est, xlim = range(c(t_rep[, est], t_obs[est])))
      abline(v = t_obs[est], col = "#d95f02", lwd = 2)
    }
    par(op); dev.off()
    cat(sprintf("  %s | %s: ppp media = %.2f | DE = %.2f | contraste = %.2f\n",
                nombre, nm, ppp["media"], ppp["DE"], ppp["contraste"]))
  }
  do.call(rbind, filas)
}

cat("(replicas desde el posterior; ppp = Pr(t(y*) <= t(y)))\n")
p33 <- ppp_modelo("3.3", obj33, FALSE)
p34 <- ppp_modelo("3.4", obj34, TRUE)
write.csv(rbind(p33, p34), file.path(dir_ver, "ppp.csv"), row.names = FALSE)

cat("\nListo. Archivos en:", dir_ver, "\n")
}
