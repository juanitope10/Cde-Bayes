# ==============================================================================
# SECCIÓN 3.4 - IMAGEN RGB: MODELO GAUSSIANO MULTIVARIADO
# ------------------------------------------------------------------------------
# Modelo:
#   Y_i | X_i, SIGMA ~ N_3(X_i, SIGMA) independientes
#   x | mu, Omega, lambda ~ N_{3P}(1_P (x) mu, Q^{-1}),
#   Q = (I_P + lambda*L) (x) Omega     [producto de Kronecker]
#   SIGMA ~ IW_3(nu_eps, S_eps) | mu ~ N_3(m_mu, V_mu)
#   Omega ~ W_3(nu_om, S_om)          | lambda ~ G(a_l, b_l)
#
# Convenciones:
#   W_p(n, V):  p(W) prop |W|^{(n-p-1)/2} exp{-tr(V^{-1} W)/2},  E[W] = n V
#   IW_p(n, S): p(S) prop |S|^{-(n+p+1)/2} exp{-tr(S Sg^{-1})/2}, E[Sg] = S/(n-p-1)
#
# Condicionales completas (ver derivación en el informe):
#   X_i | .  ~ N_3( Lam_i^{-1} h_i, Lam_i^{-1} ),
#        Lam_i = SIGMA^{-1} + (1 + lambda*d_i) Omega
#        h_i   = SIGMA^{-1} y_i + Omega ( mu + lambda*sum_{j in N_i} x_j )
#   SIGMA | . ~ IW_3( nu_eps + P, S_eps + sum_i (y_i - x_i)(y_i - x_i)^T )
#   mu | .    ~ N_3( V_n m_n, V_n ),  V_n = (V_mu^{-1} + P*Omega)^{-1},
#                    m_n = V_n ( V_mu^{-1} m_mu + Omega*sum_i x_i )
#   Omega | . ~ W_3( nu_om + P, (S_om^{-1} + M)^{-1} ),
#        M = sum_i (x_i - mu)(x_i - mu)^T + lambda*sum_{i~j} (x_i - x_j)(x_i - x_j)^T
#   lambda | . no es conjugada por |I + lambda*L|^{3/2}:
#        p(lam|.) prop lam^{a_l} exp(-b_l lam) |I + lam L|^{3/2}
#                   exp{ -lam/2 * sum_{i~j} (x_i-x_j)^T Omega (x_i-x_j) }
#        -> Metropolis-Hastings en escala log; log|I + lam L| por Cholesky disperso
#
# Actualización de x: tablero de ajedrez (ver 3.3). Como Lam_i depende de i solo
# a través de d_i en {2,3,4}, se factorizan solo TRES matrices 3x3 por iteración.
# ==============================================================================
suppressPackageStartupMessages({
  library(Matrix)
  library(png)
  library(coda)
})

# ------------------------------------------------------------------------------
# 0. CONFIGURACIÓN
# ------------------------------------------------------------------------------
M <- 128; N <- 128; P <- M * N

# --- Rutas -------------------------------------------------------------------
# Las rutas del proyecto las define el archivo de configuración del sistema
# operativo: config_linux.R (Linux) o config_windows.R (Windows), ambos en Code/.
# Allí se definen dir_raiz, dir_code, dir_img y dir_out; aquí solo se usan.
if (!exists("dir_img")) {
  archivo_config <- if (identical(.Platform$OS.type, "windows")) "config_windows.R" else "config_linux.R"
  dirs <- c(".", "Code", "..", file.path("..", "Code"))
  # junto al script que se está ejecutando (Rscript --file= y source()/Source)
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
  if (length(f) == 1L && !is.na(f)) dirs <- c(dirname(gsub("~+~", " ", f, fixed = TRUE)), dirs)
  for (i in seq_len(sys.nframe())) {
    g <- sys.frame(i)$ofile
    if (!is.null(g) && nzchar(g)) {
      dirs <- c(dirname(gsub("~+~", " ", g, fixed = TRUE)), dirs)
      break
    }
  }
  cargado <- FALSE
  for (d in unique(dirs)) {
    if (file.exists(file.path(d, archivo_config))) {
      source(file.path(d, archivo_config))
      cargado <- TRUE
      break
    }
  }
  if (!cargado) {
    stop("No se encontró '", archivo_config, "'.\n",
         "Ejecute el script desde la carpeta del proyecto o desde 'Code/',\n",
         "o defina dir_img y dir_out a mano antes de ejecutarlo.")
  }
}
dir_out <- file.path(dir_out, "3.4")   # resultados de este punto en su subcarpeta
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
cat("Imágenes en:", dir_img, "\nResultados en:", dir_out, "\n")

RAPIDO <- toupper(Sys.getenv("CDE_RAPIDO", "FALSE")) %in% c("TRUE", "1", "YES")
# "cholesky" (por defecto, estándar del curso) o "espectral" (exacto, más rápido).
# Ambos se validan numéricamente al inicio (coinciden a 1e-8); para cadenas
# largas se recomienda CDE_LOGDET=espectral, ~20x más rápido y exactamente igual.
metodo_logdet <- Sys.getenv("CDE_LOGDET", "cholesky")

config <- list(
  # Hiperparámetros: previas débilmente informativas sobre escalas de [0,1]
  nu_eps = 6, S_eps = 0.01 * diag(3),   # SIGMA ~ IW(6, 0.01 I): E[SIGMA] = 0.005 I
  m_mu = c(0.5, 0.5, 0.5), V_mu = 0.25 * diag(3),   # mu ~ N_3(0.5*1, 0.25 I)
  nu_om = 5, S_om = 0.1 * diag(3),      # Omega ~ W(5, 0.1 I): E[Omega] = nu*S = 0.5 I
                                        #   (errata corregida: antes decia 50 I; n*V, no n/V)
  a_lambda = 2, b_lambda = 0.1,         # lambda ~ G(2, 0.1): E[lambda] = 20
  # Algoritmo (B grande: la cresta casi-invariante lambda*Omega exige n_eff alto)
  B = 24000, burn = 2000,
  thin = 40,                            # cadena de píxeles adelgazada (memoria)
  s_escala = 0.5,
  adapt = TRUE, bloques_adapt = 25,
  sigma0_levels = c(0.10, 0.05, 0.20)
)
if (RAPIDO) {
  config$B <- 400; config$burn <- 200; config$thin <- 4; config$sigma0_levels <- 0.10
}
semilla <- 2026

# ------------------------------------------------------------------------------
# 1. ESTRUCTURA ESPACIAL (matrices dispersas) Y ARISTAS
# ------------------------------------------------------------------------------
construir_estructura <- function(M, N) {
  idx <- matrix(seq_len(M * N), M, N)
  i <- c(as.vector(idx[, -N]), as.vector(idx[-M, ]))   # arista i-j (cada una una vez)
  j <- c(as.vector(idx[, -1]), as.vector(idx[-1, ]))
  A <- sparseMatrix(i = c(i, j), j = c(j, i), x = 1, dims = c(M * N, M * N))
  d <- rowSums(A)
  L <- Diagonal(x = d) - A
  list(A = A, L = L, d = d, ei = i, ej = j)
}
preparar_tablero <- function(A, M, N) {
  k <- seq_len(M * N)
  r <- (k - 1) %% M + 1
  c <- (k - 1) %/% M + 1
  grupos <- list(k[(r + c) %% 2 == 0], k[(r + c) %% 2 == 1])
  stopifnot(sum(A[grupos[[1]], grupos[[1]]]) == 0, sum(A[grupos[[2]], grupos[[2]]]) == 0)
  list(idx = grupos, A_sub = lapply(grupos, function(g) A[g, , drop = FALSE]))
}

est <- construir_estructura(M, N)
A <- est$A; L <- est$L; d <- est$d; ei <- est$ei; ej <- est$ej
tab <- preparar_tablero(A, M, N)
rm(est)

# ------------------------------------------------------------------------------
# 2. LOG-DETERMINANTE DE I + lambda*L  (con factorización de Cholesky dispersa)
# ------------------------------------------------------------------------------
perm_C <- Cholesky(as(1 * Diagonal(P) + 1 * L, "dgCMatrix"), LDL = FALSE,
                   perm = TRUE, Imult = 0)@perm

log_det_IL_cholesky <- function(lam) {
  Cm <- as(1 * Diagonal(P) + lam * L, "dgCMatrix")
  ch <- Cholesky(Cm, perm = perm_C, LDL = FALSE, Imult = 0)
  dt <- determinant(ch, logarithm = TRUE, sqrt = FALSE)
  as.numeric(dt$modulus) * dt$sign
}
# Método espectral exacto (opcional), validado contra Cholesky abajo
eig_L <- outer(4 * sin(pi * (0:(M - 1)) / (2 * M))^2,
               4 * sin(pi * (0:(N - 1)) / (2 * N))^2, "+")
log_det_IL_espectral <- function(lam) sum(log(1 + lam * eig_L))

log_det_IL <- if (metodo_logdet == "espectral") log_det_IL_espectral else log_det_IL_cholesky

for (lam in c(0.5, 5, 30)) {
  a1 <- log_det_IL_cholesky(lam); a2 <- log_det_IL_espectral(lam)
  stopifnot(abs(a1 - a2) / abs(a1) < 1e-8)
}
cat(sprintf("log|I + lam*L| validado | método activo: %s\n", metodo_logdet))

# ------------------------------------------------------------------------------
# 3. MUESTREO DE WISHART / WISHART INVERSA (descomposición de Bartlett)
# ------------------------------------------------------------------------------
# W ~ W_p(n, V): p(W) prop |W|^{(n-p-1)/2} exp{-tr(V^{-1} W)/2}, E[W] = n V
# Muestreo: W = Rw^T (Z Z^T) Rw con V = Rw^T Rw
rwishart <- function(n, V) {
  p <- nrow(V)
  Rw <- chol(V)
  Z <- matrix(rnorm(p * n), p, n)
  t(Rw) %*% tcrossprod(Z) %*% Rw
}
# Sigma ~ IW_p(n, S): p(Sigma) prop |Sigma|^{-(n+p+1)/2} exp{-tr(S Sigma^{-1})/2}
# Muestreo: Sigma = R^T (Z Z^T)^{-1} R con S = R^T R
riwishart <- function(n, S) {
  p <- nrow(S)
  R <- chol(S)
  Z <- matrix(rnorm(p * n), p, n)
  G <- tcrossprod(Z)
  t(R) %*% solve(G) %*% R
}

# ------------------------------------------------------------------------------
# 4. DATOS: contaminación RGB y_i = x0_i + N_3(0, sigma0^2 I)  (sin recortar)
# ------------------------------------------------------------------------------
arr <- readPNG(file.path(dir_img, "astronaut_rgb.png"))
stopifnot(dim(arr)[1] == M, dim(arr)[2] == N, dim(arr)[3] >= 3)
X0 <- cbind(R = as.vector(arr[, , 1]), G = as.vector(arr[, , 2]), B = as.vector(arr[, , 3]))
cat(sprintf("Imagen RGB: medias por canal = (%.3f, %.3f, %.3f)\n",
            mean(X0[, 1]), mean(X0[, 2]), mean(X0[, 3]))
)
canales <- c("R", "G", "B")

# Empaquetado de matrices simétricas 3x3 (6 valores distintos, orden por columnas)
vech_sym <- function(Am) c(Am[1, 1], Am[2, 1], Am[3, 1], Am[2, 2], Am[3, 2], Am[3, 3])
unvech_sym <- function(v) {
  Am <- matrix(0, 3, 3)
  Am[1, 1] <- v[1]; Am[2, 1] <- Am[1, 2] <- v[2]; Am[3, 1] <- Am[1, 3] <- v[3]
  Am[2, 2] <- v[4]; Am[3, 2] <- Am[2, 3] <- v[5]; Am[3, 3] <- v[6]
  Am
}
etiquetas_sim <- c("11", "21", "31", "22", "32", "33")

# ------------------------------------------------------------------------------
# 5. PASOS MCMC
# ------------------------------------------------------------------------------
# 5.1 X | resto: tablero de ajedrez; solo 3 factorizaciones 3x3 por barrido
muestrear_X <- function(X, Y, mu, Omega, Sigma, lambda) {
  iSigma <- solve(Sigma)
  grados <- 2:4
  Rch <- lapply(grados, function(g) {
    Lam <- iSigma + (1 + lambda * g) * Omega
    chol(Lam)                       # Lam = t(R) %*% R
  })
  for (g in 1:2) {
    ig <- tab$idx[[g]]
    S  <- as.matrix(tab$A_sub[[g]] %*% X)            # filas: sumas de vecinos
    H  <- Y[ig, , drop = FALSE] %*% iSigma +
          (matrix(1, length(ig), 1) %*% mu + lambda * S) %*% Omega
    dg <- d[ig]
    for (k in seq_along(grados)) {
      sel <- which(dg == grados[k])                  # posiciones locales en ig
      if (!length(sel)) next
      Rk <- Rch[[k]]
      Mn <- H[sel, , drop = FALSE] %*% chol2inv(Rk)
      Zn <- matrix(rnorm(3 * length(sel)), length(sel), 3)
      En <- t(backsolve(Rk, t(Zn)))                  # filas con cov = Lam^{-1}
      X[ig[sel], ] <- Mn + En
    }
  }
  X
}

# 5.2 SIGMA | resto ~ IW(nu_eps + P, S_eps + sum (y_i - x_i)(y_i - x_i)^T)
muestrear_Sigma <- function(X, Y, nu_eps, S_eps) {
  R <- Y - X
  riwishart(nu_eps + P, S_eps + crossprod(R))
}

# 5.3 mu | resto ~ N_3(V_n m_n, V_n)
muestrear_mu <- function(X, Omega, m_mu, V_mu) {
  iVmu <- solve(V_mu)
  Vn <- solve(iVmu + P * Omega)
  mn <- Vn %*% (iVmu %*% m_mu + Omega %*% colSums(X))
  c(mn) + drop(t(chol(Vn)) %*% rnorm(3))
}

# 5.4 Omega | resto ~ W(nu_om + P, V) con V^{-1} = S_om^{-1} + M
muestrear_Omega <- function(X, mu, lambda, nu_om, S_om) {
  Z <- sweep(X, 2, mu, "-")
  D <- X[ei, , drop = FALSE] - X[ej, , drop = FALSE]
  Mm <- crossprod(Z) + lambda * crossprod(D)
  V <- solve(solve(S_om) + Mm)
  rwishart(nu_om + P, V)
}

# 5.5 lambda | resto: MH en escala log (ver encabezado)
paso_mh_lambda <- function(eta, s, a_l, b_l, logdet, quad) {
  eta_p <- eta + rnorm(1, 0, s)
  lam_p <- exp(eta_p)
  logdet_p <- log_det_IL(lam_p)
  lp   <- a_l * eta   - b_l * exp(eta) + 1.5 * logdet   - 0.5 * exp(eta) * quad
  lp_p <- a_l * eta_p - b_l * lam_p    + 1.5 * logdet_p - 0.5 * lam_p    * quad
  if (log(runif(1)) < lp_p - lp) list(eta = eta_p, lam = lam_p, logdet = logdet_p, acepta = TRUE)
  else                           list(eta = eta,   lam = exp(eta), logdet = logdet,   acepta = FALSE)
}

# ------------------------------------------------------------------------------
# 6. MCMC COMPLETO (un nivel de ruido)
# ------------------------------------------------------------------------------
ejecutar_mcmc <- function(Y, cfg, semilla, verbose = TRUE) {
  set.seed(semilla)
  B <- cfg$B; ncat <- floor(B / 10)

  # Inicialización
  X <- Y
  mu <- colMeans(Y); Sigma <- 0.01 * diag(3); Omega <- 10 * diag(3)
  lambda <- 10; eta <- log(lambda); s <- cfg$s_escala

  MU <- matrix(0, B, 3); LAM <- numeric(B)
  SIG <- matrix(0, B, 6); OME <- matrix(0, B, 6)
  ACC <- logical(B)

  n_th <- ceiling((B - cfg$burn) / cfg$thin)
  X_th <- array(0, dim = c(n_th, P, 3))
  SX <- matrix(0, P, 3); SX2 <- matrix(0, P, 3); n_post <- 0
  i_th <- 0

  logdet <- log_det_IL(lambda)
  t0 <- Sys.time()

  for (b in seq_len(B)) {
    X     <- muestrear_X(X, Y, mu, Omega, Sigma, lambda)      # 1) X | resto
    Sigma <- muestrear_Sigma(X, Y, cfg$nu_eps, cfg$S_eps)     # 2) SIGMA | resto
    mu    <- muestrear_mu(X, Omega, cfg$m_mu, cfg$V_mu)       # 3) mu | resto
    Omega <- muestrear_Omega(X, mu, lambda, cfg$nu_om, cfg$S_om)  # 4) Omega | resto

    D <- X[ei, , drop = FALSE] - X[ej, , drop = FALSE]        # 5) lambda | resto
    quad <- sum((D %*% Omega) * D)
    pm <- paso_mh_lambda(eta, s, cfg$a_lambda, cfg$b_lambda, logdet, quad)
    eta <- pm$eta; lambda <- pm$lam; ACC[b] <- pm$acepta
    if (pm$acepta) logdet <- pm$logdet

    MU[b, ] <- mu; LAM[b] <- lambda
    SIG[b, ] <- vech_sym(Sigma); OME[b, ] <- vech_sym(Omega)

    if (b > cfg$burn) {
      SX <- SX + X; SX2 <- SX2 + X^2; n_post <- n_post + 1
      if ((b - cfg$burn) %% cfg$thin == 0 && i_th < n_th) {
        i_th <- i_th + 1; X_th[i_th, , ] <- X
      }
    }

    # Adaptación de la propuesta de lambda SOLO en el calentamiento (blanco ~30%)
    if (isTRUE(cfg$adapt) && b <= cfg$burn && b %% cfg$bloques_adapt == 0) {
      acc <- mean(ACC[(b - cfg$bloques_adapt + 1):b])
      s <- s * min(max(exp(1.5 * (acc - 0.30)), 0.5), 2.0)
    }

    if (verbose && (b %% ncat == 0)) {
      cat(sprintf("%.0f%% completado...\n", 100 * b / B))
      cat(sprintf("   iter %d | lam = %.3f | mu = (%.3f, %.3f, %.3f) | acept(lam) = %.1f%% | s = %.3f | %.0fs\n",
                  b, lambda, mu[1], mu[2], mu[3], 100 * mean(ACC[1:b]), s,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
  }

  post <- (cfg$burn + 1):B
  list(mu = MU, lambda = LAM, Sigma = SIG, Omega = OME, acc = ACC,
       x_media = SX / n_post, x2_media = SX2 / n_post,
       X_th = X_th, post = post, s_escala = s,
       tiempo = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

# ------------------------------------------------------------------------------
# 7. DIAGNÓSTICOS Y MÉTRICAS
# ------------------------------------------------------------------------------
diagnostico_param <- function(cadena, nombre, etiqueta = nombre) {
  ch   <- mcmc(cadena)
  neff <- as.numeric(effectiveSize(ch))
  dsd  <- sd(cadena); media <- mean(cadena)
  emc  <- dsd / sqrt(neff)
  cvmc <- 100 * emc / abs(media)
  q    <- quantile(cadena, probs = c(0.025, 0.975))
  nota <- if (abs(media) < 0.01) "  (media ~ 0: CVMC no es interpretable)" else
          if (cvmc >= 5) "  <-- AVISO CVMC >= 5%" else ""
  cat(sprintf("%-14s Media = %9.5f | DE = %8.5f | n_eff = %6.0f | EMC = %8.5f | CVMC = %8.3f%% | IC95%% = [%.5f, %.5f]%s\n",
              nombre, media, dsd, neff, emc, cvmc, q[1], q[2], nota))
  c(Media = media, DE = dsd, EMC = emc, n_eff = neff, `CVMC(%)` = cvmc,
    IC95_inf = unname(q[1]), IC95_sup = unname(q[2]))
}

panel_traza_acf <- function(cadena, etiqueta) {
  plot(cadena, type = "l", col = "grey30", xlab = "Iteración", ylab = etiqueta,
       main = paste("Traza de", etiqueta))
  acf(cadena, lag.max = 60, main = paste("ACF de", etiqueta))
}

dibujar <- function(archivo, ancho, alto, fun) {
  png(file.path(dir_out, archivo), width = ancho, height = alto, res = 110)
  on.exit(dev.off(), add = TRUE)
  fun()
  if (interactive()) {
    dev.new(width = ancho / 110, height = alto / 110, noRStudioGD = TRUE)
    fun()
  }
}

clip_arr <- function(mat3) {
  a <- array(0, dim = c(M, N, 3))
  for (k in 1:3) a[, , k] <- matrix(pmin(pmax(mat3[, k], 0), 1), M, N)
  a
}

panel_escenario <- function(X0, Y, res, sigma0) {
  op <- par(mfrow = c(1, 4), mar = c(1, 1, 2.5, 1)); on.exit(par(op), add = TRUE)
  plot(as.raster(clip_arr(X0)), main = "Original")
  plot(as.raster(clip_arr(Y)), main = sprintf("Contaminada (s0 = %.2f)", sigma0))
  plot(as.raster(clip_arr(res$x_media)), main = "Restaurada (media posterior)")
  dsd <- sqrt(pmax(res$x2_media - res$x_media^2, 0))
  dsd_prom <- rowMeans(dsd)
  image(1:N, 1:M, t(matrix(dsd_prom, M, N)[M:1, ]), col = hcl.colors(100, "Inferno"),
        axes = FALSE, asp = 1, main = "DE posterior (promedio canales)")
}

# ------------------------------------------------------------------------------
# 8. ESCENARIO COMPLETO (un nivel de ruido)
# ------------------------------------------------------------------------------
correr_escenario <- function(X0, sigma0, cfg, semilla) {
  cat(sprintf("\n================ sigma_0 = %.2f ================\n", sigma0))
  set.seed(semilla)
  Y <- X0 + matrix(rnorm(3 * P, 0, sigma0), P, 3)   # sin recortar
  cat(sprintf("Ruido real por canal (DE): %.4f | %.4f | %.4f\n",
              sd(Y[, 1] - X0[, 1]), sd(Y[, 2] - X0[, 2]), sd(Y[, 3] - X0[, 3])))

  res <- ejecutar_mcmc(Y, cfg, semilla + 1)
  post <- res$post

  cat("\n--- Resumen posterior de los parámetros escalares (post burn-in) ---\n")
  tab_res <- NULL
  # mu (3)
  for (k in 1:3) {
    r <- diagnostico_param(res$mu[post, k], paste0("mu_", canales[k]))
    tab_res <- rbind(tab_res, `mu` = r)
  }
  rownames(tab_res) <- paste0("mu_", canales)
  # SIGMA (6 valores distintos)
  for (k in 1:6) {
    r <- diagnostico_param(res$Sigma[post, k], paste0("Sigma_", etiquetas_sim[k]))
    tab_res <- rbind(tab_res, r)
  }
  rownames(tab_res)[4:9] <- paste0("Sigma_", etiquetas_sim)
  # Omega (6 valores distintos)
  for (k in 1:6) {
    r <- diagnostico_param(res$Omega[post, k], paste0("Omega_", etiquetas_sim[k]))
    tab_res <- rbind(tab_res, r)
  }
  rownames(tab_res)[10:15] <- paste0("Omega_", etiquetas_sim)
  # lambda (1)
  r <- diagnostico_param(res$lambda[post], "lambda")
  tab_res <- rbind(tab_res, lambda = r)

  cat(sprintf("\nTasa de aceptación MH de lambda (post burn-in): %.1f%%\n",
              100 * mean(res$acc[post])))

  # SIGMA real vs posterior
  cat("\n--- SIGMA: real vs media posterior (por canal, diagonal) ---\n")
  Sig_post <- unvech_sym(colMeans(res$Sigma[post, ]))
  print(round(rbind(real = diag(sigma0^2 * diag(3)), posterior = diag(Sig_post)), 5))
  cat("Correlación de ruido real = 0 (I_3) | correlaciones posteriores:\n")
  Dp <- diag(sqrt(diag(Sig_post)))
  print(round(solve(Dp) %*% Sig_post %*% solve(Dp), 3))

  # Diagnóstico de píxeles (cadena de X_i adelgazada)
  cat("\n--- Diagnóstico de los X_i (cadena de píxeles adelgazada) ---\n")
  X_th <- res$X_th
  neff_min <- Inf; cvmc_all <- c()
  for (k in 1:3) {
    Ch <- matrix(X_th[, , k], nrow(X_th), P)
    ne <- as.numeric(effectiveSize(mcmc(Ch)))
    ne[!is.finite(ne) | ne <= 0] <- nrow(Ch)
    emc <- apply(Ch, 2, sd) / sqrt(ne)
    cv  <- 100 * emc / pmax(abs(res$x_media[, k]), 0.05)
    neff_min <- min(neff_min, ne)
    cvmc_all <- c(cvmc_all, cv)
    cat(sprintf("Canal %s: n_eff min = %.0f | mediana = %.0f | CVMC medio = %.3f%% | max = %.3f%%\n",
                canales[k], min(ne), median(ne), mean(cv), max(cv)))
  }

  # Incertidumbre: IC95% para una selección de píxeles
  x_sd <- sqrt(pmax(res$x2_media - res$x_media^2, 0))
  dsd_prom <- rowMeans(x_sd)
  idx_sel <- unique(c(1, M, (64 - 1) * M + 64, P, (20 - 1) * M + 20, (96 - 1) * M + 32,
                      order(dsd_prom, decreasing = TRUE)[1:3]))
  cat("\n--- IC95% para una selección de píxeles (canal R) ---\n")
  sel_list <- list()
  for (k in 1:3) {
    Qk <- apply(matrix(X_th[, , k], nrow(X_th), P)[, idx_sel, drop = FALSE], 2,
                quantile, probs = c(0.025, 0.975))
    sel_list[[k]] <- data.frame(
      canal = canales[k], pixel = idx_sel,
      fila = (idx_sel - 1) %% M + 1, col = (idx_sel - 1) %/% M + 1,
      x0 = X0[idx_sel, k], y = Y[idx_sel, k],
      media = res$x_media[idx_sel, k], DE = x_sd[idx_sel, k],
      IC95_inf = Qk[1, ], IC95_sup = Qk[2, ],
      dentro_IC95 = X0[idx_sel, k] >= Qk[1, ] & X0[idx_sel, k] <= Qk[2, ]
    )
  }
  sel <- do.call(rbind, sel_list)
  sel_r <- sel[sel$canal == "R", ]
  num <- vapply(sel_r, is.numeric, logical(1))
  sel_r[num] <- lapply(sel_r[num], round, 4)
  print(sel_r, row.names = FALSE)

  # Cobertura empírica por canal
  cobertura <- sapply(1:3, function(k) {
    Ch <- matrix(X_th[, , k], nrow(X_th), P)
    Qk <- apply(Ch, 2, quantile, probs = c(0.025, 0.975))
    mean(X0[, k] >= Qk[1, ] & X0[, k] <= Qk[2, ])
  })
  cat(sprintf("Cobertura empírica IC95%% por canal: R = %.1f%% | G = %.1f%% | B = %.1f%%\n",
              100 * cobertura[1], 100 * cobertura[2], 100 * cobertura[3]))

  # Calidad de la restauración: global y por canal
  mse_c <- colMeans((res$x_media - X0)^2)
  mse   <- mean((res$x_media - X0)^2)
  mse_y <- mean((Y - X0)^2)
  snr_c <- 10 * log10(colMeans(X0^2) / mse_c)
  snr   <- 10 * log10(mean(X0^2) / mse)
  psnr  <- 10 * log10(1 / mse)
  snr_y <- 10 * log10(mean(X0^2) / mse_y)
  cat(sprintf("MSE = %.5f | SNR = %.2f dB | PSNR = %.2f dB (global)\n", mse, snr, psnr))
  cat(sprintf("   por canal: MSE = (%.5f, %.5f, %.5f) | SNR = (%.2f, %.2f, %.2f) dB\n",
              mse_c[1], mse_c[2], mse_c[3], snr_c[1], snr_c[2], snr_c[3]))
  cat(sprintf("   (imagen contaminada: MSE = %.5f | SNR = %.2f dB)\n", mse_y, snr_y))

  metricas <- data.frame(
    modelo = "GMRF RGB", nivel_ruido = sigma0, tiempo_s = res$tiempo,
    neff_min = min(c(neff_min, tab_res[, "n_eff"])),
    # CVMC medio: píxeles + escalares con media apreciable (|media| >= 0.01)
    cvmc_medio = mean(c(cvmc_all,
                        tab_res[abs(tab_res[, "Media"]) >= 0.01, "CVMC(%)"])),
    mse = mse, snr = snr, psnr = psnr,
    mse_R = mse_c[1], mse_G = mse_c[2], mse_B = mse_c[3],
    snr_R = snr_c[1], snr_G = snr_c[2], snr_B = snr_c[3],
    mse_contaminada = mse_y, snr_contaminada = snr_y,
    lambda_media = mean(res$lambda[post]),
    mu_R = mean(res$mu[post, 1]), mu_G = mean(res$mu[post, 2]), mu_B = mean(res$mu[post, 3]),
    sigma2_R = Sig_post[1, 1], sigma2_G = Sig_post[2, 2], sigma2_B = Sig_post[3, 3],
    acept_mh = mean(res$acc[post]),
    cobertura_R = cobertura[1], cobertura_G = cobertura[2], cobertura_B = cobertura[3],
    stringsAsFactors = FALSE
  )

  list(sigma0 = sigma0, Y = Y, res = res, tab_res = tab_res, sel = sel,
       metricas = metricas, Sig_post = Sig_post, x_sd = x_sd)
}

# ------------------------------------------------------------------------------
# 9. EJECUCIÓN: tres niveles de ruido
# ------------------------------------------------------------------------------
resultados <- list()
for (k in seq_along(config$sigma0_levels)) {
  s0 <- config$sigma0_levels[k]
  r <- correr_escenario(X0, s0, config, semilla + 100 * k)
  resultados[[sprintf("sigma0_%.2f", s0)]] <- r

  # Diagnósticos de traza/ACF de todos los escalares del escenario
  post <- r$res$post
  dibujar(sprintf("diag_mu_%.2f.png", s0), 1300, 900, function() {
    op <- par(mfrow = c(3, 2), mar = c(4, 4, 3, 1)); on.exit(par(op), add = TRUE)
    for (j in 1:3) panel_traza_acf(r$res$mu[post, j], paste("mu", canales[j]))
  })
  dibujar(sprintf("diag_Sigma_%.2f.png", s0), 1300, 1300, function() {
    op <- par(mfrow = c(6, 2), mar = c(4, 4, 3, 1)); on.exit(par(op), add = TRUE)
    et <- c("Sigma[1,1]", "Sigma[1,2]", "Sigma[1,3]", "Sigma[2,2]", "Sigma[2,3]", "Sigma[3,3]")
    for (j in 1:6) panel_traza_acf(r$res$Sigma[post, j], et[j])
  })
  dibujar(sprintf("diag_Omega_%.2f.png", s0), 1300, 1300, function() {
    op <- par(mfrow = c(6, 2), mar = c(4, 4, 3, 1)); on.exit(par(op), add = TRUE)
    et <- c("Omega[1,1]", "Omega[1,2]", "Omega[1,3]", "Omega[2,2]", "Omega[2,3]", "Omega[3,3]")
    for (j in 1:6) panel_traza_acf(r$res$Omega[post, j], et[j])
  })
  dibujar(sprintf("diag_lambda_%.2f.png", s0), 1200, 450, function() {
    op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1)); on.exit(par(op), add = TRUE)
    panel_traza_acf(r$res$lambda[post], "lambda")
  })
  dibujar(sprintf("panel_rgb_%.2f.png", s0), 1800, 480, function() {
    panel_escenario(X0, r$Y, r$res, s0)
  })
}

dibujar("comparacion_rgb.png", 1600, 1000, function() {
  op <- par(mfrow = c(length(resultados), 3), mar = c(1, 1, 2.5, 1)); on.exit(par(op), add = TRUE)
  for (r in resultados) {
    plot(as.raster(clip_arr(X0)), main = "Original")
    plot(as.raster(clip_arr(r$Y)), main = sprintf("Contaminada s0 = %.2f", r$sigma0))
    plot(as.raster(clip_arr(r$res$x_media)), main = sprintf("Restaurada s0 = %.2f", r$sigma0))
  }
})

# ------------------------------------------------------------------------------
# 10. TABLA COMPARATIVA Y GUARDADO
# ------------------------------------------------------------------------------
resumen <- do.call(rbind, lapply(resultados, `[[`, "metricas"))
cat("\n================ COMPARACIÓN ENTRE NIVELES DE RUIDO (RGB) ================\n")
print(round(resumen[, c("nivel_ruido", "tiempo_s", "neff_min", "cvmc_medio",
                        "mse", "snr", "psnr", "mse_R", "mse_G", "mse_B",
                        "snr_R", "snr_G", "snr_B", "lambda_media",
                        "acept_mh")], 4), row.names = FALSE)
write.csv(resumen, file.path(dir_out, "3.4_resumen.csv"), row.names = FALSE)

resultados_ligero <- lapply(resultados, function(r) {
  r$res$X_th <- NULL; r
})
saveRDS(list(resumen = resumen, resultados = resultados_ligero),
        file.path(dir_out, "3.4_resultados.rds"))
cat("\nResultados guardados en:", dir_out, "\n")
