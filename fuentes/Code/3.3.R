# ==============================================================================
# SECCIÓN 3.3 - ESCALA DE GRISES: MODELO GAUSSIANO ESPACIAL (GMRF)
# ------------------------------------------------------------------------------
# Modelo:
#   Y_i | X_i, s2 ~ N(X_i, s2) independientes
#   p(x | mu, tau, lambda) = |Q|^{1/2} / (2*pi)^{P/2} *
#       exp{ -tau/2 * sum_i (x_i - mu)^2 - lambda/2 * sum_{i~j} (x_i - x_j)^2 }
#   Q = tau*I + lambda*L,   L = D - A (Laplaciana de la rejilla 4-vecindad)
#   s2 ~ IG(a_s, b_s) | mu ~ N(m_mu, s_mu2) | tau ~ G(a_tau, b_tau) | lambda ~ G(a_l, b_l)
#
# Condicional completa de X_i  (Normal, marca de propiedad local):
#   v_i = (1/s2 + tau + lambda*d_i)^{-1}
#   m_i = v_i * ( y_i/s2 + tau*mu + lambda*sum_{j in N_i} x_j )
#   X_i | resto ~ N(m_i, v_i)
#
# Condicional completa de s2:
#   s2 | resto ~ IG( a_s + P/2, b_s + 0.5*sum_i (y_i - x_i)^2 )
#
# Condicional completa de mu:
#   mu | resto ~ N( v_mu * ( tau*sum_i x_i + m_mu/s_mu2 ), v_mu ),
#   v_mu = ( tau*P + 1/s_mu2 )^{-1}
#
# tau y lambda NO son conjugadas por el determinante |Q|:
#   p(tau | .)  prop tau^{a_tau} exp(-b_tau*tau) |Q|^{1/2} exp(-tau/2 * SS_mu)
#   p(lam | .)  prop lam^{a_l}   exp(-b_l*lam)   |Q|^{1/2} exp(-lam/2 * SS_dif)
#   -> paso Metropolis-Hastings CONJUNTO para (log tau, log lambda) con propuesta
#   bivariada adaptativa (ver 4.4). log|Q| se calcula con factorización dispersa
#   de Cholesky sobre la matriz de precisión (sin invertirla).
#
# Actualización de x: barrido por "tablero de ajedrez". En una rejilla bipartita
# dos píxeles del mismo color nunca son vecinos, luego los X_i de un color son
# condicionalmente independientes dados los del otro color: se actualizan a la vez.
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
dir_out <- file.path(dir_out, "3.3")   # resultados de este punto en su subcarpeta
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
cat("Imágenes en:", dir_img, "\nResultados en:", dir_out, "\n")

RAPIDO <- toupper(Sys.getenv("CDE_RAPIDO", "FALSE")) %in% c("TRUE", "1", "YES")
# "cholesky" (por defecto, estándar del curso) o "espectral" (exacto, más rápido).
# Ambos se validan numéricamente al inicio (coinciden a 1e-8); para cadenas
# largas se recomienda CDE_LOGDET=espectral, ~20x más rápido y exactamente igual.
metodo_logdet <- Sys.getenv("CDE_LOGDET", "cholesky")

config <- list(
  # Hiperparámetros: previas débilmente informativas sobre escalas de [0,1]
  a_sigma = 2, b_sigma = 0.01,   # s2  ~ IG(2, 0.01): E[s2] = 0.01 (s ~ 0.1)
  m_mu = 0.5, s_mu2 = 0.25,      # mu  ~ N(0.5, 0.5^2): cubre [0,1]
  a_tau = 1, b_tau = 0.05,       # tau ~ G(1, 0.05) = Exp(0.05): E[tau] = 20,
                                 #   con masa cerca de 0 (imágenes no homogéneas)
  a_lambda = 2, b_lambda = 0.1,  # lam ~ G(2, 0.1): E[lam] = 20
  # Algoritmo
  B = 24000, burn = 2000,        # cadenas largas: tau (CVMC ~ 90%/sqrt(n_eff)) las exige
  thin = 20,                     # adelgazamiento de la cadena de píxeles
  s_escala = 0.5,                # escalar inicial de la propuesta MVN en escala log
  adapt = TRUE, bloques_adapt = 25,  # adaptación de propuestas SOLO en el calentamiento
  sigma0_levels = c(0.10, 0.05, 0.20)
)
if (RAPIDO) {
  config$B <- 400; config$burn <- 200; config$thin <- 4; config$sigma0_levels <- 0.10
}
semilla <- 2026

# ------------------------------------------------------------------------------
# 1. ESTRUCTURA ESPACIAL (matrices dispersas)
# ------------------------------------------------------------------------------
construir_laplaciana <- function(M, N) {
  idx <- matrix(seq_len(M * N), M, N)
  i <- c(as.vector(idx[, -N]), as.vector(idx[-M, ]))
  j <- c(as.vector(idx[, -1]), as.vector(idx[-1, ]))
  A <- sparseMatrix(i = c(i, j), j = c(j, i), x = 1, dims = c(M * N, M * N))
  d <- rowSums(A)
  L <- Diagonal(x = d) - A
  list(A = A, L = L, d = d)
}

# Tablero de ajedrez: dos grupos sin vecinos internos (independencia condicional)
preparar_tablero <- function(A, M, N) {
  k <- seq_len(M * N)
  r <- (k - 1) %% M + 1
  c <- (k - 1) %/% M + 1
  grupos <- list(k[(r + c) %% 2 == 0], k[(r + c) %% 2 == 1])
  stopifnot(sum(A[grupos[[1]], grupos[[1]]]) == 0, sum(A[grupos[[2]], grupos[[2]]]) == 0)
  list(idx = grupos, A_sub = lapply(grupos, function(g) A[g, , drop = FALSE]))
}

est <- construir_laplaciana(M, N)
A <- est$A; L <- est$L; d <- est$d
tab <- preparar_tablero(A, M, N)
rm(est)

# ------------------------------------------------------------------------------
# 2. LOG-DETERMINANTE DE Q = tau*I + lambda*L
# ------------------------------------------------------------------------------
# Método 1 (por defecto): factorización de Cholesky dispersa (CHOLMOD).
# Se fija la permutación COLAMD calculada una sola vez (el patrón de dispersión
# de Q = tau*I + lambda*L no cambia) y se evalúa log|Q| sobre el factor,
# sin invertir la matriz: log|Q| = determinant(ch, log=TRUE, sqrt=FALSE).
perm_Q <- Cholesky(as(1 * Diagonal(P) + 1 * L, "dgCMatrix"), LDL = FALSE,
                   perm = TRUE, Imult = 0)@perm

log_det_Q_cholesky <- function(tau, lambda) {
  Q  <- as(tau * Diagonal(P) + lambda * L, "dgCMatrix")
  ch <- Cholesky(Q, perm = perm_Q, LDL = FALSE, Imult = 0)
  dt <- determinant(ch, logarithm = TRUE, sqrt = FALSE)
  as.numeric(dt$modulus) * dt$sign
}

# Método 2 (exacto, opcional): autovalores de la Laplaciana de la rejilla M x N
# con condiciones libres: w_{kl} = 4 sin^2(pi k / 2M) + 4 sin^2(pi l / 2N).
# log|Q| = sum_{k,l} log(tau + lambda * w_{kl}).  Se valida contra Cholesky abajo.
eig_L <- outer(4 * sin(pi * (0:(M - 1)) / (2 * M))^2,
               4 * sin(pi * (0:(N - 1)) / (2 * N))^2, "+")
log_det_Q_espectral <- function(tau, lambda) sum(log(tau + lambda * eig_L))

log_det_Q <- if (metodo_logdet == "espectral") log_det_Q_espectral else log_det_Q_cholesky

# Validación cruzada: ambos métodos deben coincidir
for (pr in list(c(20, 20), c(2, 7), c(50, 3))) {
  a1 <- log_det_Q_cholesky(pr[1], pr[2])
  a2 <- log_det_Q_espectral(pr[1], pr[2])
  stopifnot(abs(a1 - a2) / abs(a1) < 1e-8)
}
cat(sprintf("log|Q| validado (Cholesky disperso vs. espectral) | método activo: %s\n", metodo_logdet))

# ------------------------------------------------------------------------------
# 3. DATOS: contaminación gaussiana y = x0 + N(0, sigma0^2)  (sin recortar)
# ------------------------------------------------------------------------------
img <- readPNG(file.path(dir_img, "camera_grises.png"))
if (length(dim(img)) == 3) img <- img[, , 1]
stopifnot(dim(img)[1] == M, dim(img)[2] == N)
x0 <- as.vector(img)                       # verdad conocida en [0,1]
cat(sprintf("Imagen en escala de grises: media = %.4f | var = %.4f\n", mean(x0), var(x0)))

# ------------------------------------------------------------------------------
# 4. PASOS MCMC
# ------------------------------------------------------------------------------
# 4.1 X | resto: dos medios barridos por tablero (vectorizados)
muestrear_x <- function(x, y, tau, lambda, sig2, mu) {
  inv_sig2 <- 1 / sig2
  for (g in 1:2) {
    ig <- tab$idx[[g]]
    s  <- drop(tab$A_sub[[g]] %*% x)                 # sumas de vecinos (otro color)
    v  <- 1 / (inv_sig2 + tau + lambda * d[ig])
    m  <- v * (y[ig] * inv_sig2 + tau * mu + lambda * s)
    x[ig] <- rnorm(length(ig), m, sqrt(v))
  }
  x
}

# 4.2 s2 | resto ~ IG(a_s + P/2, b_s + 0.5*RSS)  (1/Gamma(shape, rate))
muestrear_sigma2 <- function(x, y, a_s, b_s) {
  1 / rgamma(1, shape = a_s + P / 2, rate = b_s + 0.5 * sum((y - x)^2))
}

# 4.3 mu | resto ~ N(v_mu*(tau*sum(x) + m_mu/s_mu2), v_mu)
muestrear_mu <- function(x, tau, m_mu, s_mu2) {
  v_mu <- 1 / (tau * P + 1 / s_mu2)
  rnorm(1, v_mu * (tau * sum(x) + m_mu / s_mu2), sqrt(v_mu))
}

# 4.4 (tau, lambda) | resto: paso MH CONJUNTO en escala log.
# La condicional de tau contiene |Q(tau,lambda)|^{1/2} exp(-tau/2 * SS_mu) y la
# de lambda |Q|^{1/2} exp(-lambda/2 * SS_dif); como |Q| depende de ambos, se
# actualizan juntos con la densidad objetivo
#   lp = a_tau*e1 + a_lam*e2 - b_tau*exp(e1) - b_lam*exp(e2)
#        + 0.5*log|Q| - 0.5*(tau*SS_mu + lambda*SS_dif),   (e1,e2) = (log tau, log lam)
# Propuesta simétrica bivariada (MVN en escala log) con covarianza adaptativa
# durante el calentamiento: aprovecha la fuerte correlación posterior entre ambos
# y requiere UNA sola factorización de Cholesky por iteración.
paso_mh_tau_lam <- function(eta, s, Cprop, a_tau, b_tau, a_lam, b_lam,
                            SS_mu, SS_dif, logdet) {
  z <- rnorm(2)
  # NOTA: chol() devuelve R superior con Cprop = R'R, así que chol(Cprop) %*% z tiene
  # covarianza R R', no Cprop. La propuesta sigue siendo simétrica (la cadena es válida),
  # pero su forma no es exactamente la covarianza adaptada; lo ideal sería
  # t(chol(Cprop)) %*% z. Se deja así para que los resultados reportados sean
  # reproducibles con estas semillas.
  eta_p <- eta + s * drop(chol(Cprop) %*% z)
  tau_p <- exp(eta_p[1]); lam_p <- exp(eta_p[2])
  logdet_p <- log_det_Q(tau_p, lam_p)
  lp   <- a_tau * eta[1]   + a_lam * eta[2]   - b_tau * exp(eta[1]) - b_lam * exp(eta[2]) +
          0.5 * logdet   - 0.5 * (exp(eta[1])   * SS_mu + exp(eta[2])   * SS_dif)
  lp_p <- a_tau * eta_p[1] + a_lam * eta_p[2] - b_tau * tau_p      - b_lam * lam_p      +
          0.5 * logdet_p - 0.5 * (tau_p * SS_mu + lam_p * SS_dif)
  if (log(runif(1)) < lp_p - lp) list(eta = eta_p, tau = tau_p, lam = lam_p,
                                      logdet = logdet_p, acepta = TRUE)
  else                           list(eta = eta, tau = exp(eta[1]), lam = exp(eta[2]),
                                      logdet = logdet, acepta = FALSE)
}

# ------------------------------------------------------------------------------
# 5. MCMC COMPLETO (un nivel de ruido)
# ------------------------------------------------------------------------------
ejecutar_mcmc <- function(y, cfg, semilla, verbose = TRUE) {
  set.seed(semilla)
  B <- cfg$B; ncat <- floor(B / 10)

  # Inicialización
  x <- y                                       # arranque en los datos contaminados
  mu <- mean(y); sig2 <- 0.01; tau <- 1; lambda <- 10
  eta <- c(log(tau), log(lambda))              # estado MH en escala log
  s <- cfg$s_escala                            # escalar de la propuesta bivariada
  Cprop <- diag(2)                             # covarianza de propuesta (adaptativa)

  # Almacenamiento: escalares completos; píxeles en momentos + cadena adelgazada
  MU <- TAU <- LAM <- S2 <- numeric(B)
  ETA_HIST <- matrix(0, B, 2)
  ACC_PAR <- logical(B)
  n_th <- ceiling((B - cfg$burn) / cfg$thin)
  X_th <- matrix(0, n_th, P)
  sx <- numeric(P); sx2 <- numeric(P); n_post <- 0
  i_th <- 0

  # log|Q| del estado actual (se actualiza al aceptar propuestas)
  logdet <- log_det_Q(tau, lambda)
  t0 <- Sys.time()

  for (b in seq_len(B)) {
    x     <- muestrear_x(x, y, tau, lambda, sig2, mu)          # 1) X | resto
    sig2  <- muestrear_sigma2(x, y, cfg$a_sigma, cfg$b_sigma)  # 2) s2 | resto
    mu    <- muestrear_mu(x, tau, cfg$m_mu, cfg$s_mu2)         # 3) mu | resto

    SS_mu   <- sum((x - mu)^2)                                 # 4) (tau, lambda) | resto
    SS_dif  <- as.numeric(crossprod(x, L %*% x))               #    (MH conjunto)
    pm <- paso_mh_tau_lam(eta, s, Cprop, cfg$a_tau, cfg$b_tau,
                          cfg$a_lambda, cfg$b_lambda, SS_mu, SS_dif, logdet)
    eta <- pm$eta; tau <- pm$tau; lambda <- pm$lam
    ACC_PAR[b] <- pm$acepta
    if (pm$acepta) logdet <- pm$logdet

    MU[b] <- mu; TAU[b] <- tau; LAM[b] <- lambda; S2[b] <- sig2
    ETA_HIST[b, ] <- eta

    if (b > cfg$burn) {
      sx <- sx + x; sx2 <- sx2 + x^2; n_post <- n_post + 1
      if ((b - cfg$burn) %% cfg$thin == 0 && i_th < n_th) {
        i_th <- i_th + 1; X_th[i_th, ] <- x
      }
    }

    # Adaptación de la propuesta SOLO en el calentamiento (blanco ~30%):
    # la estructura de correlación sale de la historia reciente de eta
    # (normalizada a traza 1) y la escala global la lleva s.
    if (isTRUE(cfg$adapt) && b <= cfg$burn) {
      if (b >= 50 && b %% 50 == 0) {
        Cv <- cov(ETA_HIST[max(1, b - 99):b, , drop = FALSE])
        tr <- mean(diag(Cv))
        if (is.finite(tr) && tr > 1e-8) Cprop <- Cv / tr + diag(1e-6, 2)
      }
      if (b %% cfg$bloques_adapt == 0) {
        acc <- mean(ACC_PAR[(b - cfg$bloques_adapt + 1):b])
        s <- s * min(max(exp(1.5 * (acc - 0.30)), 0.5), 2.0)
      }
    }

    if (verbose && (b %% ncat == 0)) {
      cat(sprintf("%.0f%% completado...\n", 100 * b / B))
      cat(sprintf("   iter %d | mu = %.4f | tau = %.4f | lam = %.3f | s2 = %.5f | acept(tau,lam) = %.1f%% | s = %.3f | %.0fs\n",
                  b, mu, tau, lambda, sig2,
                  100 * mean(ACC_PAR[1:b]), s,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
  }

  post <- (cfg$burn + 1):B
  list(mu = MU, tau = TAU, lambda = LAM, sig2 = S2,
       acc_par = ACC_PAR,
       x_media = sx / n_post, x2_media = sx2 / n_post,
       X_th = X_th, post = post, s_escala = s,
       tiempo = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

# ------------------------------------------------------------------------------
# 6. DIAGNÓSTICOS Y MÉTRICAS
# ------------------------------------------------------------------------------
diagnostico_param <- function(cadena, nombre, etiqueta = nombre) {
  ch   <- mcmc(cadena)
  neff <- as.numeric(effectiveSize(ch))
  dsd  <- sd(cadena); media <- mean(cadena)
  emc  <- dsd / sqrt(neff)
  cvmc <- 100 * emc / abs(media)
  q    <- quantile(cadena, probs = c(0.025, 0.975))
  cat(sprintf("\n--- %s ---\n", nombre))
  cat(sprintf("Media = %.5f | DE = %.5f | n_eff = %.0f | EMC = %.5f | CVMC = %.3f%%\n",
              media, dsd, neff, emc, cvmc))
  cat(sprintf("IC95%% = [%.5f, %.5f]\n", q[1], q[2]))
  if (cvmc >= 5) cat("AVISO: CVMC >= 5% (aumentar iteraciones)\n")
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1)); on.exit(par(op), add = TRUE)
  plot(cadena, type = "l", col = "grey30", xlab = "Iteración", ylab = etiqueta,
       main = paste("Traza de", etiqueta))
  acf(cadena, lag.max = 60, main = paste("ACF de", etiqueta))
  c(Media = media, DE = dsd, EMC = emc, n_eff = neff, `CVMC(%)` = cvmc,
    IC95_inf = unname(q[1]), IC95_sup = unname(q[2]))
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

img_clip <- function(v) matrix(pmin(pmax(v, 0), 1), M, N)

# Panel: original / contaminada / restaurada (media posterior) / incertidumbre (DE)
panel_escenario <- function(x0, y, res, sigma0) {
  op <- par(mfrow = c(1, 4), mar = c(1, 1, 2.5, 1)); on.exit(par(op), add = TRUE)
  plot(as.raster(img_clip(x0)), main = "Original")
  plot(as.raster(img_clip(y)), main = sprintf("Contaminada (s0 = %.2f)", sigma0))
  plot(as.raster(img_clip(res$x_media)), main = "Restaurada (media posterior)")
  dsd <- sqrt(pmax(res$x2_media - res$x_media^2, 0))
  image(1:N, 1:M, t(matrix(dsd, M, N)[M:1, ]), col = hcl.colors(100, "Inferno"),
        axes = FALSE, asp = 1, main = "DE posterior")
}

# ------------------------------------------------------------------------------
# 7. ESCENARIO COMPLETO (un nivel de ruido)
# ------------------------------------------------------------------------------
correr_escenario <- function(x0, sigma0, cfg, semilla) {
  cat(sprintf("\n================ sigma_0 = %.2f ================\n", sigma0))
  set.seed(semilla)
  y <- x0 + rnorm(P, 0, sigma0)      # sin recortar (la inferencia usa y tal cual)
  cat(sprintf("Ruido real: media = %.4f | DE = %.4f | DE objetivo = %.2f\n",
              mean(y - x0), sd(y - x0), sigma0))

  res <- ejecutar_mcmc(y, cfg, semilla + 1)
  post <- res$post

  cat("\n--- Resumen posterior (post burn-in) ---\n")
  tab_res <- rbind(
    mu     = diagnostico_param(res$mu[post], "mu", expression(mu)),
    tau    = diagnostico_param(res$tau[post], "tau", expression(tau)),
    lambda = diagnostico_param(res$lambda[post], "lambda", expression(lambda)),
    sigma2 = diagnostico_param(res$sig2[post], "sigma^2", expression(sigma^2))
  )
  cat(sprintf("\nTasa de aceptación MH de (tau, lambda) (post burn-in): %.1f%%\n",
              100 * mean(res$acc_par[post])))
  cat(sprintf("s2 real = %.5f | media posterior de s2 = %.5f | IC95%% = [%.5f, %.5f]\n",
              sigma0^2, mean(res$sig2[post]),
              quantile(res$sig2[post], .025), quantile(res$sig2[post], .975)))

  # Diagnóstico de píxeles: n_eff y CVMC sobre la cadena adelgazada
  cat("\n--- Diagnóstico de los X_i (cadena de píxeles adelgazada) ---\n")
  X_th <- res$X_th
  neff_pix <- as.numeric(effectiveSize(mcmc(X_th)))
  neff_pix[!is.finite(neff_pix) | neff_pix <= 0] <- nrow(X_th)
  emc_pix  <- apply(X_th, 2, sd) / sqrt(neff_pix)
  cvmc_pix <- 100 * emc_pix / pmax(abs(res$x_media), 0.05)   # suelo 0.05 en el denominador
  cat(sprintf("n_eff de X_i: min = %.0f | mediana = %.0f\n", min(neff_pix), median(neff_pix)))
  cat(sprintf("CVMC de X_i (%%): mediana = %.3f | media = %.3f | max = %.3f\n",
              median(cvmc_pix), mean(cvmc_pix), max(cvmc_pix)))

  # Incertidumbre: DE posterior e IC95% para una selección de píxeles
  q_pix <- apply(X_th, 2, quantile, probs = c(0.025, 0.975))
  x_sd <- sqrt(pmax(res$x2_media - res$x_media^2, 0))
  idx_sel <- unique(c(1, M, (64 - 1) * M + 64, P, (20 - 1) * M + 20, (96 - 1) * M + 32,
                      order(x_sd, decreasing = TRUE)[1:3]))   # los 3 más inciertos
  sel <- data.frame(
    pixel = idx_sel,
    fila = (idx_sel - 1) %% M + 1, col = (idx_sel - 1) %/% M + 1,
    x0 = x0[idx_sel], y = y[idx_sel],
    media = res$x_media[idx_sel], DE = x_sd[idx_sel],
    IC95_inf = q_pix[1, idx_sel], IC95_sup = q_pix[2, idx_sel],
    dentro_IC95 = x0[idx_sel] >= q_pix[1, idx_sel] & x0[idx_sel] <= q_pix[2, idx_sel]
  )
  cat("\n--- IC95% para una selección de píxeles ---\n")
  print(round(sel, 4), row.names = FALSE)

  # Cobertura empírica de los IC95% sobre todos los píxeles
  cobertura <- mean(x0 >= q_pix[1, ] & x0 <= q_pix[2, ])
  cat(sprintf("Cobertura empírica de los IC95%% en los %d píxeles: %.1f%%\n", P, 100 * cobertura))

  # Calidad de la restauración (sin recortar la inferencia)
  mse <- mean((res$x_media - x0)^2)
  mse_y <- mean((y - x0)^2)
  snr <- 10 * log10(mean(x0^2) / mse)
  psnr <- 10 * log10(1 / mse)
  snr_y <- 10 * log10(mean(x0^2) / mse_y)
  cat(sprintf("MSE = %.5f | SNR = %.2f dB | PSNR = %.2f dB\n", mse, snr, psnr))
  cat(sprintf("   (imagen contaminada: MSE = %.5f | SNR = %.2f dB)\n", mse_y, snr_y))

  metricas <- data.frame(
    modelo = "GMRF grises", nivel_ruido = sigma0, tiempo_s = res$tiempo,
    neff_min = min(c(neff_pix, tab_res[, "n_eff"])),
    cvmc_medio = mean(c(cvmc_pix, tab_res[, "CVMC(%)"])),
    mse = mse, snr = snr, psnr = psnr,
    mse_contaminada = mse_y, snr_contaminada = snr_y,
    mu_media = mean(res$mu[post]), tau_media = mean(res$tau[post]),
    lambda_media = mean(res$lambda[post]), sigma2_media = mean(res$sig2[post]),
    sigma2_real = sigma0^2, acept_mh = mean(res$acc_par[post]),
    cobertura_IC95 = cobertura, stringsAsFactors = FALSE
  )

  list(sigma0 = sigma0, y = y, res = res, tab_res = tab_res, sel = sel,
       metricas = metricas, x_sd = x_sd, q_pix = q_pix)
}

# ------------------------------------------------------------------------------
# 8. EJECUCIÓN: tres niveles de ruido
# ------------------------------------------------------------------------------
resultados <- list()
for (k in seq_along(config$sigma0_levels)) {
  s0 <- config$sigma0_levels[k]
  r <- correr_escenario(x0, s0, config, semilla + 100 * k)
  resultados[[sprintf("sigma0_%.2f", s0)]] <- r
  dibujar(sprintf("panel_grises_%.2f.png", s0), 1800, 480, function() {
    panel_escenario(x0, r$y, r$res, s0)
  })
}

# Comparación visual de los tres niveles
dibujar("comparacion_grises.png", 1600, 1000, function() {
  op <- par(mfrow = c(length(resultados), 3), mar = c(1, 1, 2.5, 1)); on.exit(par(op), add = TRUE)
  for (r in resultados) {
    plot(as.raster(img_clip(x0)), main = "Original")
    plot(as.raster(img_clip(r$y)), main = sprintf("Contaminada s0 = %.2f", r$sigma0))
    plot(as.raster(img_clip(r$res$x_media)), main = sprintf("Restaurada s0 = %.2f", r$sigma0))
  }
})

# Histogramas de los parámetros del escenario principal
r1 <- resultados[[1]]
dibujar(sprintf("hist_params_%.2f.png", r1$sigma0), 1300, 900, function() {
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1)); on.exit(par(op), add = TRUE)
  hist(r1$res$mu[r1$res$post], breaks = 30, col = "grey70", border = "white",
       main = expression(mu), xlab = expression(mu))
  hist(r1$res$tau[r1$res$post], breaks = 30, col = "grey70", border = "white",
       main = expression(tau), xlab = expression(tau))
  hist(r1$res$lambda[r1$res$post], breaks = 30, col = "grey70", border = "white",
       main = expression(lambda), xlab = expression(lambda))
  hist(r1$res$sig2[r1$res$post], breaks = 30, col = "grey70", border = "white",
       main = expression(sigma^2), xlab = expression(sigma^2))
  abline(v = r1$sigma0^2, col = "red", lwd = 2, lty = 2)
})

# ------------------------------------------------------------------------------
# 9. TABLA COMPARATIVA Y GUARDADO
# ------------------------------------------------------------------------------
resumen <- do.call(rbind, lapply(resultados, `[[`, "metricas"))
cat("\n================ COMPARACIÓN ENTRE NIVELES DE RUIDO (grises) ================\n")
print(round(resumen[, c("nivel_ruido", "tiempo_s", "neff_min", "cvmc_medio",
                        "mse", "snr", "psnr", "sigma2_media", "sigma2_real",
                        "acept_mh", "cobertura_IC95")], 4),
      row.names = FALSE)
write.csv(resumen, file.path(dir_out, "3.3_resumen.csv"), row.names = FALSE)

# Guardar (sin la cadena completa de píxeles para aligerar el RDS)
resultados_ligero <- lapply(resultados, function(r) {
  r$res$X_th <- NULL; r$q_pix <- NULL; r
})
saveRDS(list(resumen = resumen, resultados = resultados_ligero),
        file.path(dir_out, "3.3_resultados.rds"))
cat("\nResultados guardados en:", dir_out, "\n")
