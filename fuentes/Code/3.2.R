# ==============================================================================
# SECCIÓN 3.2 - IMAGEN BINARIA: MODELO DE ISING
# Restauración bayesiana: Gibbs (X, rho) + algoritmo de intercambio aproximado (beta)
# ------------------------------------------------------------------------------
# Modelo:  y_i | x_i, rho ~ Bernoulli-flip(rho)      (y_i = -x_i con prob. rho)
#          p(x | beta) = exp{beta * S(x)} / Z(beta),  S(x) = sum_{i~j} x_i x_j
#          beta ~ Gamma(a_b, b_b) (tasa),  rho ~ Beta(a_r, b_r) truncada en (0, 1/2)
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
dir_out <- file.path(dir_out, "3.2")   # resultados de este punto en su subcarpeta
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
cat("Imágenes en:", dir_img, "\nResultados en:", dir_out, "\n")

# Modo de prueba (CDE_RAPIDO=TRUE): cadenas muy cortas, solo para validar el código
RAPIDO <- toupper(Sys.getenv("CDE_RAPIDO", "FALSE")) %in% c("TRUE", "1", "YES")

config <- list(
  # Hiperparámetros (previas débilmente informativas)
  a_beta = 2, b_beta = 4,      # beta ~ Gamma(2, 4): media 0.5, sd 0.35
  a_rho  = 1, b_rho  = 1,      # rho  ~ Beta(1, 1) truncada en (0, 0.5) = Unif(0, 0.5)
  # Algoritmo
  B      = 4000,               # iteraciones totales
  burn   = 1000,               # calentamiento
  delta  = 0.02,               # desv. estándar inicial de la propuesta de beta
  K      = 25,                 # barridos de Gibbs para X*. OJO: la media de beta baja al subir K (sesgo); verifica con K=50
  beta0  = 0.5,                # valor inicial de beta
  rho_ini = 0.10,              # valor inicial de rho
  # Adaptación de delta durante el calentamiento (blanco de aceptación ~30%)
  adapt  = TRUE, bloques_adapt = 50,
  # Cadena de píxeles adelgazada para diagnósticos (ahorro de memoria)
  thin   = 10
)
if (RAPIDO) {
  config$B <- 120; config$burn <- 40; config$K <- 5; config$thin <- 4
}
niveles_ruido <- if (RAPIDO) 0.10 else c(0.10, 0.05, 0.20)   # el principal primero (rho_0 = 0.10)
semilla <- 2026

# ------------------------------------------------------------------------------
# 1. ESTRUCTURA ESPACIAL DISPERSA
# ------------------------------------------------------------------------------
# Matriz de adyacencia A (dgCMatrix, 4 vecinos). Orden por columnas: k = (c-1)*M + r
construir_adyacencia <- function(M, N) {
  idx <- matrix(seq_len(M * N), M, N)
  i <- c(as.vector(idx[, -N]), as.vector(idx[-M, ]))
  j <- c(as.vector(idx[, -1]), as.vector(idx[-1, ]))
  A <- sparseMatrix(i = c(i, j), j = c(j, i), x = 1, dims = c(M * N, M * N))
  as(A, "generalMatrix")
}

# Concordancia S(x) = sum_{i~j} x_i x_j = (1/2) x' A x   (cada arista una sola vez)
S_concordancia <- function(x, A) as.numeric(crossprod(x, A %*% x)) / 2

# Partición "tablero de ajedrez": los píxeles de un mismo color no son vecinos
# entre sí, de modo que son condicionalmente independientes y se actualizan a la vez.
preparar_tablero <- function(A, M, N) {
  k <- seq_len(M * N)
  r <- (k - 1) %% M + 1
  c <- (k - 1) %/% M + 1
  grupos <- list(k[(r + c) %% 2 == 0], k[(r + c) %% 2 == 1])
  stopifnot(sum(A[grupos[[1]], grupos[[1]]]) == 0, sum(A[grupos[[2]], grupos[[2]]]) == 0)
  list(idx = grupos, A_sub = lapply(grupos, function(g) A[g, , drop = FALSE]))
}

# ------------------------------------------------------------------------------
# 2. CONTAMINACIÓN: y_i = -x_i con probabilidad rho_0
# ------------------------------------------------------------------------------
contaminar <- function(x0, rho0) {
  flip <- runif(length(x0)) < rho0
  y <- ifelse(flip, -x0, x0)
  y
}

# ------------------------------------------------------------------------------
# 3. MUESTREADORES
# ------------------------------------------------------------------------------
# Barrido de Gibbs vectorizado (2 medios barridos por tablero).
# Condicional completa: P(X_i=+1 | resto, y) = plogis( 2*beta*sum_{j in N_i} x_j + y_i*log((1-rho)/rho) )
# Si y = NULL se muestrea del Ising "puro" (imagen auxiliar, sin usar Y).
barrido_gibbs <- function(x, beta, tab, y = NULL, rho = NULL) {
  lrat <- if (!is.null(y)) log((1 - rho) / rho) else 0
  for (g in 1:2) {
    ig  <- tab$idx[[g]]
    eta <- 2 * beta * drop(tab$A_sub[[g]] %*% x)
    if (!is.null(y)) eta <- eta + lrat * y[ig]
    x[ig] <- 2 * (runif(length(ig)) < plogis(eta)) - 1
  }
  x
}

# rho | x, y ~ Beta(a_r + d_y, b_r + P - d_y) truncada en (0, 1/2), por inversión de la CDF
muestrear_rho <- function(x, y, a_r, b_r) {
  d_y <- sum(x != y)
  a_p <- a_r + d_y
  b_p <- b_r + length(x) - d_y
  u   <- runif(1) * pbeta(0.5, a_p, b_p)
  qbeta(u, a_p, b_p)
}

# beta por intercambio aproximado (Murray et al., 2006; K barridos de Gibbs para X*)
paso_beta <- function(beta, x, S_x, A, tab, cfg) {
  beta_p <- beta + rnorm(1, 0, cfg$delta)
  if (beta_p <= 0) return(list(beta = beta, S_x = S_x, acepta = FALSE))
  x_aux <- x                                   # arranque en caliente desde x
  for (k in seq_len(cfg$K)) x_aux <- barrido_gibbs(x_aux, beta_p, tab)
  S_aux <- S_concordancia(x_aux, A)
  log_alpha <- dgamma(beta_p, cfg$a_beta, cfg$b_beta, log = TRUE) -
    dgamma(beta,   cfg$a_beta, cfg$b_beta, log = TRUE) +
    (beta_p - beta) * (S_x - S_aux)
  if (log(runif(1)) < log_alpha) list(beta = beta_p, S_x = S_x, acepta = TRUE)
  else                           list(beta = beta,   S_x = S_x, acepta = FALSE)
}

# ------------------------------------------------------------------------------
# 4. MCMC COMPLETO
# ------------------------------------------------------------------------------
ejecutar_mcmc <- function(y, A, tab, cfg, beta0 = cfg$beta0, rho_ini = cfg$rho_ini,
                          x_ini = y, verbose = TRUE) {
  B <- cfg$B
  ncat <- floor(B / 10)                       # avance cada 10% (estándar del curso)
  x <- x_ini; beta <- beta0; rho <- rho_ini
  beta_ch <- rho_ch <- S_ch <- numeric(B)
  acepta  <- logical(B)
  conteo_mas <- numeric(length(y))             # para pi_i = P(X_i = +1 | y)
  n_post <- 0

  # Cadena de píxeles adelgazada (diagnósticos de pi_i con poco memoria)
  n_th <- ceiling((B - cfg$burn) / cfg$thin)
  X_th <- matrix(as.raw(0), n_th, length(y))
  i_th <- 0

  # delta de la propuesta de beta se adapta SOLO durante el calentamiento
  cfg_l <- cfg; delta <- cfg$delta
  t0 <- Sys.time()

  for (b in seq_len(B)) {
    x    <- barrido_gibbs(x, beta, tab, y = y, rho = rho)        # 1) X | beta, rho, y
    rho  <- muestrear_rho(x, y, cfg$a_rho, cfg$b_rho)            # 2) rho | x, y
    S_x  <- S_concordancia(x, A)
    pb   <- paso_beta(beta, x, S_x, A, tab, cfg_l)               # 3) beta | x (intercambio)
    beta <- pb$beta
    beta_ch[b] <- beta; rho_ch[b] <- rho; S_ch[b] <- S_x; acepta[b] <- pb$acepta

    if (b > cfg$burn) {
      conteo_mas <- conteo_mas + (x == 1); n_post <- n_post + 1
      if ((b - cfg$burn) %% cfg$thin == 0 && i_th < n_th) {
        i_th <- i_th + 1
        X_th[i_th, ] <- as.raw(x == 1)         # 1 si x_i = +1
      }
    }

    # Adaptación de delta en el calentamiento (blanco: aceptación ~30%)
    if (isTRUE(cfg$adapt) && b <= cfg$burn && b %% cfg$bloques_adapt == 0) {
      tasa_bloque <- mean(acepta[(b - cfg$bloques_adapt + 1):b])
      delta <- if (tasa_bloque > 0.35) delta * 1.20 else if (tasa_bloque < 0.15) delta * 0.85 else delta
      cfg_l$delta <- delta
    }

    if (verbose && (b %% ncat == 0)) {
      cat(sprintf("%.0f%% completado...\n", 100 * b / B))
      cat(sprintf("   iter %d | beta = %.4f | rho = %.4f | S/aristas = %.4f | acept(beta) = %.1f%% | delta = %.4f | %.0fs\n",
                  b, beta, rho, S_x / (M * (N - 1) + N * (M - 1)),
                  100 * mean(acepta[1:b]), delta,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
  }
  list(beta = beta_ch, rho = rho_ch, S = S_ch, acepta = acepta,
       pi = conteo_mas / n_post, x_final = x, delta = delta,
       X_th = X_th, tiempo = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

# ------------------------------------------------------------------------------
# 5. DIAGNÓSTICOS Y RESÚMENES
# ------------------------------------------------------------------------------
resumen_posterior <- function(cadena) {
  ch    <- mcmc(cadena)
  media <- mean(cadena); dsd <- sd(cadena)
  neff  <- as.numeric(effectiveSize(ch))
  emc   <- dsd / sqrt(neff)
  q     <- quantile(cadena, probs = c(0.025, 0.975))
  c(Media = media, DE = dsd, EMC = emc, n_eff = neff, `CVMC(%)` = 100 * emc / abs(media),
    IC95_inf = unname(q[1]), IC95_sup = unname(q[2]),
    Geweke_z = as.numeric(geweke.diag(ch)$z))
}

# Diagnóstico de un parámetro escalar (estándar del curso):
# traza, ACF, n_eff, EMC, CVMC (%) e IC95%. Imprime en consola y guarda PNG.
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
  invisible(c(Media = media, DE = dsd, EMC = emc, n_eff = neff, `CVMC(%)` = cvmc,
              IC95_inf = q[1], IC95_sup = q[2]))
}

# Proporción de píxeles recuperados, global y por tipo (borde / región homogénea)
precision_restauracion <- function(x_map, x0, y, A) {
  s0     <- drop(A %*% x0)
  homog  <- (s0 * x0) == rowSums(A)                       # todos los vecinos iguales al píxel
  c(MAP_global = mean(x_map == x0), MAP_bordes = mean((x_map == x0)[!homog]),
    MAP_homogeneas = mean((x_map == x0)[homog]), Y_global = mean(y == x0),
    prop_bordes = mean(!homog))
}

# ------------------------------------------------------------------------------
# 6. ESCENARIO COMPLETO (un nivel de ruido)
# ------------------------------------------------------------------------------
correr_escenario <- function(x0, rho0, A, tab, cfg, semilla) {
  set.seed(semilla)
  cat(sprintf("\n================ rho_0 = %.2f ================\n", rho0))
  y <- contaminar(x0, rho0)
  cat(sprintf("Píxeles invertidos: %d (%.2f%%)\n", sum(y != x0), 100 * mean(y != x0)))

  res <- ejecutar_mcmc(y, A, tab, cfg)
  post <- (cfg$burn + 1):cfg$B
  beta_p <- res$beta[post]; rho_p <- res$rho[post]

  tab_res <- rbind(beta = resumen_posterior(beta_p), rho = resumen_posterior(rho_p))
  cat("\n--- Resumen posterior (post burn-in) ---\n"); print(round(tab_res, 4))
  tasa <- mean(res$acepta[post])
  cat(sprintf("Tasa de aceptación de beta (post burn-in): %.1f%% (delta final = %.4f)\n",
              100 * tasa, res$delta))

  # Diagnósticos de cada parámetro escalar: traza, ACF, n_eff, EMC, CVMC, IC95%
  dibujar(sprintf("diag_beta_%.2f.png", rho0), 1200, 450,
          function() diagnostico_param(beta_p, "beta", expression(beta)))
  dibujar(sprintf("diag_rho_%.2f.png", rho0), 1200, 450,
          function() diagnostico_param(rho_p, "rho", expression(rho)))

  # Diagnóstico de pi_i = P(X_i = +1 | y) con las cadenas de píxeles adelgazadas
  cat("\n--- Diagnóstico de pi_i (cadenas de píxeles adelgazadas) ---\n")
  X_th <- matrix(as.numeric(res$X_th), nrow = nrow(res$X_th))   # 0/1
  neff_pix <- as.numeric(effectiveSize(mcmc(X_th)))
  # Cadenas constantes (pi_i en {0,1}): n_eff indefinido -> se fija en el largo
  neff_pix[!is.finite(neff_pix) | neff_pix <= 0] <- nrow(X_th)
  pi_th    <- colMeans(X_th)
  emc_pix  <- sqrt(pmax(pi_th * (1 - pi_th), 0) / neff_pix)
  cvmc_pix <- 100 * emc_pix / pmax(pi_th, 0.05)   # suelo de 0.05 en el denominador
  cat(sprintf("n_eff de pi_i: min = %.0f | mediana = %.0f\n", min(neff_pix), median(neff_pix)))
  cat(sprintf("CVMC de pi_i (%%): mediana = %.3f | media = %.3f | max = %.3f\n",
              median(cvmc_pix), mean(cvmc_pix), max(cvmc_pix)))

  x_map <- ifelse(res$pi >= 0.5, 1, -1)
  prec  <- precision_restauracion(x_map, x0, y, A)
  cat(sprintf("Píxeles correctos: imagen contaminada = %.2f%% | restaurada (MAP) = %.2f%%\n",
              100 * prec["Y_global"], 100 * prec["MAP_global"]))
  cat(sprintf("  en regiones homogéneas = %.2f%% | en bordes = %.2f%% (bordes = %.1f%% de los píxeles)\n",
              100 * prec["MAP_homogeneas"], 100 * prec["MAP_bordes"], 100 * prec["prop_bordes"]))
  cat(sprintf("rho_0 real = %.3f | IC95%% posterior de rho = [%.4f, %.4f]\n", rho0,
              quantile(rho_p, .025), quantile(rho_p, .975)))

  # Métricas de calidad en escala de imagen {0,1} (0 = negro, 1 = blanco)
  x0_img   <- ifelse(x0 == 1, 0, 1)
  y_img    <- ifelse(y == 1, 0, 1)
  xhat_img <- 1 - res$pi                       # media posterior de Pr(X_i = -1 | y)
  mse  <- mean((xhat_img - x0_img)^2)
  mse_y <- mean((y_img - x0_img)^2)
  snr  <- 10 * log10(mean(x0_img^2) / mse)     # dB
  psnr <- 10 * log10(1 / mse)                  # rango [0,1]
  snr_y <- 10 * log10(mean(x0_img^2) / mse_y)

  metricas <- data.frame(
    modelo = "Ising", nivel_ruido = rho0, tiempo_s = res$tiempo,
    neff_min = min(c(neff_pix, tab_res[, "n_eff"])),
    cvmc_medio = mean(c(cvmc_pix, tab_res[, "CVMC(%)"])),
    mse = mse, snr = snr, psnr = psnr, precision_MAP = unname(prec["MAP_global"]),
    precision_contaminada = unname(prec["Y_global"]),
    mse_contaminada = mse_y, snr_contaminada = snr_y,
    beta_media = mean(beta_p), rho_media = mean(rho_p),
    tasa_acept_beta = tasa, stringsAsFactors = FALSE
  )
  cat(sprintf("MSE = %.5f | SNR = %.2f dB | PSNR = %.2f dB\n", mse, snr, psnr))
  cat(sprintf("   (imagen contaminada: MSE = %.5f | SNR = %.2f dB)\n", mse_y, snr_y))

  c(list(rho0 = rho0, x0 = x0, y = y, x_map = x_map, tab = tab_res, tasa = tasa,
         prec = prec, post = post, metricas = metricas, xhat_img = xhat_img), res)
}

# ------------------------------------------------------------------------------
# 7. VISUALIZACIÓN
# ------------------------------------------------------------------------------
img_ising <- function(x) as.raster(matrix(ifelse(x == 1, 0, 1), M, N))   # +1 negro

mapa_prob <- function(pi_vec, titulo) {
  z <- matrix(pi_vec, M, N)
  image(1:N, 1:M, t(z[M:1, ]), zlim = c(0, 1), col = hcl.colors(100, "YlGnBu", rev = FALSE),
        axes = FALSE, asp = 1, xlab = "", ylab = "", main = titulo, useRaster = TRUE)
}

panel_escenario <- function(r) {
  op <- par(mfrow = c(2, 4), mar = c(4, 4, 3, 1)); on.exit(par(op))
  b <- r$beta[r$post]; rh <- r$rho[r$post]
  plot(r$beta, type = "l", col = "grey40", xlab = "Iteración", ylab = expression(beta),
       main = expression("Traza de " * beta)); abline(v = r$post[1], lty = 3)
  abline(h = mean(b), col = "red", lwd = 2); legend("bottomright", "media post.", col = "red", lwd = 2, bty = "n", cex = .8)
  plot(r$rho, type = "l", col = "grey40", xlab = "Iteración", ylab = expression(rho),
       main = expression("Traza de " * rho)); abline(v = r$post[1], lty = 3)
  abline(h = mean(rh), col = "red", lwd = 2); abline(h = r$rho0, col = "blue", lwd = 2, lty = 2)
  legend("topright", c("media post.", expression(rho[0])), col = c("red", "blue"), lwd = 2, lty = c(1, 2), bty = "n", cex = .8)
  acf(b,  lag.max = 60, main = expression("ACF de " * beta))
  acf(rh, lag.max = 60, main = expression("ACF de " * rho))
  plot(img_ising(r$x0),    main = "Original");
  plot(img_ising(r$y),     main = sprintf("Contaminada (rho0 = %.2f)", r$rho0))
  plot(img_ising(r$x_map), main = sprintf("Restaurada MAP (%.1f%%)", 100 * r$prec["MAP_global"]))
  mapa_prob(r$pi, expression("Mapa de " * pi[i] * " = P(X"[i] * "=+1 | y)"))
}

panel_comparacion <- function(res_list) {
  op <- par(mfrow = c(length(res_list), 4), mar = c(1, 1, 3, 1)); on.exit(par(op))
  for (r in res_list) {
    plot(img_ising(r$x0), main = "Original")
    plot(img_ising(r$y), main = sprintf("Contaminada, rho0 = %.2f", r$rho0))
    plot(img_ising(r$x_map), main = sprintf("MAP: %.1f%% correcto", 100 * r$prec["MAP_global"]))
    mapa_prob(r$pi, expression(pi[i]))
  }
}

# Guarda en PNG y, si la sesión es interactiva, también abre una ventana propia
dibujar <- function(archivo, ancho, alto, fun) {
  png(file.path(dir_out, archivo), width = ancho, height = alto, res = 110); fun(); dev.off()
  if (interactive()) { dev.new(width = ancho / 110, height = alto / 110, noRStudioGD = TRUE); fun() }
}

# ==============================================================================
# 8. EJECUCIÓN
# ==============================================================================
A   <- construir_adyacencia(M, N)
tab <- preparar_tablero(A, M, N)
stopifnot(class(A) == "dgCMatrix")

# Imagen original binaria -> Ising (+1 negro, -1 blanco)
img_h <- readPNG(file.path(dir_img, "horse_binaria.png"))
if (length(dim(img_h)) == 3) img_h <- img_h[, , 1]
x0 <- ifelse(as.vector(img_h) < 0.5, 1, -1)

n_aristas <- M * (N - 1) + N * (M - 1)
cat(sprintf("Imagen original: %d negros (+1), %d blancos (-1) | S(x0) = %d (%.4f por arista)\n",
            sum(x0 == 1), sum(x0 == -1), as.integer(S_concordancia(x0, A)),
            S_concordancia(x0, A) / n_aristas))
cat(sprintf("Beta crítico del Ising 2D en retícula cuadrada: %.4f\n", log(1 + sqrt(2)) / 2))

resultados <- list()
for (k in seq_along(niveles_ruido)) {
  r <- correr_escenario(x0, niveles_ruido[k], A, tab, config, semilla + k)
  resultados[[sprintf("rho_%.2f", niveles_ruido[k])]] <- r
  dibujar(sprintf("panel_rho_%.2f.png", niveles_ruido[k]), 1700, 800, function() panel_escenario(r))
}
dibujar("comparacion_ruidos.png", 1500, 1200, function() panel_comparacion(resultados))

# Tabla comparativa entre niveles de ruido
comp <- do.call(rbind, lapply(resultados, function(r) data.frame(
  rho0 = r$rho0, rho_post = r$tab["rho", "Media"], beta_post = r$tab["beta", "Media"],
  beta_DE = r$tab["beta", "DE"], acept_beta = 100 * r$tasa,
  neff_beta = r$tab["beta", "n_eff"], CVMC_beta = r$tab["beta", "CVMC(%)"],
  neff_rho = r$tab["rho", "n_eff"], CVMC_rho = r$tab["rho", "CVMC(%)"],
  pct_contaminada = 100 * r$prec["Y_global"], pct_MAP = 100 * r$prec["MAP_global"],
  pct_homogeneas = 100 * r$prec["MAP_homogeneas"], pct_bordes = 100 * r$prec["MAP_bordes"])))
cat("\n================ COMPARACIÓN ENTRE NIVELES DE RUIDO ================\n")
print(round(comp, 3), row.names = FALSE)
write.csv(comp, file.path(dir_out, "comparacion_ruidos.csv"), row.names = FALSE)

# Métricas unificadas (las consume 4.R para la tabla comparativa global)
resumen <- do.call(rbind, lapply(resultados, `[[`, "metricas"))
cat("\n--- Métricas unificadas (para la tabla comparativa global) ---\n")
print(round(resumen[, c("nivel_ruido", "tiempo_s", "neff_min", "cvmc_medio",
                        "mse", "snr", "precision_MAP")], 4), row.names = FALSE)
write.csv(resumen, file.path(dir_out, "3.2_resumen.csv"), row.names = FALSE)

# Guardar resultados (sin las cadenas de píxeles, para aligerar el RDS)
resultados_ligero <- lapply(resultados, function(r) { r$X_th <- NULL; r })
saveRDS(list(resumen = resumen, comp = comp, resultados = resultados_ligero),
        file.path(dir_out, "3.2_resultados.rds"))
cat("\nResultados guardados en:", dir_out, "\n")