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

dir_img <- "C:/Users/Estudiante/Downloads/caso/caso/imagenes"
if (!dir.exists(dir_img)) dir_img <- "imagenes"          # alternativa relativa
dir_out <- "resultados_ising"
dir.create(dir_out, showWarnings = FALSE)

config <- list(
  # Hiperparámetros (previas débilmente informativas)
  a_beta = 2, b_beta = 4,      # beta ~ Gamma(2, 4): media 0.5, sd 0.35
  a_rho  = 1, b_rho  = 1,      # rho  ~ Beta(1, 1) truncada en (0, 0.5) = Unif(0, 0.5)
  # Algoritmo
  B      = 4000,               # iteraciones totales
  burn   = 1000,               # calentamiento
  delta  = 0.02,               # desv. estándar de la propuesta de beta (aceptación ~26% con K=25)
  K      = 25,                 # barridos de Gibbs para X*. OJO: la media de beta baja al subir K (sesgo); verifica con K=50
  beta0  = 0.5,                # valor inicial de beta
  rho_ini = 0.10               # valor inicial de rho
)
niveles_ruido <- c(0.10, 0.05, 0.20)   # el principal primero (rho_0 = 0.10)
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
  x <- x_ini; beta <- beta0; rho <- rho_ini
  beta_ch <- rho_ch <- S_ch <- numeric(B)
  acepta  <- logical(B)
  conteo_mas <- numeric(length(y))             # para pi_i = P(X_i = +1 | y)
  n_post <- 0
  t0 <- Sys.time()
  
  for (b in seq_len(B)) {
    x    <- barrido_gibbs(x, beta, tab, y = y, rho = rho)        # 1) X | beta, rho, y
    rho  <- muestrear_rho(x, y, cfg$a_rho, cfg$b_rho)            # 2) rho | x, y
    S_x  <- S_concordancia(x, A)
    pb   <- paso_beta(beta, x, S_x, A, tab, cfg)                 # 3) beta | x (intercambio)
    beta <- pb$beta
    beta_ch[b] <- beta; rho_ch[b] <- rho; S_ch[b] <- S_x; acepta[b] <- pb$acepta
    if (b > cfg$burn) { conteo_mas <- conteo_mas + (x == 1); n_post <- n_post + 1 }
    
    if (verbose && b %% (B / 10) == 0)
      cat(sprintf("  [%3d%%] iter %5d | beta = %.4f | rho = %.4f | S/aristas = %.4f | acept(beta) = %.1f%% | %.0fs\n",
                  round(100 * b / B), b, beta, rho, S_x / (M * (N - 1) + N * (M - 1)),
                  100 * mean(acepta[1:b]), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  list(beta = beta_ch, rho = rho_ch, S = S_ch, acepta = acepta,
       pi = conteo_mas / n_post, x_final = x)
}

# ------------------------------------------------------------------------------
# 5. DIAGNÓSTICOS Y RESÚMENES
# ------------------------------------------------------------------------------
resumen_posterior <- function(cadena) {
  ch    <- mcmc(cadena)
  media <- mean(cadena); dsd <- sd(cadena)
  neff  <- as.numeric(effectiveSize(ch))
  emc   <- dsd / sqrt(neff)
  c(Media = media, DE = dsd, EMC = emc, n_eff = neff, `CVMC(%)` = 100 * emc / abs(media),
    Geweke_z = as.numeric(geweke.diag(ch)$z))
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
  cat(sprintf("Tasa de aceptación de beta (post burn-in): %.1f%%\n", 100 * tasa))
  
  x_map <- ifelse(res$pi >= 0.5, 1, -1)
  prec  <- precision_restauracion(x_map, x0, y, A)
  cat(sprintf("Píxeles correctos: imagen contaminada = %.2f%% | restaurada (MAP) = %.2f%%\n",
              100 * prec["Y_global"], 100 * prec["MAP_global"]))
  cat(sprintf("  en regiones homogéneas = %.2f%% | en bordes = %.2f%% (bordes = %.1f%% de los píxeles)\n",
              100 * prec["MAP_homogeneas"], 100 * prec["MAP_bordes"], 100 * prec["prop_bordes"]))
  cat(sprintf("rho_0 real = %.3f | IC95%% posterior de rho = [%.4f, %.4f]\n", rho0,
              quantile(rho_p, .025), quantile(rho_p, .975)))
  
  c(list(rho0 = rho0, x0 = x0, y = y, x_map = x_map, tab = tab_res, tasa = tasa,
         prec = prec, post = post), res)
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
saveRDS(resultados, file.path(dir_out, "resultados_ising.rds"))