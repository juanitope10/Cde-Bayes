# ==============================================================================
# TESTS DE CONVERGENCIA DE LAS TRAZAS (3.2, 3.3 y 3.4)
# ------------------------------------------------------------------------------
# Calcula sobre las cadenas guardadas en los .rds de cada punto:
#   - R-hat de Gelman-Rubin sobre segmentos (split) y R-hat por rangos
#     (rank-normalized split, Vehtari et al. 2021), con su IC superior (coda)
#   - ESS bulk y ESS tail (rangos)
#   - Geweke z, Heidelberger-Welch (estacionariedad y semi-amplitud),
#     Raftery-Lewis (iteraciones necesarias)
#   - ACF en retardos 1, 5, 10, 20, 50 y CVMC
#   - Incertidumbre de la asignacion: pixeles ambiguos del MAP (3.2) y
#     pixeles con DE posterior mayor que el ruido (3.3, 3.4)
# Salidas en Code/resultados/convergencia/: CSV por modelo, resumen global,
# ambiguedad.csv y PNG de trazas por segmentos + ACF de parametros representativos.
# Uso: Rscript convergencia.R
# ==============================================================================
suppressPackageStartupMessages({
  library(coda)
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
dir_conv <- file.path(dir_out, "convergencia")
dir.create(dir_conv, showWarnings = FALSE, recursive = TRUE)
cat("Resultados de convergencia en:", dir_conv, "\n")

# ------------------------------------------------------------------------------
# 1. ESTADISTICOS
# ------------------------------------------------------------------------------
rhat_cadena <- function(mat) {
  if (is.null(dim(mat)) || length(dim(mat)) != 2) return(NA_real_)
  n <- nrow(mat); m <- ncol(mat)
  if (m < 2 || n < 10) return(NA_real_)
  W <- mean(apply(mat, 2, var))
  B <- n * var(colMeans(mat))
  if (!is.finite(W) || W <= 0) return(1)
  sqrt(((n - 1) / n * W + B / n) / W)
}

rank_norm <- function(x) {
  S <- length(x)
  qnorm((rank(x, ties.method = "average") - 3/8) / (S + 1 - 2 * 3/8))
}

partir_cadena <- function(x, n_seg = 4) {
  n <- length(x)
  tam <- floor(n / n_seg)
  matrix(x[seq_len(tam * n_seg)], nrow = tam, ncol = n_seg)
}

rhat_rank_split <- function(x, n_seg = 4) {
  mat <- partir_cadena(x, n_seg)
  z <- rank_norm(as.vector(mat)); dim(z) <- dim(mat)
  r_bulk <- rhat_cadena(z)
  med <- median(x)
  zf <- rank_norm(abs(as.vector(mat) - med)); dim(zf) <- dim(mat)
  r_fold <- rhat_cadena(zf)
  max(r_bulk, r_fold, na.rm = TRUE)
}

rhat_clasico_split <- function(x, n_seg = 4) {
  rhat_cadena(partir_cadena(x, n_seg))
}

rhat_sup_coda <- function(x, n_seg = 4) {
  mat <- partir_cadena(x, n_seg)
  ml <- mcmc.list(lapply(seq_len(ncol(mat)), function(j) mcmc(mat[, j])))
  tryCatch({
    g <- gelman.diag(ml, autoburnin = FALSE, multivariate = FALSE)$psrf
    unname(g[1, "Upper C.I."])
  }, error = function(e) NA_real_)
}

ess_bulk <- function(x) {
  z <- rank_norm(x)
  v <- suppressWarnings(as.numeric(effectiveSize(mcmc(z))))
  if (!is.finite(v) || v <= 0) length(x) else v
}

ess_tail <- function(x) {
  q <- quantile(x, probs = c(0.05, 0.95))
  e1 <- suppressWarnings(as.numeric(effectiveSize(mcmc(as.numeric(x <= q[1])))))
  e2 <- suppressWarnings(as.numeric(effectiveSize(mcmc(as.numeric(x <= q[2])))))
  e <- c(e1, e2)
  e[!is.finite(e) | e <= 0] <- length(x)
  min(e)
}

diag_param <- function(x, etiqueta) {
  n <- length(x)
  media <- mean(x); dsd <- sd(x)
  neff <- suppressWarnings(as.numeric(effectiveSize(mcmc(x))))
  if (!is.finite(neff) || neff <= 0) neff <- n
  emc <- dsd / sqrt(neff)
  cvmc <- 100 * emc / max(abs(media), 1e-12)
  gw <- tryCatch(geweke.diag(mcmc(x))$z[1], error = function(e) NA_real_)
  hw <- tryCatch({
    h <- heidel.diag(mcmc(x))
    c(p_est = h[1, "pvalue"], hw_ok = unname(h[1, "htest"]))
  }, error = function(e) c(p_est = NA_real_, hw_ok = NA_real_))
  rf <- tryCatch({
    r <- raftery.diag(mcmc(x), q = 0.025, r = 0.0125, s = 0.95)
    c(N_req = unname(r$resmatrix[1, "N"]), Nmin = unname(r$resmatrix[1, "Nmin"]))
  }, error = function(e) c(N_req = NA_real_, Nmin = NA_real_))
  ac <- acf(x, lag.max = 50, plot = FALSE)$acf[-1]
  data.frame(
    parametro = etiqueta, n_muestras = n, media = media, DE = dsd,
    n_eff = neff, CVMC_pct = cvmc,
    ESS_bulk = ess_bulk(x), ESS_tail = ess_tail(x),
    Rhat_split = rhat_clasico_split(x), Rhat_rank = rhat_rank_split(x),
    Rhat_sup = rhat_sup_coda(x),
    Geweke_z = gw, Geweke_ok = abs(gw) < 1.96,
    HW_p_estacionariedad = unname(hw["p_est"]),
    HW_media_ok = unname(hw["hw_ok"]) == 1,
    Raftery_N_req = unname(rf["N_req"]), Raftery_ok = unname(rf["N_req"]) <= n,
    ACF1 = ac[1], ACF5 = ac[5], ACF10 = ac[10], ACF20 = ac[20], ACF50 = ac[50],
    stringsAsFactors = FALSE
  )
}

# Trazas divididas en segmentos (color) + ACF, por parametro representativo
dibujar_split_acf <- function(x, etiqueta, color_seg) {
  n_seg <- length(color_seg)
  mat <- partir_cadena(x, n_seg)
  plot(NULL, xlim = c(1, nrow(mat)), ylim = range(mat),
       xlab = "Iteracion (tras calentamiento, por segmento)", ylab = etiqueta,
       main = paste("Traza por segmentos:", etiqueta))
  for (j in seq_len(n_seg)) lines(mat[, j], col = color_seg[j], lwd = 0.8)
  legend("topright", legend = paste("seg", seq_len(n_seg)), col = color_seg,
         lwd = 2, bty = "n", cex = 0.8)
  acf(x, lag.max = 60, main = paste("ACF:", etiqueta))
}

correr_modelo <- function(nombre, ruta_rds, extraer, representativos) {
  if (!file.exists(ruta_rds)) {
    cat("AVISO: no existe", ruta_rds, "- se omite", nombre, "\n")
    return(NULL)
  }
  obj <- readRDS(ruta_rds)
  tablas <- list()
  for (nm in names(obj$resultados)) {
    esc <- obj$resultados[[nm]]
    cadena <- extraer(esc)
    post <- if (!is.null(esc$res$post)) esc$res$post else esc$post
    if (is.null(post)) post <- seq_along(cadena[[1]])
    filas <- lapply(names(cadena), function(p) {
      diag_param(cadena[[p]][post], p)
    })
    tab <- do.call(rbind, filas)
    nivel <- if (!is.null(esc$sigma0)) esc$sigma0 else esc$rho0
    tab <- cbind(escenario = nm, nivel_ruido = nivel, tab, stringsAsFactors = FALSE)
    tablas[[nm]] <- tab
  }
  tab <- do.call(rbind, tablas)
  write.csv(tab, file.path(dir_conv, paste0("convergencia_", nombre, ".csv")),
            row.names = FALSE)
  cat("\n====", nombre, "====\n")
  print(tab[, c("escenario", "parametro", "n_eff", "Rhat_rank", "Rhat_split",
                "Rhat_sup", "Geweke_z", "HW_p_estacionariedad", "Raftery_ok")],
        row.names = FALSE, digits = 3)

  col_seg <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
  for (nm in names(obj$resultados)) {
    esc <- obj$resultados[[nm]]
    cadena <- extraer(esc)
    post <- if (!is.null(esc$res$post)) esc$res$post else esc$post
    if (is.null(post)) post <- seq_along(cadena[[1]])
    params <- intersect(representativos, names(cadena))
    nivel <- if (!is.null(esc$sigma0)) esc$sigma0 else esc$rho0
    archivo <- file.path(dir_conv, sprintf("convergencia_%s_%.2f.png", nombre, nivel))
    png(archivo, width = 1300, height = 260 * length(params) + 80, res = 120)
    op <- par(mfrow = c(length(params), 2), mar = c(4, 4, 3, 1))
    for (p in params) dibujar_split_acf(cadena[[p]][post], p, col_seg)
    par(op); dev.off()
  }
  tab
}

# ------------------------------------------------------------------------------
# 2. MODELOS
# ------------------------------------------------------------------------------
cat("\n----- 3.2 Ising -----\n")
t32 <- correr_modelo("3.2", file.path(dir_out, "3.2", "3.2_resultados.rds"),
  function(esc) list(beta = esc$beta, rho = esc$rho),
  c("beta", "rho"))

cat("\n----- 3.3 GMRF grises -----\n")
t33 <- correr_modelo("3.3", file.path(dir_out, "3.3", "3.3_resultados.rds"),
  function(esc) list(mu = esc$res$mu, tau = esc$res$tau,
                     lambda = esc$res$lambda, sig2 = esc$res$sig2),
  c("mu", "tau", "lambda", "sig2"))

ruta34 <- file.path(dir_out, "3.4", "3.4_resultados.rds")
if (file.exists(ruta34)) {
  n34 <- length(readRDS(ruta34)$resultados[[1]]$res$lambda)
  if (n34 < 1000) {
    cat("AVISO: 3.4_resultados.rds es la corrida rapida (", n34, "iter );\n",
        "        se usa 3.4_resultados_B6000_cholesky.rds para los tests.\n")
    ruta34 <- file.path(dir_out, "3.4", "3.4_resultados_B6000_cholesky.rds")
  }
} else {
  ruta34 <- file.path(dir_out, "3.4", "3.4_resultados_B6000_cholesky.rds")
}
cat("\n----- 3.4 GMRF RGB (", basename(ruta34), ") -----\n")
t34 <- correr_modelo("3.4", ruta34,
  function(esc) {
    r <- esc$res
    out <- list(mu_R = r$mu[, 1], mu_G = r$mu[, 2], mu_B = r$mu[, 3],
                lambda = r$lambda)
    et <- c("11", "21", "31", "22", "32", "33")
    for (k in 1:6) out[[paste0("Sigma_", et[k])]] <- r$Sigma[, k]
    for (k in 1:6) out[[paste0("Omega_", et[k])]] <- r$Omega[, k]
    out
  },
  c("lambda", "mu_R", "Sigma_11", "Omega_11"))

# ------------------------------------------------------------------------------
# 3. INCERTIDUMBRE DE LA ASIGNACION (pixeles ambiguos)
# ------------------------------------------------------------------------------
cat("\n----- Incertidumbre de la asignacion -----\n")
amb <- list()

obj32 <- readRDS(file.path(dir_out, "3.2", "3.2_resultados.rds"))
for (nm in names(obj32$resultados)) {
  esc <- obj32$resultados[[nm]]
  pi_i <- esc$pi
  ambiguo <- pi_i >= 0.4 & pi_i <= 0.6
  err_map <- esc$x_map != esc$x0
  amb[[length(amb) + 1]] <- data.frame(
    modelo = "Ising", escenario = nm, nivel_ruido = esc$rho0,
    criterio = "pi_i en [0.4, 0.6]",
    pct_ambiguos = 100 * mean(ambiguo),
    pct_errores_en_ambiguos = if (any(err_map)) 100 * mean(ambiguo[err_map]) else 0,
    precision_ambiguos = 100 * mean(!err_map[ambiguo]),
    precision_no_ambiguos = 100 * mean(!err_map[!ambiguo]),
    stringsAsFactors = FALSE
  )
}

for (mod in list(list(nombre = "GMRF grises", ruta = file.path(dir_out, "3.3", "3.3_resultados.rds")),
                 list(nombre = "GMRF RGB", ruta = ruta34))) {
  if (!file.exists(mod$ruta)) next
  obj <- readRDS(mod$ruta)
  for (nm in names(obj$resultados)) {
    esc <- obj$resultados[[nm]]
    dsd <- if (mod$nombre == "GMRF grises") esc$x_sd else rowMeans(esc$x_sd)
    amb[[length(amb) + 1]] <- data.frame(
      modelo = mod$nombre, escenario = nm, nivel_ruido = esc$sigma0,
      criterio = sprintf("DE posterior (homogeneidad p90/mediana = %.3f): mediana %.4f | max %.4f",
                         quantile(dsd, 0.9) / median(dsd), median(dsd), max(dsd)),
      pct_ambiguos = NA_real_,
      pct_errores_en_ambiguos = NA_real_, precision_ambiguos = NA_real_,
      precision_no_ambiguos = NA_real_, stringsAsFactors = FALSE
    )
  }
}

amb <- do.call(rbind, amb)
write.csv(amb, file.path(dir_conv, "ambiguedad.csv"), row.names = FALSE)
print(amb, row.names = FALSE, digits = 3)

# ------------------------------------------------------------------------------
# 4. RESUMEN GLOBAL
# ------------------------------------------------------------------------------
resumen <- do.call(rbind, list(t32, t33, t34))
write.csv(resumen, file.path(dir_conv, "resumen_convergencia.csv"), row.names = FALSE)
cat("\n================ RESUMEN DE CONVERGENCIA ================\n")
cat(sprintf("Parametros con Rhat_rank > 1.01: %d de %d\n",
            sum(resumen$Rhat_rank > 1.01, na.rm = TRUE), nrow(resumen)))
cat(sprintf("Parametros con Rhat_rank > 1.05: %d de %d\n",
            sum(resumen$Rhat_rank > 1.05, na.rm = TRUE), nrow(resumen)))
cat(sprintf("Parametros con Geweke |z| > 1.96: %d de %d\n",
            sum(!resumen$Geweke_ok, na.rm = TRUE), nrow(resumen)))
cat(sprintf("Parametros que no pasan Heidelberger-Welch (estacionariedad): %d de %d\n",
            sum(resumen$HW_p_estacionariedad < 0.05, na.rm = TRUE), nrow(resumen)))
cat(sprintf("Parametros con CVMC >= 5%%: %d de %d\n",
            sum(resumen$CVMC_pct >= 5, na.rm = TRUE), nrow(resumen)))
cat("\nListo. Archivos en:", dir_conv, "\n")
