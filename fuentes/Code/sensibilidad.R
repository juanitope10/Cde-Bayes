# ==============================================================================
# ANALISIS DE SENSIBILIDAD DE LAS PREVIAS + DIC/WAIC ENTRE VARIANTES (3.3)
# ------------------------------------------------------------------------------
# Corre el escenario sigma0 = 0.10 (misma semilla de ruido que la corrida base)
# con cinco variantes de previa que conservan las medias previas del modelo base
# cuando estas existen:
#   base        s2 ~ IG(2, 0.01)     tau ~ G(1, 0.05)    lam ~ G(2, 0.10)
#   debil       s2 ~ IG(1, 0.01)     tau ~ G(0.5, 0.025) lam ~ G(0.5, 0.025)
#   informativa s2 ~ IG(5, 0.04)     tau ~ G(5, 0.25)    lam ~ G(5, 0.25)
#   jeffreys    s2 ~ 1/s2, tau ~ 1/tau, lam ~ 1/lam  (impropias de escala)
#   unit_info   s2 ~ GI(1/2, var(y)/2) (informacion unitaria; tau y lam base)
# Cada corrida acumula la verosimilitud condicional punto a punto (WAIC) y se
# compara con DIC (verosimilitud marginal, via espectral) y WAIC entre variantes.
# Uso: CDE_LOGDET=espectral Rscript sensibilidad.R
# ==============================================================================
Sys.setenv(CDE_SOLO_FUNCIONES = "1")
suppressPackageStartupMessages({
  library(Matrix); library(png); library(coda); library(mvtnorm)
})
ruta <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
if (is.na(ruta)) ruta <- "sensibilidad.R"
source(file.path(dirname(ruta), "verificaciones.R"))   # ll_marg_33 y afines
source(file.path(dirname(ruta), "3.3.R"))              # funciones del 3.3 + x0

dir_sens <- file.path(dirname(dir_out), "sensibilidad")
dir.create(dir_sens, showWarnings = FALSE, recursive = TRUE)
dir_out <- dir_sens
cat("\nResultados de sensibilidad en:", dir_sens, "\n")

reja128 <- construir_reja(M, N)

variantes <- list(
  base        = list(a_sigma = 2,   b_sigma = 0.01,  a_tau = 1,   b_tau = 0.05,
                     a_lambda = 2,   b_lambda = 0.10),
  debil       = list(a_sigma = 1,   b_sigma = 0.01,  a_tau = 0.5, b_tau = 0.025,
                     a_lambda = 0.5, b_lambda = 0.025),
  informativa = list(a_sigma = 5,   b_sigma = 0.04,  a_tau = 5,   b_tau = 0.25,
                     a_lambda = 5,   b_lambda = 0.25),
  jeffreys    = list(a_sigma = 0,   b_sigma = 0,     a_tau = 0,   b_tau = 0,
                     a_lambda = 0,   b_lambda = 0),
  unit_info   = list(a_sigma = 0.5, b_sigma = NA,    a_tau = 1,   b_tau = 0.05,
                     a_lambda = 2,   b_lambda = 0.10)   # b_sigma = var(y)/2 al correr
)

resumen <- list(); cadenas <- list(); extras <- list()
for (v in names(variantes)) {
  cat(sprintf("\n################ variante: %s ################\n", v))
  cfg <- config
  cfg[names(variantes[[v]])] <- variantes[[v]]
  cfg$sigma0_levels <- 0.10
  cfg$waic <- TRUE
  sem <- semilla + 100                     # misma semilla del escenario s0 = 0.10
  set.seed(sem)
  yv <- x0 + rnorm(P, 0, 0.10)             # ruido identico al de la corrida base
  if (v == "unit_info") cfg$b_sigma <- 0.5 * var(yv)
  r <- correr_escenario(x0, 0.10, cfg, sem)
  post <- r$res$post

  # DIC (verosimilitud marginal espectral) y WAIC (condicional punto a punto)
  Yt <- transformar(r$y, M, N)
  ll <- vapply(post, function(t) {
    ll_marg_33(r$res$mu[t], r$res$sig2[t], r$res$tau[t], r$res$lambda[t],
               Yt, reja128$w, P)
  }, numeric(1))
  ll_bar <- ll_marg_33(mean(r$res$mu[post]), mean(r$res$sig2[post]),
                       mean(r$res$tau[post]), mean(r$res$lambda[post]),
                       Yt, reja128$w, P)
  p_D <- 2 * (ll_bar - mean(ll))
  dic <- -2 * ll_bar + 2 * p_D

  wa <- r$res$waic
  Tn <- wa$n
  lppd <- sum(log(wa$s1 / Tn))
  var_lp <- (wa$s3 - wa$s2^2 / Tn) / (Tn - 1)
  p_waic <- sum(var_lp)
  waic <- -2 * lppd + 2 * p_waic

  resumen[[v]] <- data.frame(
    variante = v,
    a_sigma = cfg$a_sigma, b_sigma = cfg$b_sigma,
    a_tau = cfg$a_tau, b_tau = cfg$b_tau,
    a_lambda = cfg$a_lambda, b_lambda = cfg$b_lambda,
    mu = mean(r$res$mu[post]), tau = mean(r$res$tau[post]),
    lambda = mean(r$res$lambda[post]), sig2 = mean(r$res$sig2[post]),
    mse = mean((r$res$x_media - x0)^2),
    acept_mh = mean(r$res$acc_par[post]),
    ll_media = mean(ll), p_D = p_D, DIC = dic,
    lppd = lppd, p_WAIC = p_waic, WAIC = waic,
    stringsAsFactors = FALSE
  )
  cadenas[[v]] <- list(lambda = r$res$lambda[post], sig2 = r$res$sig2[post],
                       tau = r$res$tau[post], mu = r$res$mu[post])
  extras[[v]] <- list(waic = wa, ll = ll)
  cat(sprintf("  DIC = %.1f (p_D = %.1f) | WAIC = %.1f (p_WAIC = %.1f) | MSE = %.5f\n",
              dic, p_D, waic, p_waic, mean((r$res$x_media - x0)^2)))
}

resumen <- do.call(rbind, resumen)
write.csv(resumen, file.path(dir_sens, "sensibilidad_dic_waic.csv"), row.names = FALSE)
cat("\n================ SENSIBILIDAD DE LAS PREVIAS (sigma0 = 0.10) ================\n")
print(resumen[, c("variante", "mu", "tau", "lambda", "sig2", "mse",
                  "DIC", "p_D", "WAIC", "p_WAIC")], row.names = FALSE, digits = 5)

# Figura: densidades posteriores bajo las cinco variantes
png(file.path(dir_sens, "sensibilidad_previas.png"), width = 1300, height = 900, res = 120)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
cols <- c(base = "#1b9e77", debil = "#d95f02", informativa = "#7570b3",
          jeffreys = "#e7298a", unit_info = "#66a61e")
for (p in c("mu", "tau", "lambda", "sig2")) {
  ds <- lapply(names(cadenas), function(v) density(cadenas[[v]][[p]]))
  names(ds) <- names(cadenas)
  xr <- range(unlist(lapply(ds, function(d) d$x)))
  yr <- range(unlist(lapply(ds, function(d) d$y)))
  plot(NULL, xlim = xr, ylim = yr, xlab = p, ylab = "densidad",
       main = paste("Sensibilidad de", p))
  for (v in names(ds)) lines(ds[[v]], col = cols[v], lwd = 2)
  legend("topright", legend = names(ds), col = cols, lwd = 2, bty = "n", cex = 0.7)
}
par(op); dev.off()

saveRDS(list(resumen = resumen, cadenas = cadenas),
        file.path(dir_sens, "sensibilidad_cadenas.rds"))
cat("\nListo. Archivos en:", dir_sens, "\n")
