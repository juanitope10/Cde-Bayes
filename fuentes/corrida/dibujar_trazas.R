# Trazas y ACF de los parámetros del 3.3 a partir del .rds de la corrida B = 24,000
# (3.3.R muestra estas gráficas en pantalla pero no las guarda en PNG).
args <- commandArgs(trailingOnly = TRUE)
base <- if (length(args)) args[1] else "."
r <- readRDS(file.path(base, "3.3", "3.3_resultados.rds"))
for (nm in names(r$resultados)) {
  res <- r$resultados[[nm]]$res
  s0  <- r$resultados[[nm]]$sigma0
  post <- res$post
  png(file.path(base, sprintf("diag_params_%.2f.png", s0)), width = 1300, height = 1500, res = 120)
  op <- par(mfrow = c(4, 2), mar = c(4, 4, 2.5, 1))
  for (p in list(list(res$mu, expression(mu)), list(res$tau, expression(tau)),
                 list(res$lambda, expression(lambda)), list(res$sig2, expression(sigma^2)))) {
    plot(p[[1]][post], type = "l", col = "grey30", xlab = "Iteración (tras calentamiento)",
         ylab = "", main = p[[2]])
    acf(p[[1]][post], lag.max = 60, main = "")
    title(main = "ACF")
  }
  par(op); dev.off()
}
