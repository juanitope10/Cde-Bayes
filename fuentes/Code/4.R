# ==============================================================================
# SECCIONES 4 y 5 - TABLA COMPARATIVA GLOBAL Y DISCUSIÓN
# ------------------------------------------------------------------------------
# Consolida los resultados de 3.2.R (Ising), 3.3.R (GMRF grises) y 3.4.R
# (GMRF RGB) y produce:
#   (i)   tabla comparativa global en consola,
#   (ii)  tabla comparativa en LaTeX (tabla_global.tex),
#   (iii) análisis breve de hallazgos (consola y discusion.txt).
# ==============================================================================
suppressPackageStartupMessages(library(utils))

# --- Rutas -------------------------------------------------------------------
# Las rutas del proyecto las define el archivo de configuración del sistema
# operativo: config_linux.R (Linux) o config_windows.R (Windows), ambos en Code/.
# Allí se definen dir_raiz, dir_code, dir_img y dir_out; aquí solo se usan.
if (!exists("dir_img")) {
  archivo_config <- if (identical(.Platform$OS.type, "windows")) "config_windows.R" else "config_linux.R"
  dirs <- c(".", "Code", "..", file.path("..", "Code"))
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
dir.create(file.path(dir_out, "4"), showWarnings = FALSE, recursive = TRUE)

archivos <- c("3.2/3.2_resultados.rds", "3.3/3.3_resultados.rds", "3.4/3.4_resultados.rds")
faltan <- archivos[!file.exists(file.path(dir_out, archivos))]
if (length(faltan)) {
  stop("Faltan resultados. Ejecute primero: ",
       paste(faltan, collapse = ", "), " (carpeta ", dir_out, ")")
}

res32 <- readRDS(file.path(dir_out, "3.2", "3.2_resultados.rds"))
res33 <- readRDS(file.path(dir_out, "3.3", "3.3_resultados.rds"))
res34 <- readRDS(file.path(dir_out, "3.4", "3.4_resultados.rds"))

# ------------------------------------------------------------------------------
# 1. TABLA COMPARATIVA GLOBAL
# ------------------------------------------------------------------------------
m32 <- res32$resumen
m33 <- res33$resumen
m34 <- res34$resumen

# CVMC medio sobre los parámetros ESCALARES (criterio del curso: CVMC < 5%),
# leído de los resúmenes guardados por cada script; el CVMC de los 16,384
# píxeles se reporta aparte (CVMC medio esc.+pix.). Se excluyen las entradas
# con media ~ 0 (p. ej. covarianzas de ruido), donde el CVMC no es interpretable.
cvmc_escalar <- function(tab) {
  keep <- abs(tab[, "Media"]) >= 0.01
  mean(tab[keep, "CVMC(%)"])
}
cvmc32 <- sapply(res32$resultados, function(r) cvmc_escalar(r$tab))
cvmc33 <- sapply(res33$resultados, function(r) cvmc_escalar(r$tab_res))
cvmc34 <- sapply(res34$resultados, function(r) cvmc_escalar(r$tab_res))

# Columnas unificadas: modelo, nivel de ruido, tiempo, n_eff mín, CVMC medio,
# MSE, SNR, PSNR y precisión MAP (solo imagen binaria)
armar_fila <- function(m, modelo, nivel_nombre, cvmc_esc, precision = NA_real_) {
  data.frame(
    Modelo = modelo,
    Ruido = sprintf("%s = %.2f", nivel_nombre, m$nivel_ruido),
    `Tiempo (min)` = round(m$tiempo_s / 60, 1),
    `n_eff min` = round(m$neff_min),
    `CVMC medio (%)` = round(cvmc_esc, 3),
    MSE = round(m$mse, 5),
    `SNR (dB)` = round(m$snr, 2),
    `PSNR (dB)` = round(m$psnr, 2),
    `Precisión MAP (%)` = round(100 * precision, 2),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

tab <- rbind(
  do.call(rbind, lapply(seq_len(nrow(m32)), function(i)
    armar_fila(m32[i, ], "Ising (binaria)", "rho0", cvmc32[i], m32$precision_MAP[i]))),
  do.call(rbind, lapply(seq_len(nrow(m33)), function(i)
    armar_fila(m33[i, ], "GMRF (grises)", "sigma0", cvmc33[i]))),
  do.call(rbind, lapply(seq_len(nrow(m34)), function(i)
    armar_fila(m34[i, ], "GMRF multivariado (RGB)", "sigma0", cvmc34[i])))
)

cat("\n================ TABLA COMPARATIVA GLOBAL ================\n")
cat("(n_eff mín sobre escalares y píxeles; CVMC medio sobre parámetros escalares)\n")
print(tab, row.names = FALSE)

cat("\n--- CVMC medio incluyendo los 16,384 píxeles (referencia) ---\n")
ref <- rbind(
  data.frame(Modelo = "Ising (binaria)", Ruido = sprintf("rho0 = %.2f", m32$nivel_ruido),
             `CVMC esc.+pix. (%)` = round(m32$cvmc_medio, 3), check.names = FALSE),
  data.frame(Modelo = "GMRF (grises)", Ruido = sprintf("sigma0 = %.2f", m33$nivel_ruido),
             `CVMC esc.+pix. (%)` = round(m33$cvmc_medio, 3), check.names = FALSE),
  data.frame(Modelo = "GMRF multivariado (RGB)", Ruido = sprintf("sigma0 = %.2f", m34$nivel_ruido),
             `CVMC esc.+pix. (%)` = round(m34$cvmc_medio, 3), check.names = FALSE)
)
print(ref, row.names = FALSE)

cat("\n--- RGB: MSE y SNR por canal ---\n")
print(data.frame(
  Ruido = sprintf("sigma0 = %.2f", m34$nivel_ruido),
  MSE_R = round(m34$mse_R, 5), MSE_G = round(m34$mse_G, 5), MSE_B = round(m34$mse_B, 5),
  SNR_R = round(m34$snr_R, 2), SNR_G = round(m34$snr_G, 2), SNR_B = round(m34$snr_B, 2)
), row.names = FALSE)

# ------------------------------------------------------------------------------
# 2. TABLA EN LATEX
# ------------------------------------------------------------------------------
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
latex_filas <- apply(tab, 1, function(r) {
  paste0(r[["Modelo"]], " & ", r[["Ruido"]], " & ",
         fmt(as.numeric(r[["Tiempo (min)"]]), 1), " & ",
         fmt(as.numeric(r[["n_eff min"]]), 0), " & ",
         fmt(as.numeric(r[["CVMC medio (%)"]]), 2), " & ",
         fmt(as.numeric(r[["MSE"]]), 5), " & ",
         fmt(as.numeric(r[["SNR (dB)"]]), 2), " & ",
         fmt(as.numeric(r[["PSNR (dB)"]]), 2), " & ",
         if (is.na(as.numeric(r[["Precisión MAP (%)"]]))) "--" else
           fmt(as.numeric(r[["Precisión MAP (%)"]]), 2), " \\\\")
})

latex <- c(
  "% Tabla comparativa global - Restauración bayesiana de imágenes",
  "% Generada por Code/4.R (no editar a mano)",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Desempeño de los tres modelos de restauración bayesiana bajo distintos niveles de ruido. $n_{\\mathrm{eff}}$ mín. calculado sobre parámetros escalares y píxeles; CVMC medio sobre los parámetros escalares (criterio CVMC $< 5\\%$).}",
  "\\label{tab:comparativa}",
  "\\small",
  "\\begin{tabular}{llrrrrrrr}",
  "\\hline",
  "Modelo & Ruido & Tiempo (min) & $n_{\\mathrm{eff}}$ mín. & CVMC (\\%) & MSE & SNR (dB) & PSNR (dB) & Precisión (\\%) \\\\",
  "\\hline",
  latex_filas,
  "\\hline",
  "\\end{tabular}",
  "\\end{table}",
  "",
  "% MSE y SNR por canal RGB",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{MSE y SNR por canal del modelo GMRF multivariado.}",
  "\\label{tab:rgb_canales}",
  "\\small",
  "\\begin{tabular}{lrrrrrr}",
  "\\hline",
  "Ruido & MSE R & MSE G & MSE B & SNR R & SNR G & SNR B \\\\",
  "\\hline",
  sapply(seq_len(nrow(m34)), function(i)
    sprintf("sigma0 = %.2f & %s & %s & %s & %s & %s & %s \\\\", m34$nivel_ruido[i],
            fmt(m34$mse_R[i], 5), fmt(m34$mse_G[i], 5), fmt(m34$mse_B[i], 5),
            fmt(m34$snr_R[i], 2), fmt(m34$snr_G[i], 2), fmt(m34$snr_B[i], 2))),
  "\\hline",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(latex, file.path(dir_out, "4", "tabla_global.tex"))
cat("\nTabla LaTeX escrita en:", file.path(dir_out, "4", "tabla_global.tex"), "\n")

# ------------------------------------------------------------------------------
# 3. ANÁLISIS BREVE DE HALLAZGOS
# ------------------------------------------------------------------------------
mejora <- function(mse_cont, mse_rest) round(100 * (1 - mse_rest / mse_cont), 1)

txt <- c(
  "ANÁLISIS BREVE DE HALLAZGOS",
  "===========================",
  "",
  "1) Papel del suavizado espacial.",
  sprintf(paste0(
    "En el modelo de Ising el parámetro beta controla la probabilidad de que dos ",
    "píxeles vecinos compartan estado: con beta posterior tipicamente > 0.4 ",
    "(por encima del valor crítico log(1+sqrt(2))/2 = 0.440 de la retícula cuadrada) ",
    "la previa favorece regiones homogéneas y la regla MAP limpia el ruido ",
    "aislado preservando los bordes de la silueta (ver precisión por tipo de ",
    "píxel en la salida de 3.2.R). En los GMRF, lambda pondera las diferencias ",
    "cuadráticas entre vecinos: lambda grande produce restauraciones mas ",
    "suaves, mientras que tau (grises) y Omega (RGB) controlan la variabilidad ",
    "marginal. La imagen restaurada de escala de grises alcanza un MSE ",
    "%.5f (SNR %.2f dB) frente a %.5f del ruido puro, una mejora del %.1f%%, ",
    "y el RGB global alcanza MSE %.5f (SNR %.2f dB), mejora del %.1f%%."),
    m33$mse[1], m33$snr[1], m33$mse_contaminada[1],
    mejora(m33$mse_contaminada[1], m33$mse[1]),
    m34$mse[1], m34$snr[1], mejora(m34$mse_contaminada[1], m34$mse[1])),
  "",
  "2) Efecto de la acumulación de ruido en la incertidumbre posterior.",
  "Al pasar de sigma0 = 0.05 a 0.20 (rho0 = 0.05 a 0.20 en el binario) el MSE",
  "aumenta y la SNR cae de forma consistente en los tres modelos. La",
  "distribución posterior lo refleja: las desviaciones estándar posteriores de",
  "los píxeles crecen con el ruido. En el Ising la incertidumbre se concentra en",
  "el contorno de la silueta y los detalles finos; en los GMRF la DE posterior es",
  "casi uniforme y mayor en el marco de la imagen (píxeles con 2 o 3 vecinos),",
  "porque no depende del contenido. La comparación de sigma^2 posterior contra el",
  "valor simulado en escala de grises revela subestimación sistemática:",
  sprintf("medias posteriores %.5f, %.5f y %.5f frente a valores reales %.5f, %.5f y %.5f (sigma0 = 0.05, 0.10, 0.20).",
          m33$sigma2_media[m33$nivel_ruido == 0.05], m33$sigma2_media[m33$nivel_ruido == 0.10],
          m33$sigma2_media[m33$nivel_ruido == 0.20],
          0.05^2, 0.10^2, 0.20^2),
  "Este sesgo es más fuerte a ruido BAJO y se debe a la confusión entre la",
  "varianza del ruido y la rugosidad del campo latente: a sigma0 = 0.05 la",
  "posterior atribuye el ruido fino a la textura de X (sigma^2 colapsa y la",
  "restauración apenas mejora el dato contaminado), mientras que a ruido",
  "moderado y alto el modelo separa mejor ambos componentes y la mejora de MSE",
  "supera el 50%. Es una propiedad del modelo GMRF con varianza de ruido",
  "desconocida, no un defecto del muestreador (verificado con cadenas largas y",
  "distintas inicializaciones). En RGB, SIGMA muestra el mismo comportamiento,",
  "concentrado en la dirección de brillo (1,1,1)/sqrt(3): el ruido común a los",
  "tres canales se confunde con la imagen, mientras que el ruido de color se",
  "estima bien; de ahí las correlaciones residuales negativas entre canales.",
  "",
  "3) Ventajas y limitaciones computacionales del enfoque MCMC.",
  sprintf(paste0(
    "Los tres algoritmos procesan P = 16,384 píxeles con algebra dispersa",
    " (Matrix) y factorizaciones de Cholesky: los tiempos medianos por nivel",
    " de ruido fueron ~%.0f min (Ising), ~%.0f min (grises) y ~%.0f min (RGB).",
    " Las ventajas: inferencia completa (medias, varianzas e intervalos por",
    " píxel), sin supuestos de normalidad para la imagen binaria y con",
    " propagación natural de la incertidumbre. Las limitaciones: (i) costo que",
    " crece rapidamente con la resolucion (cada barrido toca todos los píxeles",
    " y cada paso MH exige factorizar matrices de precision); (ii) el algoritmo",
    " de intercambio para beta depende de barridos auxiliares K (aproximación),",
    " y la mezcla de beta, tau y lambda exige cadenas largas para alcanzar",
    " CVMC < 5%%; (iii) los modelos GMRF suavizan bordes agudos (lambda global),",
    " mientras que Ising, al ser discreto, los preserva mejor."),
    median(m32$tiempo_s) / 60, median(m33$tiempo_s) / 60, median(m34$tiempo_s) / 60),
  "",
  "4) Comparación entre modelos.",
  sprintf(paste0(
    "El modelo binario recupera el %.1f%% de los píxeles con rho0 = 0.10",
    " (vs %.1f%% del dato contaminado), aunque subestima rho (media posterior",
    " %.3f frente a 0.10, con el IC95%% sin contener el valor simulado): parte de",
    " los píxeles volteados se confunde con estructura real de la imagen. En",
    " escala de grises la mejora MSE es del %.1f%% y en RGB del %.1f%%, con",
    " prestaciones similares entre canales (ver tabla por canal). En conjunto,",
    " incorporar dependencia espacial mejora la restauración en todos los",
    " niveles de ruido, con ganancias relativas que en los GMRF crecen con el ruido y una",
    " incertidumbre posterior que crece de forma coherente con la contaminación."),
    100 * m32$precision_MAP[1], 100 * m32$precision_contaminada[1], m32$rho_media[1],
    mejora(m33$mse_contaminada[1], m33$mse[1]),
    mejora(m34$mse_contaminada[1], m34$mse[1]))
)

cat("\n")
writeLines(txt, file.path(dir_out, "4", "discusion.txt"))
for (l in txt) cat(l, "\n", sep = "")
cat("\nDiscusión escrita en:", file.path(dir_out, "4", "discusion.txt"), "\n")
