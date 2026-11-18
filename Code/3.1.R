# ==============================================================================
# CASO DE ESTUDIO - Preparación de imágenes y grafo de vecindad (4 vecinos)
# ==============================================================================
library(png)
library(Matrix)
library(ggplot2)

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
cat("Imágenes en:", dir_img, "\n")

# ==============================================================================
# 1. FUNCIONES DE CONVERSIÓN
# ==============================================================================
# as.vector() recorre por COLUMNAS; la Laplaciana usa el mismo orden (k = (c-1)*M + r)
matriz_a_vector <- function(mat) as.vector(mat)
vector_a_matriz <- function(vec, M = 128, N = 128) matrix(vec, nrow = M, ncol = N)

# {0 = negro, 1 = blanco} -> {+1, -1}. Umbral 0.5 por si el PNG trae valores intermedios
binario_a_ising <- function(img_bin) ifelse(img_bin < 0.5, 1, -1)
ising_a_binario <- function(ising_vec) ifelse(ising_vec == 1, 0, 1)

# Arreglo M x N x (3 o 4) -> matriz P x 3 (ignora canal alfa si existe)
rgb_a_matriz_canales <- function(arr_rgb) {
  stopifnot(length(dim(arr_rgb)) == 3, dim(arr_rgb)[3] >= 3)
  cbind(R = as.vector(arr_rgb[, , 1]),
        G = as.vector(arr_rgb[, , 2]),
        B = as.vector(arr_rgb[, , 3]))
}

matriz_canales_a_rgb <- function(mat_canales, M = 128, N = 128) {
  arr <- array(0, dim = c(M, N, 3))
  for (k in 1:3) arr[, , k] <- matrix(mat_canales[, k], nrow = M, ncol = N)
  arr
}

# ==============================================================================
# 2. LAPLACIANA DISPERSA (vectorizada, sin bucles)
# ==============================================================================
construir_laplaciana <- function(M = 128, N = 128) {
  P <- M * N
  idx <- matrix(1:P, nrow = M, ncol = N)
  
  # Aristas horizontales (r,c)-(r,c+1) y verticales (r,c)-(r+1,c)
  i <- c(as.vector(idx[, -N]), as.vector(idx[-M, ]))
  j <- c(as.vector(idx[, -1]), as.vector(idx[-1, ]))
  
  A <- sparseMatrix(i = c(i, j), j = c(j, i), x = 1, dims = c(P, P))
  d <- rowSums(A)
  L <- Diagonal(x = d) - A
  list(L = L, A = A, d = d, grid_idx = idx, n_aristas = length(i))
}

est_grafo <- construir_laplaciana(M, N)
L <- est_grafo$L; A <- est_grafo$A; d <- est_grafo$d

# ==============================================================================
# 2.2 VERIFICACIÓN
# ==============================================================================
verificar_grafo <- function(est, M = 128, N = 128) {
  d <- est$d
  cat("--- Verificación del grafo (P =", M * N, ") ---\n")
  cat("Esquinas (grado 2): obtenido", sum(d == 2), "| esperado 4\n")
  cat("Bordes   (grado 3): obtenido", sum(d == 3), "| esperado", 2 * (M - 2) + 2 * (N - 2), "\n")
  cat("Interior (grado 4): obtenido", sum(d == 4), "| esperado", (M - 2) * (N - 2), "\n")
  cat("Aristas: obtenido", est$n_aristas, "| esperado", M * (N - 1) + N * (M - 1), "\n")
  cat("L simétrica:", isSymmetric(est$L), "\n")
  cat("Filas de L suman 0:", max(abs(rowSums(est$L))) == 0, "\n")
  cat("Traza de L = suma de grados:", sum(diag(est$L)) == sum(d), "\n")
}
verificar_grafo(est_grafo, M, N)

# ==============================================================================
# 3. CARGA DE IMÁGENES (con validaciones)
# ==============================================================================
cargar_png <- function(archivo) {
  ruta <- file.path(dir_img, archivo)
  if (!file.exists(ruta)) stop("No existe: ", ruta)
  img <- readPNG(ruta)
  cat(sprintf("%-20s dim = %s | rango = [%.3f, %.3f]\n", archivo,
              paste(dim(img), collapse = " x "), min(img), max(img)))
  stopifnot(dim(img)[1] == M, dim(img)[2] == N)
  img
}

cat("\n--- Carga de imágenes ---\n")
img_horse  <- cargar_png("horse_binaria.png")
img_camera <- cargar_png("camera_grises.png")
img_astro  <- cargar_png("astronaut_rgb.png")
# Binaria -> Ising {-1, +1}
if (length(dim(img_horse)) == 3) img_horse <- img_horse[, , 1]
vec_horse_ising <- binario_a_ising(matriz_a_vector(img_horse))

# Grises en [0, 1]
if (length(dim(img_camera)) == 3) img_camera <- img_camera[, , 1]
vec_camera_grises <- matriz_a_vector(img_camera)

# RGB en [0, 1] (P x 3)
mat_astro_rgb <- rgb_a_matriz_canales(img_astro)

# Prueba de ida y vuelta: no se debe perder información
stopifnot(identical(ising_a_binario(vec_horse_ising),
                    as.vector((img_horse >= 0.5) * 1)))
stopifnot(all(vector_a_matriz(vec_camera_grises, M, N) == img_camera))
stopifnot(all(matriz_canales_a_rgb(mat_astro_rgb, M, N) == img_astro[, , 1:3]))
cat("\nConversiones ida y vuelta: OK\n")

# ==============================================================================
# 4. EXPLORACIÓN ESTADÍSTICA
# ==============================================================================
cat("\n--- Frecuencias imagen binaria (Ising) ---\n")
print(table(vec_horse_ising))
cat("Proporciones:\n"); print(round(prop.table(table(vec_horse_ising)), 4))

cat("\n--- Intensidades en escala de grises ---\n")
print(summary(vec_camera_grises))
cat("Desv. estándar:", round(sd(vec_camera_grises), 4), "\n")

cat("\n--- Canales RGB ---\n")
print(round(rbind(apply(mat_astro_rgb, 2, summary),
                  SD = apply(mat_astro_rgb, 2, sd)), 4))
cat("\nCorrelaciones entre canales:\n")
print(round(cor(mat_astro_rgb), 3))

# Imágenes originales
op <- par(mfrow = c(1, 3), mar = c(1, 1, 2, 1))
plot(as.raster(ising_a_binario(vec_horse_ising) |> vector_a_matriz(M, N)), main = "Binaria")
plot(as.raster(img_camera), main = "Escala de grises")
plot(as.raster(img_astro[, , 1:3]), main = "RGB")
par(op)

# Histograma de grises
print(
  ggplot(data.frame(x = vec_camera_grises), aes(x)) +
    geom_histogram(bins = 30, fill = "gray40", color = "white") +
    labs(title = "Distribución de intensidades (escala de grises)",
         x = expression(x[i]), y = "Frecuencia") +
    theme_minimal()
)

# Histogramas RGB
df_rgb <- data.frame(valor = as.vector(mat_astro_rgb),
                     canal = rep(c("R", "G", "B"), each = P))
df_rgb$canal <- factor(df_rgb$canal, levels = c("R", "G", "B"))
print(
  ggplot(df_rgb, aes(valor, fill = canal)) +
    geom_histogram(bins = 30, color = "white") +
    facet_wrap(~canal) +
    scale_fill_manual(values = c(R = "red", G = "forestgreen", B = "blue")) +
    labs(title = "Distribución por canal RGB", x = "Intensidad", y = "Frecuencia") +
    theme_minimal() + theme(legend.position = "none")
)

# (Bloque duplicado de figuras eliminado: repetía los mismos gráficos y usaba
#  dev.new(), que falla en sesiones no interactivas.)


