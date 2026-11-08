# Cargar librerías necesarias
library(png)      # Carga de imágenes PNG
library(Matrix)   # Matrices dispersas eficientes
library(ggplot2)  # Graficación y exploración de distribuciones

# ==============================================================================
# 1.1 CONVERSIÓN MATRIZ / VECTOR (Grises y Binario)
# ==============================================================================

# Convierte matriz M x N a vector de longitud P = M*N (por filas o columnas)
matriz_a_vector <- function(mat) {
  return(as.vector(mat))
}

# Reconstruye matriz M x N a partir de vector de longitud P
vector_a_matriz <- function(vec, M = 128, N = 128) {
  return(matrix(vec, nrow = M, ncol = N))
}

# ==============================================================================
# 1.2 MAPEO DE VALORES BINARIOS (0/1 <-> -1/+1 para Modelo de Ising)
# ==============================================================================

# Convierte valores de {0, 1} (0=negro, 1=blanco) a {-1, +1} (+1=negro, -1=blanco)
binario_a_ising <- function(img_bin) {
  # Si 0 es negro (silueta) y 1 es blanco (fondo) [2]:
  # Mapeamos: 0 -> +1 (negro), 1 -> -1 (blanco)
  ising_vec <- ifelse(img_bin == 0, 1, -1)
  return(ising_vec)
}

# Proceso inverso: {-1, +1} -> {0, 1} para visualización estándar
ising_a_binario <- function(ising_vec) {
  img_bin <- ifelse(ising_vec == 1, 0, 1)
  return(img_bin)
}

# ==============================================================================
# 1.3 CONVERSIÓN RGB (M x N x 3 <-> P x 3)
# ==============================================================================

# Convierte arreglo 3D (M x N x 3) a matriz P x 3 (R, G, B por columnas)
rgb_a_matriz_canales <- function(arr_rgb) {
  M <- dim(arr_rgb)[3]
  N <- dim(arr_rgb)[2]
  P <- M * N
  
  mat_canales <- cbind(
    R = as.vector(arr_rgb[, , 1]),
    G = as.vector(arr_rgb[, , 2]),
    B = as.vector(arr_rgb[, , 3])
  )
  return(mat_canales) # Dimensión P x 3
}

# Reconstruye el arreglo tridimensional (M x N x 3) desde la matriz P x 3
matriz_canales_a_rgb <- function(mat_canales, M = 128, N = 128) {
  arr_rgb <- array(0, dim = c(M, N, 3))
  arr_rgb[, , 1] <- matrix(mat_canales[, "R"], nrow = M, ncol = N)
  arr_rgb[, , 2] <- matrix(mat_canales[, "G"], nrow = M, ncol = N)
  arr_rgb[, , 3] <- matrix(mat_canales[, "B"], nrow = M, ncol = N)
  return(arr_rgb)
}

#### Matriz laplaciana 

# ==============================================================================
# 2. CONSTRUCCIÓN DE LA MATRIZ LAPLACIANA DISPERSA (4 VECINOS)
# ==============================================================================

construir_laplaciana <- function(M = 128, N = 128) {
  P <- M * N
  
  # Matriz para convertir coordenadas 2D (r, c) a índice 1D k
  # k = (c - 1) * M + r (indexación por columnas estándar de R)
  grid_idx <- matrix(1:P, nrow = M, ncol = N)
  
  row_indices <- c()
  col_indices <- c()
  
  # Construcción de aristas horizontales y verticales (4-vecindad)
  for (r in 1:M) {
    for (c in 1:N) {
      k_actual <- grid_idx[r, c]
      
      # Vecino Derecho (r, c+1)
      if (c < N) {
        k_derecho <- grid_idx[r, c + 1]
        row_indices <- c(row_indices, k_actual, k_derecho)
        col_indices <- c(col_indices, k_derecho, k_actual)
      }
      # Vecino Abajo (r+1, c)
      if (r < M) {
        k_abajo <- grid_idx[r + 1, c]
        row_indices <- c(row_indices, k_actual, k_abajo)
        col_indices <- c(col_indices, k_abajo, k_actual)
      }
    }
  }
  
  # Matriz de Adyacencia A (dispersa)
  A <- sparseMatrix(i = row_indices, j = col_indices, x = 1, dims = c(P, P))
  
  # Grados d_i = suma por filas de A
  d <- rowSums(A)
  D <- diagonalMatrix <- Diagonal(x = d)
  
  # Matriz Laplaciana L = D - A
  L <- D - A
  
  return(list(L = L, A = A, d = d, grid_idx = grid_idx))
}

# Construir la estructura para 128x128
est_grafo <- construir_laplaciana(128, 128)
L <- est_grafo$L
d <- est_grafo$d

# ==============================================================================
# 2.2 VERIFICACIÓN DE VECINOS (ESQUINAS, BORDES E INTERIORES)
# ==============================================================================

verificar_vecinos <- function(d, M = 128, N = 128) {
  cat("--- Verificación de Grados de Vecindad (P =", M*N, ") ---\n")
  cat("Esquinas (esperado 4) :", sum(d == 2), "\n")
  cat("Bordes   (esperado", 2*(M-2) + 2*(N-2), "):", sum(d == 3), "\n")
  cat("Interior (esperado", (M-2)*(N-2), "):", sum(d == 4), "\n")
}

verificar_vecinos(d, 128, 128)


##### Carga, Conversión y Exploración de las Imágenes


# 3. CARGA Y PROCESAMIENTO DE LAS TRES IMÁGENES

# 3.1 Imagen Binaria (Horse)
img_horse_raw <- readPNG("C:/Users/Estudiante/Downloads/caso/caso/imagenes/horse_binaria.png")
if (length(dim(img_horse_raw)) == 3) img_horse_raw <- img_horse_raw[, , 1] # Tomar solo un canal si es 2D
vec_horse_bin  <- matriz_a_vector(img_horse_raw)
vec_horse_ising <- binario_a_ising(vec_horse_bin)  # Mapeo a {-1, +1}

# 3.2 Imagen Escala de Grises (Camera)
img_camera_raw <- readPNG("C:/Users/Estudiante/Downloads/caso/caso/imagenes/camera_grises.png")
if (length(dim(img_camera_raw)) == 3) img_camera_raw <- img_camera_raw[, , 1]
vec_camera_grises <- matriz_a_vector(img_camera_raw) # Intensidades en [0, 1]

# 3.3 Imagen RGB (Astronaut)
img_astro_raw <- readPNG("C:/Users/Estudiante/Downloads/caso/caso/imagenes/astronaut_rgb.png")
mat_astro_rgb <- rgb_a_matriz_canales(img_astro_raw)  # Matriz P x 3 con valores en [0, 1]

# ==============================================================================
# 4. EXPLORACIÓN ESTADÍSTICA DE CARACTERÍSTICAS
# ==============================================================================

# Frecuencias de la imagen binaria (Ising {-1, +1})
cat("\n--- Frecuencias de Imagen Binaria (Ising) ---\n")
print(table(vec_horse_ising))

# Resumen de intensidades en Escala de Grises
cat("\n--- Resumen de Intensidades en Escala de Grises ---\n")
print(summary(vec_camera_grises))

# Resumen por Canales RGB
cat("\n--- Resumen por Canales RGB ---\n")
print(apply(mat_astro_rgb, 2, summary))

# Histograma de Escala de Grises
hist(vec_camera_grises, breaks = 30, col = "gray40", main = "Distribución de Intensidades (Escala de Grises)",
     xlab = "Intensidad x_i", ylab = "Frecuencia")

# Histogramas por Canal RGB
par(mfrow = c(1, 3))
hist(mat_astro_rgb[, "R"], col = "red", main = "Canal Rojo (R)", xlab = "Intensidad")
hist(mat_astro_rgb[, "G"], col = "green", main = "Canal Verde (G)", xlab = "Intensidad")
hist(mat_astro_rgb[, "B"], col = "blue", main = "Canal Azul (B)", xlab = "Intensidad")
par(mfrow = c(1, 1))
