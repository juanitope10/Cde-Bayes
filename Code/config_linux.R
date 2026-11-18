# ==============================================================================
# CONFIGURACIÓN DE RUTAS - Linux
# ==============================================================================
# Este archivo define las variables de rutas del proyecto:
#   dir_raiz : carpeta del proyecto (la que contiene caso.Rmd e imagenes/)
#   dir_code : carpeta Code/ (donde están los scripts)
#   dir_img  : carpeta imagenes/ (imágenes de referencia del caso)
#   dir_out  : carpeta Code/resultados/ (salidas de los scripts)
#
# Los scripts (3.1.R, 3.2.R, 3.3.R, ...) cargan este archivo automáticamente y
# usan esas variables. En Windows se usa el equivalente config_windows.R.
# Suba ambos archivos a GitHub: cada usuario ejecutará el de su sistema.
#
# La detección es automática (este archivo vive en Code/). Solo si fallara,
# defina dir_raiz a mano descomentando la línea marcada [EDITAR].
# ==============================================================================

# [EDITAR solo si la detección automática falla] --------------------------------
# dir_raiz <- "/home/USUARIO/Desktop/Cde bayes"
# ------------------------------------------------------------------------------

if (!exists("dir_raiz")) {
  # 1) Si este archivo se cargó con source(), su ruta está en el frame de source()
  .ruta_cfg <- NULL
  for (i in rev(seq_len(sys.nframe()))) {
    g <- sys.frame(i)$ofile
    if (!is.null(g) && nzchar(g)) { .ruta_cfg <- g; break }
  }
  if (!is.null(.ruta_cfg)) {
    # .../Cde bayes/Code/config_linux.R -> .../Cde bayes
    dir_raiz <- dirname(dirname(normalizePath(.ruta_cfg)))
  } else {
    # 2) Pegado en la consola: subir desde getwd() hasta encontrar imagenes/
    q <- normalizePath(getwd())
    dir_raiz <- NULL
    for (k in 0:4) {
      if (dir.exists(file.path(q, "imagenes"))) { dir_raiz <- q; break }
      p <- dirname(q)
      if (identical(p, q)) break
      q <- p
    }
    if (is.null(dir_raiz)) {
      stop("No se encontró la carpeta del proyecto (carpeta 'imagenes/').\n",
           "Defina dir_raiz a mano al inicio de config_linux.R.")
    }
  }
  rm(.ruta_cfg)
}

dir_code <- file.path(dir_raiz, "Code")
dir_img  <- file.path(dir_raiz, "imagenes")
dir_out  <- file.path(dir_code, "resultados")

# Validación: las imágenes deben estar donde se esperan
if (!file.exists(file.path(dir_img, "horse_binaria.png"))) {
  stop("No se encontró 'horse_binaria.png' en: ", dir_img,
       "\nRevise dir_raiz en config_linux.R.")
}
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

cat("Rutas configuradas (Linux)\n",
    "  dir_raiz: ", dir_raiz, "\n",
    "  dir_img:  ", dir_img, "\n",
    "  dir_out:  ", dir_out, "\n", sep = "")
