# Cde-Bayes

Restauración bayesiana de imágenes — caso de estudio completo (3.1–3.4 y discusión):
`3.1.R` representación de imágenes y ruido, `3.2.R` modelo de Ising,
`3.3.R` GMRF en escala de grises, `3.4.R` GMRF multivariado RGB y `4.R` tabla
comparativa y discusión. Documentación en `fuentes/explicacion_completa.qmd` (renderizar
con `quarto render fuentes/explicacion_completa.qmd` y copiar el HTML a la raíz).

## Estructura

- Raíz: los HTML renderizados (`caso.html`, `explicacion_completa.html`,
  `gmrf_grises.html`, `gmrf_rgb.html`, `ising.html`, `representacion_y_discusion.html`).
- `fuentes/`: documentos fuente (`.qmd` y `.Rmd`), `Cuadernos/`, `Code/`,
  `imagenes/` y las salidas de apoyo (`corrida/`, `figuras/`).

## Rutas (Linux y Windows)

Las rutas se definen una sola vez por sistema operativo, dentro de `fuentes/Code/`:

- `fuentes/Code/config_linux.R` — para Linux
- `fuentes/Code/config_windows.R` — para Windows

Cada script carga automáticamente el archivo de configuración de su sistema
operativo y usa las variables `dir_raiz`, `dir_code`, `dir_img` y `dir_out`.
La detección es automática; solo si fallara, abra el archivo de configuración
correspondiente y defina `dir_raiz` en la línea marcada `[EDITAR]`, por ejemplo:

```r
dir_raiz <- "/home/usuario/Desktop/cde bayescompleto"    # Linux
dir_raiz <- "C:/Users/usuario/Desktop/cde bayescompleto" # Windows
```

Cada script deja sus salidas en su propia subcarpeta de `fuentes/Code/resultados/`:
`3.2/`, `3.3/`, `3.4/`, `convergencia/`, `verificacion/`, `sensibilidad/` y `4/`.

## Cómo ejecutar

Desde la carpeta `fuentes/Code/`:

```bash
Rscript 3.1.R   # exploración, grafo y verificaciones (~10 s)
Rscript 3.2.R   # modelo de Ising, 3 niveles de ruido (~30 min)
Rscript 3.3.R   # GMRF en escala de grises (~25 min con CDE_LOGDET=espectral)
Rscript 3.4.R   # GMRF RGB (~60 min con CDE_LOGDET=espectral)
Rscript convergencia.R   # tests de convergencia de las trazas guardadas (~2 min)
Rscript verificaciones.R # log-verosimilitud/DIC, predictivas y contrastes (~30 min)
Rscript sensibilidad.R   # 5 variantes de previa con DIC/WAIC (~40 min con espectral)
Rscript 4.R     # tabla comparativa global + LaTeX + discusión
```

En R/RStudio también sirve `source("fuentes/Code/3.2.R")` o el botón *Source*. Paquetes
necesarios: `png`, `Matrix`, `ggplot2`, `coda`.

Variables de entorno opcionales:

| Variable | Efecto |
|:--|:--|
| `CDE_RAPIDO=TRUE` | Cadenas cortas de prueba (valida el código en ~2 min) |
| `CDE_LOGDET=espectral` | Log-determinante exacto y más rápido (`3.3.R`, `3.4.R`) |
