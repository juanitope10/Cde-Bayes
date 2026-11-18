# Cde-Bayes

Restauración bayesiana de imágenes — caso de estudio hasta el punto 3.3
(`3.1.R` representación de imágenes y ruido, `3.2.R` modelo de Ising,
`3.3.R` GMRF en escala de grises). Documentación en `explicacion_completa.qmd`
(renderizar con `quarto render explicacion_completa.qmd`).

## Rutas (Linux y Windows)

Las rutas se definen una sola vez por sistema operativo, dentro de `Code/`:

- `Code/config_linux.R` — para Linux
- `Code/config_windows.R` — para Windows

Cada script carga automáticamente el archivo de configuración de su sistema
operativo y usa las variables `dir_raiz`, `dir_code`, `dir_img` y `dir_out`.
La detección es automática; solo si fallara, abra el archivo de configuración
correspondiente y defina `dir_raiz` en la línea marcada `[EDITAR]`, por ejemplo:

```r
dir_raiz <- "/home/usuario/Desktop/Cde bayes"    # Linux
dir_raiz <- "C:/Users/usuario/Desktop/Cde bayes" # Windows
```

## Cómo ejecutar

Desde la carpeta `Code/`:

```bash
Rscript 3.1.R   # exploración, grafo y verificaciones (~10 s)
Rscript 3.2.R   # modelo de Ising, 3 niveles de ruido (~30 min)
Rscript 3.3.R   # GMRF en escala de grises (~25 min)
```

En R/RStudio también sirve `source("Code/3.2.R")` o el botón *Source*. Paquetes
necesarios: `png`, `Matrix`, `ggplot2`, `coda`.

Variables de entorno opcionales:

| Variable | Efecto |
|:--|:--|
| `CDE_RAPIDO=TRUE` | Cadenas cortas de prueba (valida el código en ~2 min) |
| `CDE_LOGDET=espectral` | Log-determinante exacto y más rápido (`3.3.R`) |

Todas las salidas se guardan en `Code/resultados/`.
