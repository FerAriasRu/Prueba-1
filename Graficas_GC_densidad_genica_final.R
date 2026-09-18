# ============================================================
# GC (%) + DENSIDAD GÉNICA POR CROMOSOMA
# FIGURA FINAL EN GRILLA, ESTILO PUBLICACIÓN
# (mismo lenguaje visual que el lollipop de GC global:
#  Helvetica, ejes finos, título en negrita centrado,
#  mismos colores azul/rojo)
#
# Este script reconoce automáticamente el tamaño de ventana
# (20kb, 50kb, 100kb, etc.) a partir del nombre del archivo de
# entrada (o, si no lo puede leer del nombre, lo calcula de los
# propios datos). Ese tamaño de ventana se usa para:
#   1) el título de la figura final
#   2) una carpeta de salida separada por ventana, para que las
#      corridas con distintos archivos no se pisen entre sí
# ============================================================

library(ggplot2)
library(patchwork)
library(grid)
library(gtable)


# ============================================================
# 0. ARCHIVO DE ENTRADA
#    Para correr con otro tamaño de ventana, cambiá SOLO esta
#    línea (por ejemplo a "tabla_20kb.tsv" o "tabla_100kb.tsv").
#    Todo lo demás (carpeta de salida, título de la figura) se
#    ajusta solo.
# ============================================================

archivo_tabla <- "~/Documentos/Tarea1/GCA_033216535.1/genes/tabla_50kb.tsv"

tabla <- read.delim(archivo_tabla, header = TRUE)


# ============================================================
# 0bis. DETECTAR EL TAMAÑO DE VENTANA
#    Primero intenta leerlo del nombre del archivo (busca un
#    patrón tipo "50kb", "20Kb", "100KB", etc.). Si no lo
#    encuentra, lo calcula directamente de los datos (tamaño de
#    la primera ventana = end - start).
# ============================================================

nombre_archivo <- basename(archivo_tabla)
match_ventana  <- regmatches(
  nombre_archivo,
  regexpr("[0-9]+ *[kKmM][bB]", nombre_archivo)
)

if (length(match_ventana) == 1) {
  # Normaliza el formato: "50Kb", "50 kb", "50KB" -> "50kb"
  numero_ventana <- regmatches(match_ventana, regexpr("[0-9]+", match_ventana))
  unidad_ventana <- tolower(gsub("[0-9 ]", "", match_ventana))
  ventana_label  <- paste0(numero_ventana, unidad_ventana)
} else {
  # No se pudo leer del nombre -> se calcula de los datos
  ventana_bp <- (tabla$end - tabla$start)[1]
  ventana_label <- if (ventana_bp >= 1e6) {
    paste0(ventana_bp / 1e6, "mb")
  } else {
    paste0(ventana_bp / 1e3, "kb")
  }
}

cat("Tamaño de ventana detectado:", ventana_label, "\n")


# ============================================================
# 0ter. CARPETA DE SALIDA (una por tamaño de ventana)
# ============================================================

carpeta_base   <- "~/Documentos/Tarea1/GCA_033216535.1/genes/salidaR/"
carpeta_salida <- paste0(carpeta_base, "figuras_cromosomas_", ventana_label, "/")

dir.create(carpeta_salida, showWarnings = FALSE, recursive = TRUE)


# ============================================================
# 1. TEMA BASE "TIPO PUBLICACIÓN"
#    Se reutiliza en todos los paneles para que la figura sea
#    visualmente consistente con el gráfico de GC global.
# ============================================================

font_family <- "Helvetica"   # cambiar a "Arial" si está disponible en el sistema

# Mismos colores que el gradiente del lollipop (low / high)
color_GC    <- "#2166AC"   # azul
color_genes <- "#B2182B"   # rojo

tema_publicacion <- theme_classic(base_size = 10, base_family = font_family) +
  theme(
    axis.line          = element_line(color = "grey20", linewidth = 0.4),
    axis.ticks         = element_line(color = "grey40", linewidth = 0.3),
    axis.text          = element_text(color = "black", size = 8),
    plot.title         = element_text(face = "bold", hjust = 0.5, size = 11,
                                      margin = margin(b = 3)),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    legend.position    = "none"
  )


# ============================================================
# 2. CROMOSOMAS Y CONFIGURACIÓN DE LA GRILLA
# ============================================================

cromosomas <- unique(tabla$chrom)
n_chr      <- length(cromosomas)

# Número de columnas de la grilla final -> 1 columna, todos los
# cromosomas apilados verticalmente, con un único eje X al pie
# (igual que en la imagen de referencia)
ncol_grid <- 1
n_filas   <- ceiling(n_chr / ncol_grid)

# Fila de la grilla en la que cae cada cromosoma, para saber
# cuáles van en la última fila (solo esos muestran los números
# del eje X; el título "Genomic position (Mb)" va una sola vez
# al pie de toda la figura, como en la imagen de referencia)
fila_de           <- ceiling(seq_along(cromosomas) / ncol_grid)
ultima_fila_chrs  <- cromosomas[fila_de == n_filas]


# ============================================================
# 3. LONGITUD DE CADA CROMOSOMA
# ============================================================

longitudes <- aggregate(end ~ chrom, data = tabla, FUN = max)
colnames(longitudes)[2] <- "length_bp"
longitudes$length_Mb <- longitudes$length_bp / 1000000

max_length_Mb <- max(longitudes$length_Mb)


# ============================================================
# 3bis. ESCALAS GLOBALES (ESTANDARIZADAS PARA TODOS LOS PANELES)
#    Se calculan UNA sola vez, a partir de TODA la tabla, para
#    que el eje Y de GC (y el de densidad génica) sea el mismo
#    en todos los cromosomas y se puedan comparar entre sí a
#    simple vista, sin que cada uno tenga su propia escala.
# ============================================================

min_GC_global <- floor(min(tabla$GC_percent, na.rm = TRUE))
max_GC_global <- ceiling(max(tabla$GC_percent, na.rm = TRUE))
breaks_GC_global <- pretty(c(min_GC_global, max_GC_global), n = 4)
breaks_GC_global <- breaks_GC_global[breaks_GC_global >= min_GC_global &
                                       breaks_GC_global <= max_GC_global]

max_genes_global       <- max(tabla$genes_per_Mb, na.rm = TRUE)
max_genes_break_global <- ceiling(max_genes_global / 100) * 100
breaks_genes_global    <- seq(0, max_genes_break_global, by = 100)


# ============================================================
# 4. CONSTRUIR CADA PANEL (GC arriba + densidad génica abajo)
# ============================================================

figuras_final <- list()

for (chr in cromosomas) {
  
  # ----------------------------------------------------------
  # DATOS DEL CROMOSOMA
  # ----------------------------------------------------------
  
  datos_chr <- tabla[tabla$chrom == chr, ]
  datos_chr$position_Mb <- datos_chr$start / 1000000
  
  longitud_chr <- longitudes$length_Mb[longitudes$chrom == chr]
  titulo_chr   <- paste0(chr, " (", sprintf("%.2f", longitud_chr), " Mb)")
  
  es_ultima_fila <- chr %in% ultima_fila_chrs
  
  # ==========================================================
  # GRÁFICA SUPERIOR: GC
  #    Usa la escala global (min_GC_global / max_GC_global /
  #    breaks_GC_global) en vez de una escala propia por
  #    cromosoma, para que todos los paneles sean comparables.
  # ==========================================================
  
  p_GC <- ggplot(datos_chr, aes(x = position_Mb, y = GC_percent)) +
    geom_line(linewidth = 0.6, color = color_GC) +
    labs(title = titulo_chr, x = NULL, y = NULL) +
    scale_x_continuous(limits = c(0, max_length_Mb), expand = c(0, 0)) +
    scale_y_continuous(
      breaks = breaks_GC_global,
      limits = c(min_GC_global, max_GC_global),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    tema_publicacion +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin  = margin(t = 2, r = 6, b = 10, l = 6)
    )
  
  # ==========================================================
  # GRÁFICA INFERIOR: DENSIDAD GÉNICA
  #    También usa la escala global (breaks_genes_global /
  #    max_genes_break_global), y solo la última fila de la
  #    grilla muestra los números del eje X.
  # ==========================================================
  
  p_genes <- ggplot(datos_chr, aes(x = position_Mb, y = genes_per_Mb)) +
    geom_line(linewidth = 0.6, color = color_genes) +
    labs(x = NULL, y = NULL) +
    scale_x_continuous(limits = c(0, max_length_Mb), expand = c(0, 0)) +
    scale_y_continuous(
      breaks = breaks_genes_global,
      limits = c(0, max_genes_break_global),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    tema_publicacion +
    theme(
      axis.text.x  = if (es_ultima_fila) element_text(size = 8) else element_blank(),
      axis.ticks.x = if (es_ultima_fila) element_line(color = "grey40", linewidth = 0.3) else element_blank(),
      plot.margin  = margin(t = 6, r = 6, b = 2, l = 6)
    )
  
  # ==========================================================
  # UNIR GC + DENSIDAD (por cromosoma)
  # ==========================================================
  
  figura_chr <- p_GC / p_genes + plot_layout(heights = c(1, 1))
  figuras_final[[chr]] <- figura_chr
  
  print(figura_chr)
  
  # ==========================================================
  # GUARDAR FIGURA INDIVIDUAL
  #    El nombre de archivo incluye la ventana (ej. "_50kb")
  #    para que no se pisen entre corridas con otros tamaños.
  # ==========================================================
  
  ggsave(
    filename = paste0(carpeta_salida, chr, "_GC_gene_density_", ventana_label, ".pdf"),
    plot     = figura_chr,
    width    = 7,
    height   = 5,
    units    = "in",
    device   = cairo_pdf
  )
}


# ============================================================
# 5. FIGURA CON TODOS LOS CROMOSOMAS, EN GRILLA
# ============================================================

figura_todos <- wrap_plots(figuras_final, ncol = ncol_grid)


# ============================================================
# 6. TÍTULO ÚNICO DEL EJE X, AL PIE DE TODA LA FIGURA
# ============================================================

eje_x_texto <- wrap_elements(full = textGrob(
  "Genomic position (Mb)",
  gp = gpar(fontfamily = font_family, fontsize = 10, col = "black")
))


# ============================================================
# 7. CONSTRUIR UNA ÚNICA LEYENDA (GC + DENSIDAD GÉNICA)
# ============================================================

legend_plot <- ggplot(
  data.frame(
    x    = c(1, 2),
    y    = c(1, 1),
    tipo = factor(
      c("GC content (%)", "Gene density (genes/Mb)"),
      levels = c("GC content (%)", "Gene density (genes/Mb)")
    )
  ),
  aes(x = x, y = y, color = tipo)
) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(
    values = c(
      "GC content (%)"          = color_GC,
      "Gene density (genes/Mb)" = color_genes
    )
  ) +
  labs(color = NULL) +
  theme_void(base_family = font_family) +
  theme(
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.text        = element_text(size = 10, color = "black"),
    legend.key.width   = unit(1.3, "cm"),
    legend.background  = element_rect(color = "black", fill = "white", linewidth = 0.4),
    legend.margin      = margin(6, 14, 6, 14),
    legend.key         = element_blank()
  )

# Extraer solo el grob de la leyenda

g <- ggplotGrob(legend_plot)
legend_index <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
legend_grob  <- g$grobs[[legend_index]]


# ============================================================
# 8. FIGURA FINAL: GRILLA DE CROMOSOMAS + EJE X + LEYENDA
#    El título general incluye el tamaño de ventana detectado
#    (ej. "... (50kb windows)"), para identificar de un vistazo
#    con qué archivo se generó la figura.
# ============================================================

titulo_general <- paste0(
  "GC content and gene density per chromosome (", ventana_label, " windows)"
)

figura_completa <-
  (wrap_elements(full = figura_todos) / eje_x_texto / wrap_elements(full = legend_grob)) +
  plot_layout(heights = c(n_filas, 0.15, 0.35)) &
  theme(plot.margin = margin(5, 15, 5, 15))

figura_completa <- figura_completa +
  plot_annotation(
    title = titulo_general,
    theme = theme(
      plot.title = element_text(
        face   = "bold",
        hjust  = 0.5,
        size   = 14,
        family = font_family,
        margin = margin(b = 10)
      )
    )
  )

print(figura_completa)


# ============================================================
# 9. GUARDAR PDF FINAL (ancho/alto se ajustan solos a la grilla)
#    El nombre de archivo también incluye la ventana detectada.
# ============================================================

ancho_fig <- 8                    # ancho fijo, legible con una sola columna
alto_fig  <- 2.6 * n_filas + 1.4  # alto proporcional a la cantidad de cromosomas

ggsave(
  filename = paste0(carpeta_salida, "ALL_chromosomes_GC_gene_density_", ventana_label, ".pdf"),
  plot     = figura_completa,
  width    = ancho_fig,
  height   = alto_fig,
  units    = "in",
  device   = cairo_pdf
)


# ============================================================
# 10. GUARDAR TIFF FINAL (600 dpi, listo para publicación)
# ============================================================

ggsave(
  filename    = paste0(carpeta_salida, "ALL_chromosomes_GC_gene_density_", ventana_label, ".tiff"),
  plot        = figura_completa,
  width       = ancho_fig,
  height      = alto_fig,
  units       = "in",
  dpi         = 600,
  compression = "lzw"
)


# ============================================================
# 11. (OPCIONAL) GUARDAR TABLA RESUMEN CON LOS DATOS USADOS
#     Útil para comparar distintas ventanas más adelante.
# ============================================================

write.csv(
  tabla,
  file      = paste0(carpeta_salida, "tabla_usada_", ventana_label, ".csv"),
  row.names = FALSE
)
