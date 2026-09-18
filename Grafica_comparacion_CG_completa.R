# ============================================================
# CORRELACIÓN DE SPEARMAN ENTRE GC Y DENSIDAD GÉNICA
# POR CROMOSOMA + SCATTER PLOT
#
# Igual que el script de las figuras por cromosoma: reconoce
# automáticamente el tamaño de ventana (20kb, 50kb, 100kb, etc.)
# a partir del nombre del archivo de entrada, y usa eso para:
#   1) los títulos de las gráficas
#   2) los nombres de los archivos de salida
#   3) una carpeta de salida separada por ventana, para que las
#      corridas con distintos archivos no se pisen entre sí
# ============================================================

library(ggplot2)


# ============================================================
# 0. ARCHIVO DE ENTRADA
#    Para correr con otro tamaño de ventana, cambiá SOLO esta
#    línea (por ejemplo a "tabla_20kb.tsv" o "tabla_100kb.tsv").
#    Todo lo demás se ajusta solo.
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
# 0ter. CARPETA DE SALIDA (misma convención que el script de
#    las figuras por cromosoma: una carpeta por tamaño de
#    ventana, así todo lo generado para un mismo archivo de
#    entrada queda junto)
# ============================================================

carpeta_base   <- "~/Documentos/Tarea1/GCA_033216535.1/genes/salidaR/"
carpeta_salida <- paste0(carpeta_base, "figuras_cromosomas_", ventana_label, "/")

dir.create(carpeta_salida, showWarnings = FALSE, recursive = TRUE)


# ------------------------------------------------------------
# 1. Preparar datos
# ------------------------------------------------------------

tabla_cor <- tabla[
  is.finite(tabla$GC_percent) &
    is.finite(tabla$genes_per_Mb) &
    !is.na(tabla$chrom) &
    tabla$chrom != "",
]


# ------------------------------------------------------------
# 2. Calcular Spearman para cada cromosoma
# ------------------------------------------------------------

resultados_spearman <- do.call(
  rbind,
  lapply(
    split(tabla_cor, tabla_cor$chrom),
    function(datos_chr) {
      
      # Correlación
      prueba <- cor.test(
        datos_chr$GC_percent,
        datos_chr$genes_per_Mb,
        method = "spearman",
        exact = FALSE
      )
      
      data.frame(
        chrom = unique(datos_chr$chrom),
        n_ventanas = nrow(datos_chr),
        rho = unname(prueba$estimate),
        p_value = prueba$p.value
      )
    }
  )
)


# ------------------------------------------------------------
# 3. Corregir por múltiples comparaciones
# ------------------------------------------------------------

resultados_spearman$p_ajustado <-
  p.adjust(
    resultados_spearman$p_value,
    method = "BH"
  )


# ------------------------------------------------------------
# 4. Añadir interpretación de la dirección
# ------------------------------------------------------------

resultados_spearman$direccion <-
  
  ifelse(
    resultados_spearman$rho > 0,
    "Positiva",
    ifelse(
      resultados_spearman$rho < 0,
      "Negativa",
      "Sin asociación"
    )
  )


# ------------------------------------------------------------
# 5. Ordenar por cromosoma
# ------------------------------------------------------------

resultados_spearman <-
  resultados_spearman[
    order(resultados_spearman$chrom),
  ]


# ------------------------------------------------------------
# 6. Mostrar y guardar resultados
# ------------------------------------------------------------

print(
  resultados_spearman
)

write.csv(
  resultados_spearman,
  file      = paste0(carpeta_salida, "Spearman_GC_gene_density_", ventana_label, ".csv"),
  row.names = FALSE
)


# ============================================================
# GRÁFICA DE CORRELACIÓN DE SPEARMAN POR CROMOSOMA
# ============================================================

# Ordenar cromosomas por rho
resultados_spearman$chrom <- factor(
  resultados_spearman$chrom,
  levels = resultados_spearman$chrom[
    order(resultados_spearman$rho)
  ]
)

p_spearman <- ggplot(
  resultados_spearman,
  aes(x = rho, y = chrom)
) +
  
  # Línea de referencia: ausencia de correlación
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey55",
    linewidth = 0.55
  ) +
  
  # Línea desde 0 hasta rho
  geom_segment(
    aes(
      x = 0,
      xend = rho,
      y = chrom,
      yend = chrom
    ),
    color = "grey65",
    linewidth = 0.9
  ) +
  
  # Punto de rho
  geom_point(
    aes(color = rho),
    size = 4.5
  ) +
  
  # Valor de rho
  geom_text(
    aes(
      label = sprintf("%.2f", rho),
      hjust = ifelse(rho >= 0, -0.4, 1.4)
    ),
    size = 3.2,
    color = "grey20",
    family = "Helvetica"
  ) +
  
  scale_color_gradient2(
    low = "#2166AC",
    mid = "grey80",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    guide = "none"
  ) +
  
  scale_x_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.2),
    expand = expansion(mult = c(0.04, 0.10))
  ) +
  
  labs(
    title = "Spearman correlation between GC content and gene density",
    subtitle = paste0(ventana_label, " windows"),
    x = expression("Spearman's " * rho),
    y = NULL
  ) +
  
  theme_minimal(
    base_size = 10,
    base_family = "Helvetica"
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 13,
      margin = margin(b = 4)
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10,
      color = "grey40",
      margin = margin(b = 14)
    ),
    
    axis.text.y = element_text(
      color = "black",
      size = 9.5
    ),
    
    axis.text.x = element_text(
      color = "black",
      size = 9
    ),
    
    axis.title.x = element_text(
      color = "black",
      size = 10,
      margin = margin(t = 8)
    ),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    panel.grid.major.x = element_line(
      color = "grey92",
      linewidth = 0.4
    ),
    
    axis.ticks.y = element_blank(),
    
    axis.ticks.x = element_line(
      color = "grey40",
      linewidth = 0.3
    ),
    
    plot.margin = margin(18, 35, 12, 12)
  )

print(p_spearman)


# ============================================================
# GUARDAR GRÁFICA DE SPEARMAN EN PDF Y TIFF
# ============================================================

alto_spearman <- max(4, 0.42 * nrow(resultados_spearman) + 1.6)

ggsave(
  filename = paste0(carpeta_salida, "Spearman_GC_gene_density_", ventana_label, ".pdf"),
  plot     = p_spearman,
  width    = 8,
  height   = alto_spearman,
  units    = "in",
  device   = cairo_pdf
)

ggsave(
  filename    = paste0(carpeta_salida, "Spearman_GC_gene_density_", ventana_label, ".tiff"),
  plot        = p_spearman,
  width       = 8,
  height      = alto_spearman,
  units       = "in",
  dpi         = 600,
  compression = "lzw"
)


# ============================================================
# SCATTER PLOT: GC vs DENSIDAD GÉNICA
# CADA PUNTO = UNA VENTANA
# UN PANEL POR CROMOSOMA
# ============================================================

tabla_scatter <- tabla[
  is.finite(tabla$GC_percent) &
    is.finite(tabla$genes_per_Mb) &
    !is.na(tabla$chrom) &
    tabla$chrom != "",
]

# Ordenar cromosomas
tabla_scatter$chrom <- factor(
  tabla_scatter$chrom,
  levels = unique(tabla_scatter$chrom)
)

p_scatter <- ggplot(
  tabla_scatter,
  aes(
    x = GC_percent,
    y = genes_per_Mb
  )
) +
  
  # Cada ventana
  geom_point(
    size = 1.5,
    alpha = 0.65
  ) +
  
  # Tendencia lineal
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.7
  ) +
  
  facet_wrap(
    ~ chrom,
    scales = "free_y",
    ncol = 1
  ) +
  
  labs(
    title = paste0("GC content and gene density across ", ventana_label, " windows"),
    x = "GC content (%)",
    y = "Gene density (genes/Mb)"
  ) +
  
  theme_minimal(
    base_size = 10,
    base_family = "Helvetica"
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 13,
      margin = margin(b = 14)
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 9,
      color = "black"
    ),
    
    strip.background = element_blank(),
    
    axis.text = element_text(
      color = "black",
      size = 8
    ),
    
    axis.title = element_text(
      color = "black",
      size = 10
    ),
    
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.35
    ),
    
    panel.grid.minor = element_blank(),
    
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.3
    ),
    
    panel.border = element_rect(
      color = "grey70",
      fill = NA,
      linewidth = 0.4
    ),
    
    plot.margin = margin(10, 15, 10, 15)
  )

print(p_scatter)


# ============================================================
# GUARDAR SCATTER PLOT EN PDF Y TIFF
# ============================================================

ggsave(
  filename = paste0(carpeta_salida, "Scatter_GC_gene_density_", ventana_label, ".pdf"),
  plot     = p_scatter,
  width    = 8,
  height   = max(6, 2.8 * length(unique(tabla_scatter$chrom))),
  units    = "in",
  device   = cairo_pdf
)

ggsave(
  filename    = paste0(carpeta_salida, "Scatter_GC_gene_density_", ventana_label, ".tiff"),
  plot        = p_scatter,
  width       = 8,
  height      = max(6, 2.8 * length(unique(tabla_scatter$chrom))),
  units       = "in",
  dpi         = 600,
  compression = "lzw"
)


# ============================================================
# GUARDAR LA TABLA DE ENTRADA USADA (nombre con la ventana,
# no más un "20kb" hardcodeado que no correspondía)
# ============================================================

write.table(
  tabla,
  file.path(carpeta_salida, paste0("tabla_usada_", ventana_label, ".tsv")),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)
