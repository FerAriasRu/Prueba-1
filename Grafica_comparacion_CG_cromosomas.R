# ============================================================
# COMPARACIÓN DEL GC (%) GLOBAL ENTRE CROMOSOMAS
# GRÁFICO DE PUNTOS + LÍNEA (ORDENADO), ESTILO PUBLICACIÓN
# ============================================================
#
# Requiere el paquete ggplot2. Si no lo tenés instalado, correr
# UNA sola vez (no hace falta repetirlo cada sesión):
#
#   install.packages("ggplot2")
#
# ============================================================
setwd("/home/luisafernanda/Documentos/Tarea1")
library(ggplot2)


# ============================================================
# 1. CREAR CARPETA DE SALIDA
# ============================================================

dir.create("figuras_cromosomas", showWarnings = FALSE)


# ============================================================
# 2. CALCULAR GC GLOBAL POR CROMOSOMA
#
#    Se pondera por el largo de cada ventana (end - start),
#    en vez de un promedio simple, para que ventanas más largas
#    pesen más en el GC global del cromosoma (más preciso que
#    un mean() simple si las ventanas no son todas del mismo
#    tamaño).
# ============================================================

tabla$window_length <- tabla$end - tabla$start

GC_por_cromosoma <- aggregate(
  cbind(GC_weighted = GC_percent * window_length, total_length = window_length)
  ~ chrom,
  data = tabla,
  FUN  = sum
)

GC_por_cromosoma$GC_global <-
  GC_por_cromosoma$GC_weighted / GC_por_cromosoma$total_length

GC_por_cromosoma <- GC_por_cromosoma[, c("chrom", "GC_global")]


# ============================================================
# 3. ORDENAR CROMOSOMAS DE MAYOR A MENOR GC
# ============================================================

GC_por_cromosoma <- GC_por_cromosoma[
  order(-GC_por_cromosoma$GC_global),
]

GC_por_cromosoma$chrom <- factor(
  GC_por_cromosoma$chrom,
  levels = GC_por_cromosoma$chrom
)


# ============================================================
# 4. PROMEDIO GLOBAL DEL GENOMA (línea de referencia)
# ============================================================

GC_promedio_genoma <- weighted.mean(
  tabla$GC_percent,
  w = tabla$window_length
)


# ============================================================
# 5. LÍMITES DEL EJE Y ("ZOOM" A LA VARIACIÓN REAL)
# ============================================================

rango_GC <- range(GC_por_cromosoma$GC_global)
margen    <- max(0.3, diff(rango_GC) * 0.25)

y_min <- rango_GC[1] - margen
y_max <- rango_GC[2] + margen


# ============================================================
# 6. ORDEN PARA EL EJE Y
#    (el cromosoma con mayor GC arriba del todo, por eso se
#    invierte el orden de los niveles del factor)
# ============================================================

GC_por_cromosoma$chrom <- factor(
  GC_por_cromosoma$chrom,
  levels = rev(levels(GC_por_cromosoma$chrom))
)

# Posición horizontal de las etiquetas de texto:
# a la derecha del punto si el valor es mayor al promedio,
# a la izquierda si es menor (para que nunca se pisen con el
# segmento ni con la línea del promedio)

GC_por_cromosoma$label_hjust <- ifelse(
  GC_por_cromosoma$GC_global >= GC_promedio_genoma,
  -0.35,
  1.35
)


# ============================================================
# 7. GRÁFICO: LOLLIPOP HORIZONTAL CON GRADIENTE DE COLOR
# ============================================================

p_GC_global <- ggplot(
  GC_por_cromosoma,
  aes(x = GC_global, y = chrom)
) +
  
  # Línea vertical de referencia: promedio del genoma
  geom_vline(
    xintercept = GC_promedio_genoma,
    linetype   = "dashed",
    color      = "grey55",
    linewidth  = 0.55
  ) +
  
  # "Palito" del lollipop: desde el promedio hasta el valor
  geom_segment(
    aes(x = GC_promedio_genoma, xend = GC_global, yend = chrom),
    color     = "grey65",
    linewidth = 0.9
  ) +
  
  # Punto coloreado según GC (gradiente diverging centrado en el promedio)
  geom_point(
    aes(color = GC_global),
    size = 4.6
  ) +
  
  scale_color_gradient2(
    low      = "#2166AC",
    mid      = "grey80",
    high     = "#B2182B",
    midpoint = GC_promedio_genoma,
    guide    = "none"
  ) +
  
  # Etiqueta con el valor exacto, a un costado del punto
  geom_text(
    aes(label = sprintf("%.2f%%", GC_global), hjust = label_hjust),
    size   = 3.1,
    color  = "grey20",
    family = "Helvetica"
  ) +
  
  # Anotación del promedio del genoma (arriba de la línea punteada)
  annotate(
    "text",
    x      = GC_promedio_genoma,
    y      = Inf,
    label  = paste0("Genome mean: ", sprintf("%.2f", GC_promedio_genoma), "%"),
    vjust  = -0.9,
    size   = 3.1,
    color  = "grey35",
    family = "Helvetica",
    fontface = "italic"
  ) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    title = "Global GC content comparison across chromosomes",
    x     = "GC content (%)",
    y     = NULL
  ) +
  
  scale_x_continuous(
    limits = c(y_min, y_max),
    expand = expansion(mult = c(0.06, 0.06))
  ) +
  
  theme_minimal(base_size = 11, base_family = "Helvetica") +
  
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 13,
                                    margin = margin(b = 14, t = 6)),
    axis.text.y      = element_text(color = "black", size = 9.5,
                                    face = "plain", margin = margin(r = 4)),
    axis.text.x      = element_text(color = "black", size = 9),
    axis.title.x     = element_text(color = "black", size = 10,
                                    margin = margin(t = 8)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.4),
    axis.ticks.y     = element_blank(),
    axis.ticks.x     = element_line(color = "grey40", linewidth = 0.3),
    plot.margin      = margin(18, 30, 12, 12)
  )

print(p_GC_global)

# Alto de la figura proporcional a la cantidad de cromosomas
altura_fig <- max(4, 0.42 * nrow(GC_por_cromosoma) + 1.6)


# ============================================================
# 7. GUARDAR PDF (vectorial, ideal para figuras de manuscrito)
# ============================================================

ggsave(
  filename = "figuras_cromosomas/GC_global_por_cromosoma.pdf",
  plot     = p_GC_global,
  width    = 8,
  height   = altura_fig,
  units    = "in",
  device   = cairo_pdf
)


# ============================================================
# 8. GUARDAR TIFF (600 dpi, formato típico exigido por revistas)
# ============================================================

ggsave(
  filename    = "figuras_cromosomas/GC_global_por_cromosoma.tiff",
  plot        = p_GC_global,
  width       = 8,
  height      = altura_fig,
  units       = "in",
  dpi         = 600,
  compression = "lzw"
)


# ============================================================
# 9. (OPCIONAL) EXPORTAR TABLA CON LOS VALORES USADOS
# ============================================================

write.csv(
  GC_por_cromosoma,
  file      = "figuras_cromosomas/GC_global_por_cromosoma.csv",
  row.names = FALSE
)
