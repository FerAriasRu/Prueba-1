# ============================================================
# GC (%) + DENSIDAD GÉNICA POR CROMOSOMA
# FIGURA FINAL EN GRILLA, ESTILO PUBLICACIÓN
# (mismo lenguaje visual que el lollipop de GC global:
#  Helvetica, ejes finos, título en negrita centrado,
#  mismos colores azul/rojo)
# ============================================================

library(ggplot2)
library(patchwork)
library(grid)
library(gtable)


# ============================================================
# 0. TEMA BASE "TIPO PUBLICACIÓN"
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
# 1. CROMOSOMAS Y CONFIGURACIÓN DE LA GRILLA
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
# 2. CREAR CARPETA DE SALIDA
# ============================================================

dir.create("figuras_cromosomas", showWarnings = FALSE)


# ============================================================
# 3. LONGITUD DE CADA CROMOSOMA
# ============================================================

longitudes <- aggregate(end ~ chrom, data = tabla, FUN = max)
colnames(longitudes)[2] <- "length_bp"
longitudes$length_Mb <- longitudes$length_bp / 1000000

max_length_Mb <- max(longitudes$length_Mb)


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
  # ESCALA DE GC
  # ==========================================================
  
  min_GC <- floor(min(datos_chr$GC_percent, na.rm = TRUE))
  max_GC <- ceiling(max(datos_chr$GC_percent, na.rm = TRUE))
  breaks_GC <- c(min_GC, max_GC)
  
  # ==========================================================
  # ESCALA DE DENSIDAD GÉNICA
  # ==========================================================
  
  max_genes <- max(datos_chr$genes_per_Mb, na.rm = TRUE)
  max_genes_break <- ceiling(max_genes / 100) * 100
  breaks_genes <- seq(0, max_genes_break, by = 100)
  
  # ==========================================================
  # GRÁFICA SUPERIOR: GC
  # ==========================================================
  
  p_GC <- ggplot(datos_chr, aes(x = position_Mb, y = GC_percent)) +
    geom_line(linewidth = 0.6, color = color_GC) +
    labs(title = titulo_chr, x = NULL, y = NULL) +
    scale_x_continuous(limits = c(0, max_length_Mb), expand = c(0, 0)) +
    scale_y_continuous(
      breaks = breaks_GC,
      limits = c(min_GC, max_GC),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    tema_publicacion +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin  = margin(2, 6, 0, 6)
    )
  
  # ==========================================================
  # GRÁFICA INFERIOR: DENSIDAD GÉNICA
  #    (solo la última fila de la grilla muestra los números
  #     del eje X, igual que en la imagen de referencia)
  # ==========================================================
  
  p_genes <- ggplot(datos_chr, aes(x = position_Mb, y = genes_per_Mb)) +
    geom_line(linewidth = 0.6, color = color_genes) +
    labs(x = NULL, y = NULL) +
    scale_x_continuous(limits = c(0, max_length_Mb), expand = c(0, 0)) +
    scale_y_continuous(
      breaks = breaks_genes,
      limits = c(0, max_genes_break),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    tema_publicacion +
    theme(
      axis.text.x  = if (es_ultima_fila) element_text(size = 8) else element_blank(),
      axis.ticks.x = if (es_ultima_fila) element_line(color = "grey40", linewidth = 0.3) else element_blank(),
      plot.margin  = margin(0, 6, 2, 6)
    )
  
  # ==========================================================
  # UNIR GC + DENSIDAD (por cromosoma)
  # ==========================================================
  
  figura_chr <- p_GC / p_genes + plot_layout(heights = c(1, 1))
  figuras_final[[chr]] <- figura_chr
  
  print(figura_chr)
  
  # ==========================================================
  # GUARDAR FIGURA INDIVIDUAL
  # ==========================================================
  
  ggsave(
    filename = paste0("figuras_cromosomas/", chr, "_GC_gene_density.pdf"),
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
# ============================================================

figura_completa <-
  (wrap_elements(full = figura_todos) / eje_x_texto / wrap_elements(full = legend_grob)) +
  plot_layout(heights = c(n_filas, 0.15, 0.35)) &
  theme(plot.margin = margin(5, 15, 5, 15))

print(figura_completa)


# ============================================================
# 9. GUARDAR PDF FINAL (ancho/alto se ajustan solos a la grilla)
# ============================================================

ancho_fig <- 8                  # ancho fijo, legible con una sola columna
alto_fig  <- 2.6 * n_filas + 1.2  # alto proporcional a la cantidad de cromosomas

ggsave(
  filename = "figuras_cromosomas/ALL_chromosomes_GC_gene_density.pdf",
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
  filename    = "figuras_cromosomas/ALL_chromosomes_GC_gene_density.tiff",
  plot        = figura_completa,
  width       = ancho_fig,
  height      = alto_fig,
  units       = "in",
  dpi         = 600,
  compression = "lzw"
)
