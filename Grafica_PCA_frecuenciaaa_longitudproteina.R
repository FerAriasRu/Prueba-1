# ============================================================
# COMPOSICIÓN DE AMINOÁCIDOS + PCA A PARTIR DE UN FASTA DE
# PROTEÍNAS
#
# Este script:
#   1) Lee el FASTA de proteínas
#   2) Calcula la longitud de cada proteína
#   3) Calcula la frecuencia relativa de los 20 aminoácidos
#      estándar en cada proteína
#   4) Identifica proteínas con caracteres ambiguos (X, B, Z, J,
#      U, O, *, etc.) o con longitud inusualmente corta, y
#      documenta el criterio usado para cada filtro
#   5) Corre una PCA (una proteína por fila, las 20 frecuencias
#      relativas como variables) sobre las proteínas que pasan
#      el filtro de calidad
#   6) Guarda una figura de PCA de calidad de publicación, y
#      todas las tablas intermedias, para que quede documentado
#      qué se filtró y por qué
#
# Paquetes necesarios (instalar UNA sola vez si no los tenés):
#   install.packages(c("seqinr", "ggplot2"))
# ============================================================
#install.packages("seqinr")

library(seqinr)
library(ggplot2)


# ============================================================
# 0. ARCHIVO DE ENTRADA Y CARPETA DE SALIDA
# ============================================================

archivo_fasta  <- "/home/luisafernanda/Documentos/Tarea1/Proteinas/proteins.fasta"
carpeta_salida <- "/home/luisafernanda/Documentos/Tarea1/Proteinas/PCA_aminoacidos/"

dir.create(carpeta_salida, showWarnings = FALSE, recursive = TRUE)

font_family <- "Helvetica"
color_bajo  <- "#2166AC"   # mismo azul usado en las otras figuras
color_alto  <- "#B2182B"   # mismo rojo usado en las otras figuras


# ============================================================
# 1. CRITERIOS DE FILTRADO (documentados acá, en un solo lugar,
#    para que sea fácil de revisar y de ajustar)
# ============================================================

# Longitud mínima para considerar una proteína "inusualmente
# corta". 30 aa es un umbral convencional en trabajos de
# composición de aminoácidos/PCA: por debajo de eso, la
# frecuencia relativa de cada aminoácido está calculada sobre
# muy pocos residuos y es muy ruidosa (un solo residuo ya
# representa >3% de la secuencia).
LONGITUD_MINIMA <- 30

# Los 20 aminoácidos estándar (código de una letra)
AA_20 <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
           "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

# Caracteres considerados "ambiguos" / no estándar en una
# secuencia de proteína (además de cualquier símbolo que no
# esté en AA_20): X = desconocido, B = Asx (Asn o Asp),
# Z = Glx (Gln o Glu), J = Leu o Ile, U = Selenocisteína,
# O = Pirrolisina, * = codón de stop, - = gap
CARACTERES_AMBIGUOS <- c("X", "B", "Z", "J", "U", "O", "*", "-")


# ============================================================
# 2. LEER EL FASTA
# ============================================================

secuencias_raw <- read.fasta(
  archivo_fasta,
  seqtype    = "AA",
  as.string  = TRUE,
  forceDNAtolower = FALSE
)

proteinas <- data.frame(
  id       = names(secuencias_raw),
  sequence = toupper(sapply(secuencias_raw, function(x) as.character(x[1]))),
  stringsAsFactors = FALSE
)

# Quitar un "*" final (codón de stop) si está presente, para
# que no cuente como parte de la longitud ni de la composición
proteinas$sequence <- sub("\\*$", "", proteinas$sequence)

cat("Proteínas leídas del FASTA:", nrow(proteinas), "\n")


# ============================================================
# 3. LONGITUD DE CADA PROTEÍNA
# ============================================================

proteinas$length_aa <- nchar(proteinas$sequence)


# ============================================================
# 4. FRECUENCIA RELATIVA DE LOS 20 AMINOÁCIDOS POR PROTEÍNA
#    (conteo de cada aminoácido / longitud total de la
#    secuencia; si hay caracteres ambiguos, las 20 frecuencias
#    no van a sumar 1 -- eso es intencional, y es justamente lo
#    que se usa más abajo para detectarlos)
# ============================================================

calcular_frecuencias <- function(seq_string) {
  residuos <- strsplit(seq_string, "")[[1]]
  conteos  <- table(factor(residuos, levels = AA_20))
  as.numeric(conteos) / length(residuos)
}

matriz_frecuencias <- t(sapply(proteinas$sequence, calcular_frecuencias))
colnames(matriz_frecuencias) <- AA_20

proteinas <- cbind(proteinas, as.data.frame(matriz_frecuencias))


# ============================================================
# 5. IDENTIFICAR PROTEÍNAS CON CARACTERES AMBIGUOS
#    (cualquier letra que no sea una de las 20 estándar)
# ============================================================

detectar_ambiguos <- function(seq_string) {
  residuos  <- unique(strsplit(seq_string, "")[[1]])
  ambiguos  <- residuos[!(residuos %in% AA_20)]
  paste(sort(ambiguos), collapse = ",")
}

proteinas$caracteres_ambiguos <- sapply(proteinas$sequence, detectar_ambiguos)
proteinas$tiene_ambiguos      <- proteinas$caracteres_ambiguos != ""


# ============================================================
# 6. IDENTIFICAR PROTEÍNAS INUSUALMENTE CORTAS
# ============================================================

proteinas$muy_corta <- proteinas$length_aa < LONGITUD_MINIMA


# ============================================================
# 7. TABLA DE CONTROL DE CALIDAD (para documentar qué se
#    filtró y por qué, ANTES de tocar los datos)
# ============================================================

qc <- data.frame(
  id                  = proteinas$id,
  length_aa           = proteinas$length_aa,
  tiene_ambiguos      = proteinas$tiene_ambiguos,
  caracteres_ambiguos = proteinas$caracteres_ambiguos,
  muy_corta           = proteinas$muy_corta
)
qc$se_excluye_de_PCA <- qc$tiene_ambiguos | qc$muy_corta

write.csv(
  qc,
  file      = paste0(carpeta_salida, "control_calidad_proteinas.csv"),
  row.names = FALSE
)

n_ambiguos <- sum(proteinas$tiene_ambiguos)
n_cortas   <- sum(proteinas$muy_corta)
n_excluidas <- sum(qc$se_excluye_de_PCA)
n_finales   <- nrow(proteinas) - n_excluidas

resumen_filtrado <- c(
  "==================================================================",
  "CRITERIOS DE FILTRADO APLICADOS ANTES DE LA PCA",
  "==================================================================",
  paste0("Total de proteínas leídas del FASTA:            ", nrow(proteinas)),
  "",
  paste0("1) Caracteres ambiguos/no estándar (", paste(CARACTERES_AMBIGUOS, collapse = ", "), ", u otro"),
  paste0("   símbolo fuera de los 20 aminoácidos estándar):"),
  paste0("   Proteínas con al menos un carácter ambiguo:   ", n_ambiguos),
  "",
  paste0("2) Longitud inusualmente corta (umbral: < ", LONGITUD_MINIMA, " aa):"),
  paste0("   Proteínas por debajo del umbral:              ", n_cortas),
  "",
  paste0("Total de proteínas excluidas de la PCA:          ", n_excluidas,
         " (una proteína puede cumplir ambos criterios a la vez)"),
  paste0("Proteínas usadas en la PCA:                      ", n_finales),
  "==================================================================",
  "El detalle completo, proteína por proteína, está en:",
  "  control_calidad_proteinas.csv",
  "=================================================================="
)

writeLines(resumen_filtrado, paste0(carpeta_salida, "criterios_de_filtrado.txt"))
cat(paste(resumen_filtrado, collapse = "\n"), "\n")


# ============================================================
# 8. FILTRAR ANTES DE LA PCA
#    (se excluyen las proteínas con caracteres ambiguos y/o
#    con longitud menor al umbral definido arriba)
# ============================================================

proteinas_pca <- proteinas[!proteinas$tiene_ambiguos & !proteinas$muy_corta, ]

if (nrow(proteinas_pca) < 3) {
  stop("Quedan muy pocas proteínas después del filtrado como para correr una PCA. Revisá 'criterios_de_filtrado.txt'.")
}


# ============================================================
# 9. TABLA FINAL: UNA PROTEÍNA POR FILA, 20 FRECUENCIAS COMO
#    VARIABLES (esto es lo que entra a la PCA)
# ============================================================

tabla_pca <- proteinas_pca[, c("id", AA_20)]

write.csv(
  tabla_pca,
  file      = paste0(carpeta_salida, "frecuencias_aminoacidos_para_PCA.csv"),
  row.names = FALSE
)


# ============================================================
# 10. PCA
#     Se estandarizan las variables (center + scale) porque
#     algunos aminoácidos son mucho más frecuentes que otros en
#     promedio (p.ej. Leucina vs. Triptófano), y sin
#     estandarizar la PCA quedaría dominada por esa diferencia
#     de escala en vez de por los patrones de composición.
# ============================================================

matriz_pca <- as.matrix(tabla_pca[, AA_20])
rownames(matriz_pca) <- tabla_pca$id

pca <- prcomp(matriz_pca, center = TRUE, scale. = TRUE)

# Varianza explicada por cada componente
varianza <- (pca$sdev^2) / sum(pca$sdev^2)
var_PC1  <- round(varianza[1] * 100, 1)
var_PC2  <- round(varianza[2] * 100, 1)

# Guardar el resumen completo de la PCA (importancia de
# cada componente) y los loadings (contribución de cada
# aminoácido a cada componente)
sink(paste0(carpeta_salida, "PCA_resumen.txt"))
cat("Resumen de la PCA (importancia de cada componente):\n\n")
print(summary(pca))
sink()

write.csv(
  as.data.frame(pca$rotation),
  file      = paste0(carpeta_salida, "PCA_loadings_aminoacidos.csv"),
  row.names = TRUE
)

# Tabla de scores (coordenadas de cada proteína en la PCA),
# con la longitud de la proteína como columna extra, útil para
# colorear la figura
scores <- as.data.frame(pca$x)
scores$id        <- tabla_pca$id
scores$length_aa <- proteinas_pca$length_aa

write.csv(
  scores,
  file      = paste0(carpeta_salida, "PCA_scores_proteinas.csv"),
  row.names = FALSE
)


# ============================================================
# 11. FIGURA DE PCA, CALIDAD DE PUBLICACIÓN
#     Mismo lenguaje visual que el resto de las figuras del
#     proyecto: Helvetica, tema minimalista, gradiente
#     azul/rojo (acá coloreando por longitud de la proteína).
# ============================================================

p_pca <- ggplot(
  scores,
  aes(x = PC1, y = PC2, color = length_aa)
) +
  
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.4) +
  
  geom_point(size = 2.2, alpha = 0.8) +
  
  scale_color_gradient2(
    low      = color_bajo,
    mid      = "grey80",
    high     = color_alto,
    midpoint = median(scores$length_aa),
    name     = "Protein\nlength (aa)"
  ) +
  
  labs(
    title    = "PCA of amino acid composition",
    subtitle = paste0("Relative frequencies of the 20 standard amino acids  |  n = ",
                      nrow(scores), " proteins"),
    x = paste0("PC1 (", var_PC1, "% of variance)"),
    y = paste0("PC2 (", var_PC2, "% of variance)")
  ) +
  
  theme_minimal(base_size = 11, base_family = font_family) +
  
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 13,
                                 margin = margin(b = 4)),
    plot.subtitle = element_text(hjust = 0.5, size = 9.5, color = "grey40",
                                 margin = margin(b = 14)),
    axis.title    = element_text(color = "black", size = 10),
    axis.text     = element_text(color = "black", size = 9),
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_line(color = "grey92", linewidth = 0.35),
    legend.title       = element_text(size = 9),
    legend.text        = element_text(size = 8),
    legend.position    = "right",
    plot.margin        = margin(16, 20, 12, 12)
  )

print(p_pca)


# ============================================================
# 12. GUARDAR LA FIGURA (PDF vectorial + TIFF 600 dpi)
# ============================================================

ggsave(
  filename = paste0(carpeta_salida, "PCA_composicion_aminoacidos.pdf"),
  plot     = p_pca,
  width    = 7.5,
  height   = 6,
  units    = "in",
  device   = cairo_pdf
)

ggsave(
  filename    = paste0(carpeta_salida, "PCA_composicion_aminoacidos.tiff"),
  plot        = p_pca,
  width       = 7.5,
  height      = 6,
  units       = "in",
  dpi         = 600,
  compression = "lzw"
)

cat("\nListo. Todos los archivos quedaron en:\n", carpeta_salida, "\n")