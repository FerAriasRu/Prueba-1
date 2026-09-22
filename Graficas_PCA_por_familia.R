setwd("/home/luisafernanda/Documentos/Tarea1/Genoma_anotado")

```r
# ============================================================
# PROTEÍNAS Y FAMILIAS Pfam DESDE UN GFF3
# ============================================================

# Archivo GFF3
gff_file <- "GCA_033216535.1_TgRH_pasteur.gff3"


# ============================================================
# 1. LEER EL GFF3
# ============================================================

gff <- readLines(gff_file)

# Eliminar comentarios y líneas vacías
gff <- gff[
  !grepl("^#", gff) &
    nchar(trimws(gff)) > 0
]

# Separar las 9 columnas del GFF3
gff_split <- strsplit(gff, "\t", fixed = TRUE)

# Convertir a data.frame
gff_df <- do.call(
  rbind,
  lapply(gff_split, function(x) {
    length(x) <- 9
    x
  })
)

gff_df <- as.data.frame(gff_df, stringsAsFactors = FALSE)

colnames(gff_df) <- c(
  "seqid",
  "source",
  "type",
  "start",
  "end",
  "score",
  "strand",
  "phase",
  "attributes"
)


# ============================================================
# 2. SELECCIONAR LAS ANOTACIONES Pfam
# ============================================================

pfam <- gff_df[
  gff_df$seqid %in% c(
    "CM065576.1",
    "CM065577.1",
    "CM065578.1",
    "CM065579.1",
    "CM065580.1",
    "CM065581.1",
    "CM065582.1",
    "CM065583.1",
    "CM065584.1",
    "CM065585.1",
    "CM065586.1",
    "CM065587.1",
    "CM065588.1"
  ) &
    gff_df$source == "Pfam" &
    gff_df$type == "protein_match",
]

cat("Número de anotaciones Pfam:", nrow(pfam), "\n")


# ============================================================
# 3. EXTRAER EL IDENTIFICADOR DE LA PROTEÍNA
# ============================================================

# Ejemplo:
# Parent=TgRH_000005700.1:pep

pfam$protein <- sub(
  ".*Parent=([^;:]+).*",
  "\\1",
  pfam$attributes
)


# ============================================================
# 4. EXTRAER EL IDENTIFICADOR Pfam
# ============================================================

# Ejemplo:
# Name=PF00240

pfam$Pfam <- sub(
  ".*Name=([^;]+).*",
  "\\1",
  pfam$attributes
)


# ============================================================
# 5. EXTRAER LA DESCRIPCIÓN DE LA FAMILIA
# ============================================================

# Ejemplo:
# signature_desc=Ubiquitin family

pfam$family <- sub(
  ".*signature_desc=([^;]+).*",
  "\\1",
  pfam$attributes
)


# ============================================================
# 6. CREAR LA TABLA FINAL
# ============================================================

protein_family <- pfam[
  ,
  c("protein", "Pfam", "family")
]

# Eliminar posibles duplicados
protein_family <- unique(protein_family)

# Ordenar
protein_family <- protein_family[
  order(protein_family$protein, protein_family$Pfam),
]


# ============================================================
# 7. VER LOS PRIMEROS RESULTADOS
# ============================================================

head(protein_family, 20)


# ============================================================
# 8. GUARDAR LA TABLA
# ============================================================

write.table(
  protein_family,
  file = "proteins_Pfam_families.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nTabla guardada como: proteins_Pfam_families.tsv\n"
)



# ============================================================
# PCA DE COMPOSICIÓN DE AMINOÁCIDOS
# PASO 1: FRECUENCIAS DE LOS 20 AMINOÁCIDOS
# + FAMILIA Pfam
# ============================================================

library(Biostrings)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)


# ============================================================
# 1. ARCHIVOS
# ============================================================

fasta_file <- "proteins_chromosomes.fasta"

pfam_file <- "proteins_Pfam_families.tsv"


# ============================================================
# 2. LEER LAS PROTEÍNAS
# ============================================================

proteins <- readAAStringSet(fasta_file)

protein_ids <- names(proteins)

# En caso de que los nombres tengan información adicional
protein_ids <- sub(" .*", "", protein_ids)

names(proteins) <- protein_ids


# ============================================================
# 3. AMINOÁCIDOS ESTÁNDAR
# ============================================================

aa <- c(
  "A", "C", "D", "E", "F",
  "G", "H", "I", "K", "L",
  "M", "N", "P", "Q", "R",
  "S", "T", "V", "W", "Y"
)


# ============================================================
# 4. CALCULAR FRECUENCIAS DE AMINOÁCIDOS
# ============================================================

aa_frequency <- function(sequence) {
  
  sequence <- toupper(as.character(sequence))
  
  # contar aminoácidos
  counts <- table(
    factor(
      strsplit(sequence, "")[[1]],
      levels = aa
    )
  )
  
  # frecuencia relativa
  frequencies <- counts / sum(counts)
  
  return(as.numeric(frequencies))
}


# Aplicar a todas las proteínas

freq_matrix <- t(
  sapply(
    proteins,
    aa_frequency
  )
)

colnames(freq_matrix) <- aa

freq_df <- as.data.frame(freq_matrix)

freq_df$protein <- names(proteins)


# ============================================================
# 5. LONGITUD DE LAS PROTEÍNAS
# ============================================================

freq_df$length_aa <- width(proteins)


# ============================================================
# 6. REORDENAR COLUMNAS
# ============================================================

freq_df <- freq_df %>%
  select(
    protein,
    length_aa,
    all_of(aa)
  )


# ============================================================
# 7. LEER LAS FAMILIAS Pfam
# ============================================================

pfam <- read.delim(
  pfam_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# 8. UNA FAMILIA POR PROTEÍNA
#
# Algunas proteínas tienen varios Pfam.
# Para este primer análisis conservamos la primera
# familia encontrada.
# ============================================================

protein_family <- pfam %>%
  filter(
    !is.na(protein),
    !is.na(Pfam),
    !is.na(family)
  ) %>%
  distinct(protein, Pfam, family) %>%
  group_by(protein) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    protein,
    Pfam,
    family
  )


# ============================================================
# 9. UNIR FRECUENCIAS + FAMILIA
# ============================================================

pca_data <- freq_df %>%
  left_join(
    protein_family,
    by = "protein"
  )


# ============================================================
# 10. PROTEÍNAS SIN FAMILIA Pfam
# ============================================================

pca_data$family[
  is.na(pca_data$family)
] <- "No Pfam annotation"


# ============================================================
# 11. REVISAR RESULTADO
# ============================================================

head(pca_data)

cat(
  "\nNúmero total de proteínas:",
  nrow(pca_data),
  "\n"
)

cat(
  "Proteínas con familia Pfam:",
  sum(pca_data$family != "No Pfam annotation"),
  "\n"
)

cat(
  "Proteínas sin familia Pfam:",
  sum(pca_data$family == "No Pfam annotation"),
  "\n"
)


# ============================================================
# 12. TABLA DE NÚMERO DE PROTEÍNAS POR FAMILIA
# ============================================================

family_summary <- pca_data %>%
  count(
    family,
    sort = TRUE
  )

print(family_summary)


# ============================================================
# 13. GUARDAR MATRIZ
# ============================================================

write.table(
  pca_data,
  file = "protein_aminoacid_frequencies_Pfam.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nArchivo generado:\n",
  "protein_aminoacid_frequencies_Pfam.tsv\n"
)


# ============================================================
# PCA DE COMPOSICIÓN DE AMINOÁCIDOS POR FAMILIA Pfam
# ============================================================

library(ggplot2)
library(dplyr)

# ------------------------------------------------------------
# 1. Variables para el PCA
# ------------------------------------------------------------

aa <- c(
  "A", "C", "D", "E", "F",
  "G", "H", "I", "K", "L",
  "M", "N", "P", "Q", "R",
  "S", "T", "V", "W", "Y"
)

X <- pca_data %>%
  select(all_of(aa))

# ------------------------------------------------------------
# 2. PCA
# ------------------------------------------------------------

pca <- prcomp(
  X,
  center = TRUE,
  scale. = TRUE
)

# Porcentaje de varianza
variance <- pca$sdev^2
variance_percent <- variance / sum(variance) * 100

var_PC1 <- round(variance_percent[1], 1)
var_PC2 <- round(variance_percent[2], 1)

# ------------------------------------------------------------
# 3. Scores del PCA
# ------------------------------------------------------------

scores <- as.data.frame(pca$x)

scores$protein <- pca_data$protein
scores$family <- pca_data$family
scores$length_aa <- pca_data$length_aa

# ============================================================
# 4. ELIMINAR PROTEÍNAS SIN ANOTACIÓN Pfam
# ============================================================

scores_pfam <- scores %>%
  filter(
    !is.na(family),
    family != "No Pfam annotation"
  )


# ============================================================
# 5. SELECCIONAR LAS 15 FAMILIAS MÁS ABUNDANTES
# ============================================================

top_n <- 15

top_families <- scores_pfam %>%
  count(family, sort = TRUE) %>%
  slice_head(n = top_n) %>%
  pull(family)


# ============================================================
# 6. CONSERVAR SOLAMENTE ESAS FAMILIAS
# ============================================================

scores_pfam <- scores_pfam %>%
  filter(family %in% top_families)

scores_pfam$family_plot <- factor(
  scores_pfam$family,
  levels = top_families
)


# ============================================================
# 7. PALETA DE COLORES
# ============================================================

n_colors <- length(top_families)

family_colors <- scales::hue_pal(
  h = c(15, 375),
  c = 90,
  l = 65
)(n_colors)

names(family_colors) <- top_families


# ============================================================
# 8. PCA
# ============================================================

p_pca <- ggplot(
  scores_pfam,
  aes(
    x = PC1,
    y = PC2,
    color = family_plot
  )
) +
  
  geom_hline(
    yintercept = 0,
    color = "grey85",
    linewidth = 0.4
  ) +
  
  geom_vline(
    xintercept = 0,
    color = "grey85",
    linewidth = 0.4
  ) +
  
  geom_point(
    size = 2.2,
    alpha = 0.8
  ) +
  
  scale_color_manual(
    values = family_colors,
    name = "Protein family"
  ) +
  
  labs(
    title = "PCA of amino acid composition",
    
    subtitle = paste0(
      "Relative frequencies of the 20 standard amino acids | ",
      "Top 15 Pfam families | n = ",
      nrow(scores_pfam),
      " proteins"
    ),
    
    x = paste0(
      "PC1 (",
      var_PC1,
      "% of variance)"
    ),
    
    y = paste0(
      "PC2 (",
      var_PC2,
      "% of variance)"
    )
  ) +
  
  theme_minimal(
    base_size = 11,
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
      size = 9.5,
      color = "grey40",
      margin = margin(b = 14)
    ),
    
    axis.title = element_text(
      color = "black",
      size = 10
    ),
    
    axis.text = element_text(
      color = "black",
      size = 9
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.grid.major = element_line(
      color = "grey92",
      linewidth = 0.35
    ),
    
    legend.title = element_text(
      size = 9
    ),
    
    legend.text = element_text(
      size = 7.5
    ),
    
    legend.key.height = unit(
      0.45,
      "cm"
    ),
    
    legend.position = "right",
    
    plot.margin = margin(
      16, 20, 12, 12
    )
  )

print(p_pca)


# ============================================================
# 9. GUARDAR
# ============================================================

ggsave(
  "PCA_aminoacid_composition_top15_Pfam.png",
  plot = p_pca,
  width = 9,
  height = 6,
  units = "in",
  dpi = 600
)

ggsave(
  "PCA_aminoacid_composition_top15_Pfam.pdf",
  plot = p_pca,
  width = 9,
  height = 6,
  units = "in"
)


# ============================================================
# RESUMEN DE PROTEÍNAS Y FAMILIAS Pfam
# ============================================================

# Número total de proteínas
cat(
  "Número total de proteínas:",
  nrow(scores),
  "\n"
)

# Número de proteínas con Pfam
cat(
  "Proteínas con Pfam:",
  sum(
    !is.na(scores$family) &
      scores$family != "No Pfam annotation"
  ),
  "\n"
)

# Número de proteínas sin Pfam
cat(
  "Proteínas sin Pfam:",
  sum(
    is.na(scores$family) |
      scores$family == "No Pfam annotation"
  ),
  "\n\n"
)


# Tabla de proteínas por familia
family_summary <- scores %>%
  mutate(
    family = ifelse(
      is.na(family) |
        family == "No Pfam annotation",
      "No Pfam annotation",
      family
    )
  ) %>%
  count(
    family,
    name = "n_proteins"
  ) %>%
  arrange(desc(n_proteins))

print(family_summary)


