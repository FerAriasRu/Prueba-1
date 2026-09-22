#install.packages("BiocManager")
#BiocManager::install("Biostrings")
library(Biostrings)
packageVersion("Biostrings")
library(tidyverse)
#install.packages("dplyr")
library(dplyr)
library(stringr)
library(ggplot2)
# ============================================================
# 1. ARCHIVOS
# ============================================================

setwd("/home/luisafernanda/Documentos/Tarea1/Proteinas")

fasta_file <- "proteins.fasta"

gtf_file <- "/home/luisafernanda/Documentos/Tarea1/GCA_033216535.1/genes/GCA_033216535.1_TgRH_pasteur.augustus.gtf"

# Leer FASTA
proteinas <- readAAStringSet(fasta_file)

# IDs de las proteínas
protein_id <- names(proteinas)

# Eliminar descripción después de un espacio, si existe
protein_id <- sub(" .*", "", protein_id)

# Secuencias
secuencias <- as.character(proteinas)

# Longitud de cada proteína
longitud <- width(proteinas)

# ============================================================
# 2. LEER GTF
# ============================================================

gtf <- read.delim(
  gtf_file,
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "#",
  fill = TRUE,
  stringsAsFactors = FALSE
)

# Nombres de las columnas
colnames(gtf) <- c(
  "chromosome",
  "source",
  "feature",
  "start",
  "end",
  "score",
  "strand",
  "frame",
  "attributes"
)

head(gtf)

# ============================================================
# 3. EXTRAER IDENTIFICADORES
# ============================================================

gtf <- gtf %>%
  mutate(
    transcript_id = str_extract(
      attributes,
      '(?<=transcript_id ")[^"]+'
    ),
    
    gene_id = str_extract(
      attributes,
      '(?<=gene_id ")[^"]+'
    )
  )

# Revisar
head(
  gtf %>%
    select(
      chromosome,
      feature,
      transcript_id,
      gene_id
    )
)


# ============================================================
# 4. ANOTACIÓN GENÓMICA POR TRANSCRIPT
# ============================================================

gtf_cds <- gtf %>%
  filter(
    feature == "CDS",
    !is.na(transcript_id)
  )

anotacion <- gtf_cds %>%
  group_by(transcript_id) %>%
  summarise(
    chromosome = paste(unique(chromosome), collapse = ";"),
    gene_id = first(na.omit(gene_id)),
    start = min(start, na.rm = TRUE),
    end = max(end, na.rm = TRUE),
    strand = first(strand),
    
    # Longitud total de las regiones CDS
    CDS_length = sum(end - start + 1),
    
    .groups = "drop"
  )

head(anotacion)



# ============================================================
# 5. FRECUENCIAS DE AMINOÁCIDOS
# ============================================================

aa20 <- c(
  "A", "C", "D", "E", "F",
  "G", "H", "I", "K", "L",
  "M", "N", "P", "Q", "R",
  "S", "T", "V", "W", "Y"
)

frecuencias <- matrix(
  0,
  nrow = length(secuencias),
  ncol = length(aa20)
)

colnames(frecuencias) <- aa20

for (i in seq_along(secuencias)) {
  
  caracteres <- strsplit(secuencias[i], "")[[1]]
  
  conteos <- table(
    factor(
      caracteres,
      levels = aa20
    )
  )
  
  frecuencias[i, ] <- as.numeric(conteos) / longitud[i]
  
}

frecuencias <- as.data.frame(frecuencias)

frecuencias <- cbind(
  protein_id = protein_id,
  length = longitud,
  frecuencias
)


# ============================================================
# 6. UNIR INFORMACIÓN
# ============================================================

tabla_maestra <- frecuencias %>%
  left_join(
    anotacion,
    by = c("protein_id" = "transcript_id")
  )

# ============================================================
# 7. CONTROL DE CORRESPONDENCIA
# ============================================================

cat(
  "Total de proteínas:",
  nrow(tabla_maestra),
  "\n"
)

cat(
  "Proteínas con cromosoma:",
  sum(!is.na(tabla_maestra$chromosome)),
  "\n"
)

cat(
  "Proteínas sin cromosoma:",
  sum(is.na(tabla_maestra$chromosome)),
  "\n"
)

# ============================================================
# 8. GUARDAR
# ============================================================

write.table(
  tabla_maestra,
  "tabla_maestra_proteinas.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# 9. PCA
# ============================================================

# Seleccionar únicamente proteínas con anotación cromosómica
# y sin caracteres ambiguos, si ya has creado ese filtro.

datos_pca <- tabla_maestra %>%
  filter(!is.na(chromosome))

matriz_pca <- datos_pca %>%
  select(all_of(aa20))

pca <- prcomp(
  matriz_pca,
  center = TRUE,
  scale. = TRUE
)

# ============================================================
# 10. SCORES + CROMOSOMA
# ============================================================

scores <- as.data.frame(pca$x)

scores$protein_id <- datos_pca$protein_id

scores$chromosome <- datos_pca$chromosome

scores$length <- datos_pca$length

datos_cromosomas <- tabla_maestra %>%
  filter(str_detect(chromosome, "^CM"))

scores_CM <- scores %>%
  filter(str_detect(chromosome, "^CM"))
# ============================================================
# 11. GRÁFICA PCA POR CROMOSOMA_ feA
# ============================================================

ggplot(
  scores_CM,
  aes(
    x = PC1,
    y = PC2,
    color = chromosome
  )
) +
  
  geom_point(
    size = 1.5,
    alpha = 0.6
  ) +
  
  theme_classic() +
  
  labs(
    title = "PCA de composición de aminoácidos por cromosoma",
    x = "PC1",
    y = "PC2",
    color = "Cromosoma"
  )

# ============================================================
# 12. GRÁFICA PCA POR CROMOSOMA_Publicacion
# ============================================================
# ============================================================
# 1. PORCENTAJE DE VARIANZA EXPLICADA
# ============================================================
var_PC1 <- round(
  100 * summary(pca)$importance[2, "PC1"],
  1
)

var_PC2 <- round(
  100 * summary(pca)$importance[2, "PC2"],
  1
)


# ============================================================
# 2. FIGURA PCA
# ============================================================

p_pca <- ggplot(
  scores_CM,
  aes(
    x = PC1,
    y = PC2,
    color = chromosome
  )
) +
  
  # Líneas de referencia en los ejes
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
  
  # Puntos
  geom_point(
    size = 2.2,
    alpha = 0.8
  ) +
  
  # Escala de colores por cromosoma
  scale_color_viridis_d(
    name = "Chromosome",
    option = "turbo"
  ) +
  
  # Títulos y ejes
  labs(
    title = "PCA of amino acid composition",
    
    subtitle = paste0(
      "Chromosomes CM | Relative frequencies of the 20 ",
      "standard amino acids | n = ",
      nrow(scores_CM),
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
  
  # Tema minimalista
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
      size = 8
    ),
    
    legend.position = "right",
    
    plot.margin = margin(
      16, 20, 12, 12
    )
  )

# Mostrar figura
print(p_pca)


colnames(tabla_maestra)


