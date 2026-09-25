# ============================================================
# ANÁLISIS DE CORRESPONDENCIA (COA) DE USO DE CODONES
# Coloreado por familias funcionales curadas (mismas que script 1)
# Metodología de referencia: Charif et al. 2005 (seqinr / ade4)
#   http://pbil.univ-lyon1.fr/datasets/charif04/
#
# CORRECCIÓN RESPECTO A TU SCRIPT ANTERIOR:
#   El script viejo buscaba los patrones de familia (SRS, ROP, GRA...)
#   en la columna `attributes` de las líneas mRNA del GFF3. Esa columna
#   normalmente solo trae ID=/Parent=/Name=, NO la descripción funcional
#   (`product=`), así que casi todo caía en "Otras / No clasificadas".
#   Aquí se reutiliza directamente `products.tsv` (la misma fuente que
#   usaste, correctamente, en el script de PCA fisicoquímico) y se le
#   aplica el MISMO clasificador curado, para que las familias sean
#   consistentes entre ambos análisis y el join por ID funcione.
#
# Input esperado:
#   - cds_transcritos_chromosomes.fasta  -> generado con gffread -x
#   - products.tsv                       -> mismo archivo que en el
#                                            script de PCA fisicoquímico
#                                            (protein_id <TAB> product)
# ============================================================

# ---- 0. Paquetes ----
# install.packages(c("seqinr","dplyr","stringr","readr","tibble","ade4","ggplot2","ggrepel","RColorBrewer"))
library(seqinr)
library(dplyr)
library(stringr)
library(readr)
library(tibble)
library(ade4)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

set.seed(123)
setwd("/home/larias/Descargas/Transcritos")




# ============================================================
# FUNCIÓN: Normalizar IDs (idéntica a la del script de PCA fisicoquímico)
# ============================================================
normalizar_id <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- sub("\\..*", "", x)   # quita ".1:pep" o ".1"
  x <- sub(":.*", "", x)     # quita ":pep" si quedó
  x <- str_trim(x)
  return(x)
}

# ============================================================
# 1. FAMILIAS: leer products.tsv y aplicar el MISMO clasificador
#    curado usado en el script de PCA fisicoquímico
# ============================================================
cat("=== 1. LEYENDO products.tsv Y ASIGNANDO FAMILIAS ===\n")

products <- read.delim("products.tsv", header = FALSE,
                       col.names = c("protein_id", "product"),
                       stringsAsFactors = FALSE) %>%
  mutate(
    protein_id = normalizar_id(protein_id),
    product = URLdecode(product),
    product = str_remove(product, "^term="),
    product = str_split_i(product, ";", 1),
    product = str_remove(product, ", putative$"),
    product = str_remove(product, " putative$"),
    product = str_trim(product)
  )

familias <- products %>%
  mutate(familia = case_when(
    # ===== APICOMPLEXA-ESPECÍFICAS =====
    str_detect(product, regex("SAG-related|\\bSRS[0-9]*\\b|surface antigen", ignore_case = TRUE)) ~ "SAG/SRS (superficie)",
    str_detect(product, regex("Toxoplasma gondii family", ignore_case = TRUE)) ~ "Familias Tg (A-E)",
    str_detect(product, regex("dense granule|\\bGRA[0-9]*\\b", ignore_case = TRUE)) ~ "GRA (gránulos densos)",
    str_detect(product, regex("rhoptry neck|\\bRON[0-9]*\\b", ignore_case = TRUE)) ~ "RON (cuello de roptria)",
    str_detect(product, regex("rhoptry|\\bROP[0-9]*\\b", ignore_case = TRUE)) ~ "ROP (roptrias)",
    str_detect(product, regex("microneme|\\bMIC[0-9]*\\b", ignore_case = TRUE)) ~ "MIC (micronemas)",
    str_detect(product, regex("inner membrane complex|\\bIMC[0-9]*\\b", ignore_case = TRUE)) ~ "IMC (complejo membrana)",
    str_detect(product, regex("apical cap|apical cone|apical neck", ignore_case = TRUE)) ~ "Apical (cap/cone)",
    str_detect(product, regex("alveolin", ignore_case = TRUE)) ~ "Alveolinas",
    str_detect(product, regex("\\bKRUF[0-9]*\\b", ignore_case = TRUE)) ~ "KRUF",
    
    # ===== Tráfico / señalización de membrana =====
    str_detect(product, regex("SNARE|vesicle transport|clathrin|coatomer|dynamin|adaptin|\\bVPS[0-9]+\\b", ignore_case = TRUE)) ~ "Tráfico vesicular / SNARE",
    str_detect(product, regex("\\bRab\\b|Rho GTPase|Ras-related|ADP-ribosylation factor|\\bseptin", ignore_case = TRUE)) ~ "GTPasas pequeñas",
    str_detect(product, regex("apicoplast|\\bplastid\\b|ferredoxin|iron-sulfur|Fe-S cluster", ignore_case = TRUE)) ~ "Apicoplasto / Fe-S",
    
    # ===== REGULADORES =====
    str_detect(product, regex("AP2 domain transcription", ignore_case = TRUE)) ~ "AP2 factores transcripción",
    str_detect(product, regex("\\bcyclin\\b|cell cycle|\\bcentrin\\b|aurora kinase", ignore_case = TRUE)) ~ "Ciclo celular",
    str_detect(product, regex("\\bkinase\\b|kinase-like|\\bNEK[0-9]*\\b", ignore_case = TRUE)) ~ "Quinasas",
    str_detect(product, regex("phosphatase", ignore_case = TRUE)) ~ "Fosfatasas",
    str_detect(product, regex("subtilisin|protease|peptidase|proteasome", ignore_case = TRUE)) ~ "Proteasas",
    
    # ===== Replicación / reparación de DNA =====
    str_detect(product, regex("DNA polymerase|primase|\\bMCM[0-9]+\\b|mismatch repair|replication factor|replication protein", ignore_case = TRUE)) ~ "Replicación/reparación DNA",
    
    # ===== METABOLISMO =====
    str_detect(product, regex("ribosomal|ribosome", ignore_case = TRUE)) ~ "Ribosomales",
    str_detect(product, regex("elongation factor|initiation factor|translation|tRNA", ignore_case = TRUE)) ~ "Traducción/tRNA",
    str_detect(product, regex("RNA polymerase|mediator complex|transcription", ignore_case = TRUE)) ~ "RNA polimerasa / transcripción",
    str_detect(product, regex("cytochrome|ATP synthase|oxidase|reductase|dehydrogenase", ignore_case = TRUE)) ~ "Metabolismo energético",
    str_detect(product, regex("transporter|permease|carrier|channel|ABC transporter|facilitator", ignore_case = TRUE)) ~ "Transportadores",
    str_detect(product, regex("glycosyltransferase|GPI anchor|glycosylphosphatidylinositol", ignore_case = TRUE)) ~ "Glicosilación / GPI",
    str_detect(product, regex("synthase|transferase|isomerase|lyase|ligase|mutase|hydrolase|phosphodiesterase|palmitoyltransferase", ignore_case = TRUE)) ~ "Metabolismo general",
    
    # ===== OTROS =====
    str_detect(product, regex("histone|methyltransferase.*histone|acetyltransferase|chromatin|SWI2/SNF2", ignore_case = TRUE)) ~ "Modificadores de cromatina",
    str_detect(product, regex("ubiquitin|HECT|OTU domain|carboxyl-terminal hydrolase", ignore_case = TRUE)) ~ "Ubiquitina / UPS",
    str_detect(product, regex("DnaJ|heat shock|Hsp|chaperone|thioredoxin", ignore_case = TRUE)) ~ "Chaperonas / Heat shock",
    str_detect(product, regex("tubulin|kinesin|dynein|myosin|actin|centrosomal|basal complex", ignore_case = TRUE)) ~ "Citoesqueleto / Motores",
    str_detect(product, regex("helicase|DEAD/DEAH", ignore_case = TRUE)) ~ "Helicasas",
    str_detect(product, regex("EF hand|calcium-dependent|calmodulin", ignore_case = TRUE)) ~ "Señalización calcio",
    str_detect(product, regex("MORN|kelch|RCC1|TBC domain|RAP domain|endonuclease|exonuclease", ignore_case = TRUE)) ~ "Otras enzimáticas",
    
    # ===== REPETICIONES (función desconocida) =====
    str_detect(product, regex("WD domain|RNA recognition motif|tetratricopeptide|leucine rich|HEAT repeat|ankyrin|zinc finger", ignore_case = TRUE)) ~ "Repeticiones (función desconocida)",
    
    # ===== HIPOTÉTICAS =====
    str_detect(product, regex("hypothetical protein, conserved|conserved hypothetical", ignore_case = TRUE)) ~ "Hipotéticas conservadas",
    str_detect(product, regex("hypothetical", ignore_case = TRUE)) ~ "Hipotéticas",
    
    TRUE ~ "Otras"
  )) %>%
  select(protein_id, familia)

cat("Transcritos con producto/familia:", nrow(familias), "\n")
cat("Distribución de familias:\n")
print(familias %>% count(familia, sort = TRUE))
cat("\n")

# ============================================================
# 2. LEER SECUENCIAS CDS
# ============================================================
cat("=== 2. LEYENDO CDS ===\n")
cds <- read.fasta("cds_transcritos_chromosomes.fasta",
                  seqtype = "DNA",
                  as.string = FALSE,
                  set.attributes = FALSE)

names(cds) <- normalizar_id(names(cds))
cat("Total de secuencias leídas:", length(cds), "\n")

# Diagnóstico de coincidencia ANTES de seguir (para no repetir el bug)
n_match <- sum(names(cds) %in% familias$protein_id)
cat("IDs del fasta con familia asignada:", n_match, "de", length(cds), "\n")
if (n_match == 0) {
  cat("⚠️  Ejemplos names(cds):\n"); print(head(names(cds), 5))
  cat("⚠️  Ejemplos familias$protein_id:\n"); print(head(familias$protein_id, 5))
  stop("Los IDs no coinciden. Revisa normalizar_id() vs el formato real de tus headers.")
}
cat("\n")

# ============================================================
# 3. FILTRADO DE CALIDAD
# ============================================================
cat("=== 3. FILTRADO DE CALIDAD ===\n")
long_ok   <- sapply(cds, length) %% 3 == 0
n_codones <- sapply(cds, length) / 3

stops <- c("taa", "tag", "tga")
stop_interno <- sapply(cds, function(s) {
  cod <- seqinr::splitseq(s, frame = 0, word = 3)
  if (length(cod) < 2) return(TRUE)
  any(cod[-length(cod)] %in% stops)
})

filtro <- long_ok & !stop_interno & n_codones >= 100
cds_f  <- cds[filtro]
cat("Secuencias tras filtrado:", length(cds_f), "de", length(cds), "\n\n")

# ============================================================
# 4. CONTEOS DE CODONES (para el COA se usan conteos crudos,
#    no RSCU -> así es el método original de Charif et al. 2005)
# ============================================================
cat("=== 4. CALCULANDO CONTEOS DE CODONES ===\n")
conteos_list <- lapply(cds_f, uco, index = "eff", as.data.frame = FALSE)
conteos_mat  <- do.call(rbind, conteos_list)
rownames(conteos_mat) <- names(cds_f)

# Se excluyen los 3 codones de stop y los 2 codones únicos por aa
# (Met=atg, Trp=tgg), que no aportan señal de sesgo sinónimo
codones_excluir <- c("taa", "tag", "tga", "atg", "tgg")
conteos_mat <- conteos_mat[, !(colnames(conteos_mat) %in% codones_excluir)]

# ============================================================
# 5. UNIR CONTEOS CON FAMILIAS
# ============================================================
cat("=== 5. UNIENDO CON FAMILIAS ===\n")
df_ids <- data.frame(protein_id = rownames(conteos_mat), stringsAsFactors = FALSE) %>%
  inner_join(familias, by = "protein_id")

# Reordenar conteos_mat según df_ids y quitar los que no tienen familia
conteos_mat <- conteos_mat[df_ids$protein_id, ]
conteos_mat <- conteos_mat[, colSums(conteos_mat) > 0]  # quita codones sin uso

cat("Genes con familia y conteo de codones disponibles:", nrow(conteos_mat), "\n\n")
if (nrow(conteos_mat) == 0) stop("El join sigue dando 0. Revisa el diagnóstico del paso 2.")

# ============================================================
# 6. ANÁLISIS DE CORRESPONDENCIA (COA) - ade4::dudi.coa
# ============================================================
cat("=== 6. CALCULANDO COA ===\n")
coa <- dudi.coa(as.data.frame(conteos_mat), scannf = FALSE, nf = 2)
var_exp_coa <- (coa$eig / sum(coa$eig)) * 100
cat("Inercia explicada Eje1:", round(var_exp_coa[1], 1), "%  |  Eje2:",
    round(var_exp_coa[2], 1), "%\n\n")

familia_fac <- factor(df_ids$familia)
n_fam <- nlevels(familia_fac)
paleta_fam <- colorRampPalette(brewer.pal(12, "Paired"))(n_fam)

# ============================================================
# 7. GRÁFICA ESTILO CHARIF ET AL. 2005 (ade4 clásico)
#    Puntos = genes coloreados por familia (con elipses)
#    Cajas  = codones (coa$co), en su posición real
#    -> Esta es la réplica exacta del estilo de la imagen de referencia
# ============================================================
cat("=== 7. GRÁFICA ESTILO ade4 (genes + codones superpuestos) ===\n")

png("COA_codones_por_familia_ade4.png", width = 3000, height = 3000, res = 300)
s.class(coa$li, fac = familia_fac, col = paleta_fam,
        cpoint = 1, pch = 16, clabel = 0.8, cellipse = 1,
        sub = "COA de uso de codones coloreado por familia funcional",
        possub = "topleft")
s.label(coa$co, add.plot = TRUE, boxes = TRUE, clabel = 0.7)
dev.off()

pdf("COA_codones_por_familia_ade4.pdf", width = 10, height = 10)
s.class(coa$li, fac = familia_fac, col = paleta_fam,
        cpoint = 1, pch = 16, clabel = 0.8, cellipse = 1,
        sub = "COA de uso de codones coloreado por familia funcional",
        possub = "topleft")
s.label(coa$co, add.plot = TRUE, boxes = TRUE, clabel = 0.7)
dev.off()

cat("Generado: COA_codones_por_familia_ade4.png / .pdf\n\n")

# ============================================================
# 8. VERSIÓN ggplot2 (más editable / calidad de publicación)
#    Con los codones MÁS REPRESENTATIVOS resaltados: se define
#    "representativo" como el codón con mayor distancia al origen
#    (mayor contribución conjunta a Eje1+Eje2), que es lo que en
#    ade4 se llama "contribución relativa" (coa$co + coa$cos2)
# ============================================================
cat("=== 8. GRÁFICA ggplot2 CON CODONES MÁS REPRESENTATIVOS ===\n")

scores_genes <- coa$li %>%
  rownames_to_column("protein_id") %>%
  left_join(df_ids, by = "protein_id")

scores_codones <- coa$co %>%
  rownames_to_column("codon") %>%
  rename(Axis1 = Comp1, Axis2 = Comp2) %>%
  mutate(dist_origen = sqrt(Axis1^2 + Axis2^2))

# Top N codones más representativos (mayor distancia al origen = mayor
# contribución al patrón de sesgo de codones capturado por los 2 ejes)
N_TOP <- 20
codones_top <- scores_codones %>% arrange(desc(dist_origen)) %>% slice(1:N_TOP)

cat("Codones más representativos (mayor contribución a Ejes 1-2):\n")
print(codones_top %>% select(codon, Axis1, Axis2, dist_origen))

p_coa <- ggplot() +
  geom_point(data = scores_genes, aes(Axis1, Axis2, color = familia),
             size = 1.8, alpha = 0.75, stroke = 0) +
  stat_ellipse(data = scores_genes, aes(Axis1, Axis2, color = familia, group = familia),
               type = "t", level = 0.95, linewidth = 0.4, show.legend = FALSE) +
  geom_point(data = codones_top, aes(Axis1, Axis2),
             color = "black", size = 1.6, shape = 15) +
  geom_label_repel(data = codones_top, aes(Axis1, Axis2, label = codon),
                   size = 3, fontface = "bold", max.overlaps = 30,
                   label.padding = unit(0.12, "lines"), segment.size = 0.3) +
  scale_color_manual(values = paleta_fam, name = "Familia funcional") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  labs(
    title = "COA de uso de codones por familia funcional",
    subtitle = paste0(nrow(scores_genes), " genes | ", n_fam, " familias | ",
                      "codones más representativos (top ", N_TOP, ") resaltados"),
    x = paste0("Eje 1 (", round(var_exp_coa[1], 1), "% inercia)"),
    y = paste0("Eje 2 (", round(var_exp_coa[2], 1), "% inercia)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    legend.text = element_text(size = 7),
    legend.title = element_text(face = "bold", size = 9),
    legend.key.size = unit(0.35, "cm")
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3, alpha = 1)))

ggsave("COA_codones_por_familia_ggplot.png", p_coa, width = 10, height = 7.5, dpi = 600)
ggsave("COA_codones_por_familia_ggplot.pdf", p_coa, width = 10, height = 7.5, dpi = 600)
print(p_coa)

cat("\n✓ Listo. Archivos generados:\n",
    " - familias.csv (si querés guardarla, ver línea comentada abajo)\n",
    " - COA_codones_por_familia_ade4.png / .pdf   (réplica estilo Charif/PBIL)\n",
    " - COA_codones_por_familia_ggplot.png / .pdf  (versión editable, codones top resaltados)\n")

# write_csv(familias, "familias.csv")

# ============================================================
# 9. FACET: primer panel = SOLO las letras de los codones (sin
#    color, sin nube de genes) a modo de referencia. Los paneles
#    siguientes = nube general (gris) + familia de ese panel
#    resaltada en color, SIN codones superpuestos.
#    Mismos límites de ejes en todos los paneles para que sean
#    comparables (igual criterio que el script de PCA fisicoquímico,
#    sección 11).
# ============================================================
cat("=== 9. FACET: panel de codones + paneles por familia ===\n")

# Límites de ejes comunes a TODOS los paneles
rango_x <- range(scores_genes$Axis1) * 1.05
rango_y <- range(scores_genes$Axis2) * 1.05

familias_niveles <- sort(unique(scores_genes$familia))
nombre_panel_codones <- "Codones (referencia)"
niveles_panel <- c(nombre_panel_codones, familias_niveles)

# --- Panel 1: SOLO codones, sin color, sin nube de genes ---
df_codones_facet <- codones_top %>%
  mutate(panel = factor(nombre_panel_codones, levels = niveles_panel))

# --- Paneles siguientes: nube general (gris) replicada en cada
#     panel de familia (pero NO en el panel de codones) ---
df_bg_facet <- scores_genes %>%
  select(Axis1, Axis2) %>%
  tidyr::crossing(panel = factor(familias_niveles, levels = niveles_panel))

# --- Paneles siguientes: puntos de la familia de cada panel,
#     coloreados (sin codones) ---
df_family_facet <- scores_genes %>%
  mutate(panel = factor(familia, levels = niveles_panel))

p_facet <- ggplot() +
  
  # Paneles de familia: nube general en gris de fondo
  geom_point(data = df_bg_facet, aes(Axis1, Axis2),
             color = "grey82", size = 0.6, alpha = 0.5, stroke = 0) +
  
  # Panel de codones: solo texto, negro, sin puntos ni nube
#  geom_text(data = df_codones_facet, aes(Axis1, Axis2, label = codon),
#            size = 1.8, fontface = "bold", color = "black") +
  
  ggrepel::geom_text_repel(
    data = df_codones_facet,
    aes(Axis1, Axis2, label = codon),
    size = 2.5,
    fontface = "bold",
    color = "black",
    box.padding = 0.4,
    point.padding = 0.3,
    force = 2,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey50",
    segment.size = 0.25
  ) +
  
  # Paneles de familia: puntos de esa familia, coloreados
  geom_point(data = df_family_facet, aes(Axis1, Axis2, color = familia),
             size = 1.0, alpha = 0.9, stroke = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  facet_wrap(~ panel, ncol = 5) +                   # scales = "fixed" -> mismos ejes en todos
  scale_color_manual(values = paleta_fam, guide = "none") +
  coord_equal(xlim = rango_x, ylim = rango_y) +
  labs(
    title = "COA de uso de codones: referencia de codones + nube por familia funcional",
    subtitle = paste0("Panel 1 = ubicación de los ", nrow(codones_top),
                      " codones más representativos | resto = nube general (gris) + familia resaltada"),
    x = paste0("Eje 1 (", round(var_exp_coa[1], 1), "% inercia)"),
    y = paste0("Eje 2 (", round(var_exp_coa[2], 1), "% inercia)")
  ) +
  theme_bw(base_size = 9) +
  theme(
    strip.text       = element_text(face = "bold", size = 7),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle    = element_text(size = 8, hjust = 0.5, color = "grey30"),
    axis.text        = element_text(size = 6)
  )

n_filas_facet <- ceiling(length(niveles_panel) / 5)
ggsave("COA_codones_facet_por_familia.pdf", p_facet,
       width = 15, height = n_filas_facet * 2.8, dpi = 600, limitsize = FALSE)
print(p_facet)

cat("Generado: COA_codones_facet_por_familia.pdf (",
    length(niveles_panel), "paneles: 1 de codones +", n_fam, "de familias, ejes comparables)\n")
