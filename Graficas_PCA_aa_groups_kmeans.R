# ============================================================
# PCA de propiedades fisicoquímicas (en vez de aa crudos)
# Toxoplasma gondii RH Pasteur - proteoma predicho
# ============================================================
# Input esperado: data.frame `df` con columnas:
# protein_id, length, A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,
# chromosome, gene_id, start, end, strand, CDS_length
# (los 20 aminoácidos como % relativo, ej. df$A = 7.3)

library(dplyr)
library(ggplot2)
library(factoextra) # para fviz_pca / visualizacion opcional

# Ajusta esta línea a como cargues tu tabla real:
df <- tabla_maestra

# ------------------------------------------------------------
# 1. Índices fisicoquímicos derivados de la composición de aa
# ------------------------------------------------------------
# Valores de hidropatía Kyte-Doolittle (para GRAVY)
kd <- c(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5,
        G=-0.4, H=-3.2, I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8,
        P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)

aa_cols <- c("A","C","D","E","F","G","H","I","K","L",
             "M","N","P","Q","R","S","T","V","W","Y")

df <- df %>%
  mutate(
    # GRAVY: hidrofobicidad promedio ponderada por frecuencia (%/100)
    GRAVY = rowSums(sapply(aa_cols, function(a) (get(a)/100) * kd[a])),
    
    # Carga neta aproximada (residuos básicos - ácidos)
    net_charge = (K + R) - (D + E),
    
    # Fracción de residuos cargados (proxy de superficie polar)
    charged_frac = K + R + D + E,
    
    # Fracción hidrofóbica "core" (A,V,L,I,M,F,W,C)
    hydrophobic_frac = A + V + L + I + M + F + W + C,
    
    # Aromaticidad (F,W,Y,H)
    aromatic_frac = F + W + Y + H,
    
    # Baja complejidad / desorden (S,P,Q,N) — típico de efectores secretados
    low_complexity_frac = S + P + Q + N,
    
    # Índice alifático simplificado (A,V,L,I ponderados; escala Ikai simplificada)
    aliphatic_index = A + 2.9 * V + 3.9 * (L + I)
  )

# ------------------------------------------------------------
# 2. PCA sobre las variables derivadas (no sobre los 20 aa crudos)
# ------------------------------------------------------------
feat_cols <- c("GRAVY", "net_charge", "charged_frac", "hydrophobic_frac",
               "aromatic_frac", "low_complexity_frac", "aliphatic_index")

pca_data <- df %>% select(all_of(feat_cols)) %>% scale()

pca_res <- prcomp(pca_data, center = FALSE, scale. = FALSE) # ya escalado arriba

var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
var_exp

scores <- as.data.frame(pca_res$x[, 1:2])
colnames(scores) <- c("PC1", "PC2")
scores$protein_id <- df$protein_id
scores$length <- df$length
scores$GRAVY <- df$GRAVY
scores$net_charge <- df$net_charge

# ------------------------------------------------------------
# 3. Clustering sobre el espacio PCA (para explorar grupos)
# ------------------------------------------------------------
set.seed(42)

# Método del codo para elegir k
wss <- sapply(2:10, function(k) {
  kmeans(pca_res$x[, 1:2], centers = k, nstart = 10)$tot.withinss
})
# plot(2:10, wss, type="b") # revisa el "codo" antes de fijar k

k_final <- 4 # AJUSTA según el codo/silhouette que observes
km <- kmeans(pca_res$x[, 1:2], centers = k_final, nstart = 25)
scores$cluster <- factor(km$cluster)



# ------------------------------------------------------------
# 4. Visualización
# ------------------------------------------------------------
ggplot(scores, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(alpha = 0.6, size = 1.2) +
  labs(
    title = "PCA de índices fisicoquímicos derivados de composición de aa",
    subtitle = paste0("PC1 (", var_exp[1], "% var) vs PC2 (", var_exp[2], "% var) | n = ", nrow(df)),
    x = paste0("PC1 (", var_exp[1], "% de varianza)"),
    y = paste0("PC2 (", var_exp[2], "% de varianza)")
  ) +
  theme_minimal()

# Colorea también por GRAVY y por carga neta para ver si los clusters
# corresponden a proteínas hidrofóbicas (membrana) vs básicas (DNA-binding, ej. AP2)
ggplot(scores, aes(x = PC1, y = PC2, color = GRAVY)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_gradient2(low = "blue", mid = "grey90", high = "red", midpoint = 0) +
  labs(title = "PCA coloreado por hidrofobicidad (GRAVY)",
       x = paste0("PC1 (", var_exp[1], "% de varianza)"),
       y = paste0("PC2 (", var_exp[2], "% de varianza)") ) +
  theme_minimal()

ggplot(scores, aes(x = PC1, y = PC2, color = net_charge)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_gradient2(low = "blue", mid = "grey90", high = "red", midpoint = 0) +
  labs(title = "PCA coloreado por carga neta (posible señal de proteínas AP2/DNA-binding)") +
  theme_minimal()

# ------------------------------------------------------------
# 5. Guardar tabla enriquecida para cruzar luego con InterProScan
# ------------------------------------------------------------
write.csv(df, "proteins_with_physicochemical_features.csv", row.names = FALSE)
write.csv(scores, "pca_scores_with_clusters.csv", row.names = FALSE)

# ------------------------------------------------------------
# NOTA IMPORTANTE:
# Estos clusters reflejan clases estructurales/fisicoquímicas
# (membrana, básicas, desordenadas, etc.), NO función biológica
# específica en el sentido de GO/Pfam. Una vez tengas la salida
# de InterProScan, cruza `protein_id` con las anotaciones de
# dominios/GO y usa esas categorías (no aa composition) para
# colorear el PCA y validar si los clusters fisicoquímicos
# coinciden con alguna familia funcional conocida.
# ============================================================