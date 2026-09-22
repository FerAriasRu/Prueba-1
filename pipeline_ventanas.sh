#!/bin/bash

# ============================================================
# PIPELINE AUTOMATIZADO:
# GC + N + DENSIDAD GÉNICA POR VENTANAS
#
# Uso:
#   ./pipeline_ventanas.sh 20000
#   ./pipeline_ventanas.sh 50000
#
# Genera:
#   resultados_20kb/
#   resultados_50kb/
# ============================================================


# ============================================================
# 1. ARCHIVOS DE ENTRADA
# ============================================================

FASTA="../GCA_033216535.1/GCA_033216535.1.fa"
GTF="../GCA_033216535.1/genes/GCA_033216535.1_TgRH_pasteur.augustus.gtf"


# ============================================================
# 2. TAMAÑO DE VENTANA
# ============================================================

WINDOW=$1

if [ -z "$WINDOW" ]; then

    echo ""
    echo "ERROR: debes indicar el tamaño de la ventana."
    echo ""
    echo "Ejemplos:"
    echo "  ./pipeline_ventanas.sh 20000"
    echo "  ./pipeline_ventanas.sh 50000"
    echo ""

    exit 1
fi


WINDOW_KB=$((WINDOW / 1000))

OUTDIR="resultados_${WINDOW_KB}kb"

mkdir -p "$OUTDIR"


# ============================================================
# 3. COMPROBAR PROGRAMAS
# ============================================================

if ! command -v infoseq >/dev/null 2>&1; then

    echo ""
    echo "ERROR: no se encontró 'infoseq'."
    echo ""
    echo "Activa tu ambiente Conda:"
    echo "  conda activate genomica"
    echo ""
    echo "o instala EMBOSS:"
    echo "  conda install -c bioconda -c conda-forge emboss"
    echo ""

    exit 1
fi


if ! command -v bedtools >/dev/null 2>&1; then

    echo ""
    echo "ERROR: no se encontró 'bedtools'."
    echo ""
    echo "Activa tu ambiente Conda:"
    echo "  conda activate genomica"
    echo ""
    echo "o instala bedtools:"
    echo "  conda install -c bioconda -c conda-forge bedtools"
    echo ""

    exit 1
fi


# ============================================================
# 4. COMPROBAR ARCHIVOS
# ============================================================

if [ ! -f "$FASTA" ]; then

    echo ""
    echo "ERROR: no se encontró el FASTA:"
    echo "$FASTA"
    echo ""

    exit 1
fi


if [ ! -f "$GTF" ]; then

    echo ""
    echo "ERROR: no se encontró el GTF:"
    echo "$GTF"
    echo ""

    exit 1
fi


# ============================================================
# 5. INFORMACIÓN DE LOS CROMOSOMAS
# ============================================================

echo ""
echo "============================================================"
echo "Obteniendo tamaños de los cromosomas"
echo "============================================================"
echo ""

infoseq \
-sequence "$FASTA" \
-only \
-name \
-length \
-outfile "$OUTDIR/chrom.sizes"


# ============================================================
# 6. CREAR chrom.sizes PARA BEDTOOLS
#
# IMPORTANTE:
# - elimina encabezado Name Length
# - conserva solamente cromosomas CM...
# - elimina los contigs JAVJAM...
# ============================================================

awk '
NR > 1 && $1 ~ /^CM/ {
    print $1 "\t" $2
}
' "$OUTDIR/chrom.sizes" > "$OUTDIR/chrom.sizes3"


# Comprobar que chrom.sizes3 no esté vacío

if [ ! -s "$OUTDIR/chrom.sizes3" ]; then

    echo ""
    echo "ERROR: chrom.sizes3 está vacío."
    echo ""
    echo "Contenido de chrom.sizes:"
    cat "$OUTDIR/chrom.sizes"
    echo ""

    exit 1
fi


echo "Cromosomas utilizados:"
cat "$OUTDIR/chrom.sizes3"

echo ""


# ============================================================
# 7. CREAR VENTANAS
# ============================================================

echo "============================================================"
echo "Creando ventanas de ${WINDOW_KB} kb"
echo "============================================================"
echo ""

bedtools makewindows \
-g "$OUTDIR/chrom.sizes3" \
-w "$WINDOW" \
> "$OUTDIR/windows_${WINDOW_KB}kb.bed"


if [ ! -s "$OUTDIR/windows_${WINDOW_KB}kb.bed" ]; then

    echo ""
    echo "ERROR: no se generaron las ventanas."
    echo ""

    exit 1
fi


echo "Primeras ventanas:"
head "$OUTDIR/windows_${WINDOW_KB}kb.bed"

echo ""


# ============================================================
# 8. CALCULAR COMPOSICIÓN NUCLEOTÍDICA
# ============================================================

echo "============================================================"
echo "Calculando composición GC y N"
echo "============================================================"
echo ""

bedtools nuc \
-fi "$FASTA" \
-bed "$OUTDIR/windows_${WINDOW_KB}kb.bed" \
> "$OUTDIR/nuc_${WINDOW_KB}kb.txt"


if [ ! -s "$OUTDIR/nuc_${WINDOW_KB}kb.txt" ]; then

    echo ""
    echo "ERROR: bedtools nuc no produjo resultados."
    echo ""

    exit 1
fi


# ============================================================
# 9. CREAR TABLA DE COMPOSICIÓN
# ============================================================

awk '
BEGIN {
    OFS="\t"
}

NR == 1 {

    print "chrom", "start", "end", "GC_percent", "N_count", "N_percent"

    next
}

{

    print $1, $2, $3, $5 * 100, $10, ($10 / $12) * 100
}

' "$OUTDIR/nuc_${WINDOW_KB}kb.txt" \
> "$OUTDIR/composicion_${WINDOW_KB}kb.tsv"


# ============================================================
# 10. EXTRAER GENES DEL GTF
#
# Cada gene_id se convierte en una sola región:
#
# inicio = coordenada mínima
# final  = coordenada máxima
# ============================================================

echo "============================================================"
echo "Extrayendo genes"
echo "============================================================"
echo ""

awk '
BEGIN {
    OFS="\t"
}

{

    gene=""

    if (match($0, /gene_id "[^"]+"/)) {
        gene = substr($0, RSTART+9, RLENGTH-10)
    }

    if (gene != "") {

        key = $1 SUBSEP gene

        if (!(key in min)) {
            min[key] = $4
            max[key] = $5
            chr[key] = $1
        }
        else {
            if ($4 < min[key]) min[key] = $4
            if ($5 > max[key]) max[key] = $5
        }

    }

}

END {

    for (key in min) {
        split(key, a, SUBSEP)
        print chr[key], min[key]-1, max[key], a[2]
    }

}

' "$GTF" \
> "$OUTDIR/genes.bed"


if [ ! -s "$OUTDIR/genes.bed" ]; then

    echo ""
    echo "ERROR: genes.bed está vacío."
    echo ""

    exit 1
fi


# ============================================================
# 11. CONTAR GENES POR VENTANA
# ============================================================

echo "============================================================"
echo "Contando genes por ventana"
echo "============================================================"
echo ""

bedtools intersect \
-a "$OUTDIR/windows_${WINDOW_KB}kb.bed" \
-b "$OUTDIR/genes.bed" \
-c \
> "$OUTDIR/windows_genes_${WINDOW_KB}kb.bed"


if [ ! -s "$OUTDIR/windows_genes_${WINDOW_KB}kb.bed" ]; then

    echo ""
    echo "ERROR: no se pudo calcular el número de genes."
    echo ""

    exit 1
fi


# ============================================================
# 12. CALCULAR DENSIDAD GÉNICA
#
# genes/Mb =
#
# número de genes /
# longitud de la ventana en Mb
#
# Se utiliza la longitud REAL de cada ventana.
# Esto es importante para las ventanas finales
# de cada cromosoma.
# ============================================================

awk '
BEGIN {
    OFS="\t"
}

{

    window_length = $3 - $2

    genes_per_Mb = $4 / (window_length / 1000000)

    print $1, $2, $3, $4, genes_per_Mb

}

' "$OUTDIR/windows_genes_${WINDOW_KB}kb.bed" \
> "$OUTDIR/gene_density_${WINDOW_KB}kb.tsv"


# ============================================================
# 13. CREAR TABLA FINAL
#
# columnas:
#
# chrom
# start
# end
# GC_percent
# N_count
# N_percent
# genes
# genes_per_Mb
# ============================================================

echo "============================================================"
echo "Construyendo tabla final"
echo "============================================================"
echo ""


awk '
BEGIN {
    OFS="\t"
}


NR == FNR {

    if (FNR > 1)
        genes[$1 FS $2 FS $3] = $4

    next
}


FNR == 1 {

    print "chrom", "start", "end", "GC_percent", "N_count", "N_percent", "genes", "genes_per_Mb"

    next
}


{

    key = $1 FS $2 FS $3

    if (key in genes) {

        window_length = $3 - $2

        genes_per_Mb = genes[key] / (window_length / 1000000)

        print $1, $2, $3, $5 * 100, $10, ($10 / $12) * 100, genes[key], genes_per_Mb

    }

}

' \
"$OUTDIR/gene_density_${WINDOW_KB}kb.tsv" \
"$OUTDIR/nuc_${WINDOW_KB}kb.txt" \
> "$OUTDIR/tabla_${WINDOW_KB}kb.tsv"


# ============================================================
# 14. COMPROBAR TABLA FINAL
# ============================================================

if [ ! -s "$OUTDIR/tabla_${WINDOW_KB}kb.tsv" ]; then

    echo ""
    echo "ERROR: la tabla final está vacía."
    echo ""

    exit 1
fi


# ============================================================
# 15. RESUMEN
# ============================================================

echo ""
echo "============================================================"
echo "PIPELINE TERMINADO"
echo "============================================================"
echo ""

echo "Tamaño de ventana: ${WINDOW_KB} kb"
echo ""

echo "Directorio:"
echo "$OUTDIR"
echo ""

echo "Archivos generados:"
echo ""

ls -lh "$OUTDIR"

echo ""

echo "Número de ventanas:"
tail -n +2 "$OUTDIR/tabla_${WINDOW_KB}kb.tsv" | wc -l

echo ""

echo "Primeras filas de la tabla:"
echo ""

head "$OUTDIR/tabla_${WINDOW_KB}kb.tsv"

echo ""
echo "============================================================"
echo ""
