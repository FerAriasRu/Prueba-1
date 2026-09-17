import subprocess
import re
from collections import defaultdict


# ==========================================
# ARCHIVOS
# ==========================================

FASTA = "../../GCA_033216535.1.fa"
GTF = "../GCA_033216535.1_TgRH_pasteur.augustus.gtf"

FASTA_CROMOSOMAS = "cromosomas.fasta"
SALIDA = "tabla_nueva3"


# ==========================================
# 1. SELECCIONAR LOS 13 CROMOSOMAS CON SEQKIT
# ==========================================

subprocess.run([
    "seqkit",
    "grep",
    "-r",
    "-p", "^CM",
    FASTA,
    "-o", FASTA_CROMOSOMAS
], check=True)


# ==========================================
# 2. OBTENER LONGITUD Y GC% CON SEQKIT
# ==========================================

resultado = subprocess.run(
    [
        "seqkit",
        "fx2tab",
        "-n",
        "-l",
        "-g",
        FASTA_CROMOSOMAS
    ],
    capture_output=True,
    text=True,
    check=True
)


estadisticas = {}

for linea in resultado.stdout.strip().split("\n"):

    cromosoma, longitud, gc = linea.split("\t")

    estadisticas[cromosoma] = {
        "longitud": int(longitud),
        "gc": float(gc)
    }


# ==========================================
# 3. CONTAR BASES AMBIGUAS (N)
# ==========================================

n_count = defaultdict(int)

with open(FASTA_CROMOSOMAS) as archivo:

    cromosoma = None

    for linea in archivo:

        linea = linea.strip()

        if linea.startswith(">"):

            cromosoma = linea[1:].split()[0]

        else:

            n_count[cromosoma] += linea.upper().count("N")


# ==========================================
# 4. CONTAR GENES ÚNICOS POR CROMOSOMA
# ==========================================

genes = defaultdict(set)

with open(GTF) as archivo:

    for linea in archivo:

        if linea.startswith("#"):
            continue

        campos = linea.rstrip().split("\t")

        cromosoma = campos[0]

        # Solo cromosomas que empiezan por CM
        if not cromosoma.startswith("CM"):
            continue

        atributos = campos[8]

        # Extraer gene_id
        encontrado = re.search(
            r'gene_id\s+"([^"]+)"',
            atributos
        )

        if encontrado:

            gene_id = encontrado.group(1)

            # Guardamos el gene_id en un conjunto
            # para no contar el mismo gen varias veces
            genes[cromosoma].add(gene_id)


# ==========================================
# 5. CREAR LA TABLA
# ==========================================


with open(SALIDA, "w") as salida:

    # Encabezados
    salida.write(
        "Cromosoma     "
        "Longitud_bp     "
        "N     "
        "GC_porcentaje     "
        "Genes     "
        "Genes_por_Mb\n"
    )

    # Datos
    for cromosoma in estadisticas:
        longitud = estadisticas[cromosoma]["longitud"]
        gc = estadisticas[cromosoma]["gc"]
        n = n_count[cromosoma]
        numero_genes = len(genes[cromosoma])

        # Genes por megabase
        densidad_genica = numero_genes / (longitud / 1_000_000)

        salida.write(
            f"{cromosoma}     "
            f"{longitud}     "
            f"{n}     "
            f"{gc:.2f}     "
            f"{numero_genes}     "
            f"{densidad_genica:.2f}\n"
        )


print("Tabla creada:", SALIDA)
