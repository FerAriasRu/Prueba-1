# Prueba-1

# 1. Acceder a genomas desde SRA/UCSC
   # Descargar genoma desde NCBI
   datasets download genome accession GCA_033216535.1
   # Descargar genoma desde USCS, contiene anotaciones realizadas por ANGUS Y XENOREFGENE
   wget --timestamping -m -nH -x --cut-dirs=6 -e robots=off -np -k     --reject "index.html*" -P "GCA_033216535.1"        https://hgdownload.soe.ucsc.edu/hubs/GCA/033/216/535/GCA_033216535.1/
   # En este caso se uso el genoma de USCS, con el fin de determinar si es diferente al de NCBI, se compraron los genomas con MuMmer con el siguiente comando:
   nucmer --prefix=comparacion GCA_033216535.1/GCA_033216535.1.fa ncbi/ncbi_dataset/data/GCA_033216535.1/GCA_033216535.1_TgRH_pasteur_genomic.fna
   # El archivo de salida en este caso es comparacion.delta
   
# 2. Secuencias CDS
   # Se debe instalar el programa  gffread
   conda install -c bioconda gffread
   #Para generar los CDS se usa el comando
   gffread ../GCA_033216535.1.fa -o CDS.gff -g ../GCA_033216535.1.fa -x CDS.fasta
   # primero debe llamarse el archivo fasta que contiene el genoma
   # -o es el archivo de salida en formato gff (Anotado)
   # -g nuevamente es el genoma 
   # -x archivo de salida en formato fasta

# 3. Secuencias de proteinas
  # Se utiliza el mismo promgara gffread
  gffread  GCA_033216535.1_TgRH_pasteur.augustus.gtf -g ../GCA_033216535.1.fa -y proteins.fasta
  # En este caso se debe llamar el archivo con las anotaciones del genoma. 
  # -g  es el archivo del genoma
  # -y es el archivo de salida que contiene secuencias de proteinas en formato fasta. 
  
# 4. Para comprara los nombres de los cromosomas y contings del archivo del genoma y las anotaciones se pueden usar los siguientes comandos de bash on el fin de  guardar los nombres en archivos txt y luego compararlos 
# generar el archivo txt del genoma 
grep '^>' ../GCA_033216535.1.fa | cut -d' ' -f1 | sed 's/^>//' | sort > fasta_ids.txt
# generar el archivo txt de la anotacion 
cut -f1 GCA_033216535.1_TgRH_pasteur.augustus.gtf | grep -v '^#' | sort -u > gtf_ids.txt
#comparar 
#  diferencia entre archivos 
diff diff fasta_ids.txt gtf_ids.txt 

# Información general del genoma


 
