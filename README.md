# Prueba-1

## 1. Acceder a genomas desde SRA/UCSC
   Descargar genoma desde NCBI
   ### datasets download genome accession GCA_033216535.1
   Descargar genoma desde USCS, contiene anotaciones realizadas por ANGUS Y XENOREFGENE
   ### wget --timestamping -m -nH -x --cut-dirs=6 -e robots=off -np -k     --reject "index.html*" -P "GCA_033216535.1"        https://hgdownload.soe.ucsc.edu/hubs/GCA/033/216/535/GCA_033216535.1/
   En este caso se uso el genoma de USCS, con el fin de determinar si es diferente al de NCBI, se compraron los genomas con MuMmer con el siguiente comando:
   ### nucmer --prefix=comparacion GCA_033216535.1/GCA_033216535.1.fa ncbi/ncbi_dataset/data/GCA_033216535.1/GCA_033216535.1_TgRH_pasteur_genomic.fna
   El archivo de salida en este caso es comparacion.delta
   
## 2. Secuencias CDS
   Se debe instalar el programa  gffread
   ### conda install -c bioconda gffread
   Para generar los CDS se usa el comando
   ### gffread ../GCA_033216535.1.fa -o CDS.gff -g ../GCA_033216535.1.fa -x CDS.fasta
   primero debe llamarse el archivo fasta que contiene el genoma<br>
   -o es el archivo de salida en formato gff (Anotado)<br>
   -g nuevamente es el genoma<br> 
   -x archivo de salida en formato fasta<br>

## 3. Secuencias de proteinas
  Se utiliza el mismo promgara gffread
  ### gffread  GCA_033216535.1_TgRH_pasteur.augustus.gtf -g ../GCA_033216535.1.fa -y proteins.fasta
  En este caso se debe llamar el archivo con las anotaciones del genoma.<br> 
  -g  es el archivo del genoma<br>
  -y es el archivo de salida que contiene secuencias de proteinas en formato fasta.<br> 
  
## 4. Para comprara los nombres de los cromosomas y contings del archivo del genoma y las anotaciones se pueden usar los siguientes comandos de bash on el fin de  guardar los nombres en archivos txt y luego compararlos 
  generar el archivo txt del genoma<br> 
  ### grep '^>' ../GCA_033216535.1.fa | cut -d' ' -f1 | sed 's/^>//' | sort > fasta_ids.txt<br>
  generar el archivo txt de la anotacion<br> 
  ### cut -f1 GCA_033216535.1_TgRH_pasteur.augustus.gtf | grep -v '^#' | sort -u > gtf_ids.txt<br>
  comparar diferencia entre archivos<br> 
  ### diff diff fasta_ids.txt gtf_ids.txt<br> 

## 5. Información general del genoma
  En este punto se desea una tabla que mencione:
   Cromosoma;
   longitud;
   cantidad de bases ambiguas (N);
   porcentaje de GC;
   número de genes;
   densidad génica, expresada como genes por megabase.

   Por lo cual se utilizo el strip [tabla.py](tabla.py), el cual necesita como entrada la anotacion del genoma y el genoma. Entrega un archivo que contiene la tabla con las 6 columnas y otro con la secunecia fasta de los cromosomas. 

## 6. Distribución del contenido GC en ventanas de 50Kb, se debe contruir un archivo con las ventanas y su respectiva información en este caso la salida sera [tabla_50kb.tsv](tabla_50kb.tsv)


1. Cromosomas y longitd (Salen secuencias de scafolds)

### infoseq -sequence ../GCA_033216535.1.fa -only -name -length -outfile chrom.sizes
 Se corta con head y sale el archivo chrom.sizes2
 
2. Se debe convertir a un formato bet 
 ### awk '{print $1"\t"$2}' chrom.sizes2 > chrom.sizes3

3. Se hace un archivo de ventanas (se debe quitar  la primera fila del archivo chrom.sizes3  que es el nombre y la longitud 
### bedtools makewindows -g chrom.sizes3 -w 50000 > windows_50kb.bed^C

4. Luego se usa bedtools para calcular lo que se necesita 
### bedtools nuc  -fi ../GCA_033216535.1.fa -bed windows_50kb.bed  > nuc_50kb.txt

5. Se calcula el % de bases ambiguas usando
   
awk 'BEGIN{OFS="\t"}
NR==1 {
    print "chrom","start","end","GC_percent","N_count","N_percent"
    next
}
{
    print $1,$2,$3,$5*100,$10,($10/$12)*100
}' nuc_50kb.txt > composicion_50kb.tsv

6. Para número de genes por ventana densidad génica (Se debe saber como el archivo de anotaciones identifica los genes, en este caso

### grep -v '^#' GCA_033216535.1_TgRH_pasteur.augustus.gtf | head -3


* 6.1 Se crea un bed genes. Para cada gene_id, guarda el inicio mínimo y el final máximo.

Por ejemplo, si el GTF tiene:

g1    exon    1322-2538
g1    CDS     1340-2500
g1    exon    5000-6000

 awk 'BEGIN{OFS="\t"}{gene=""
    if (match($0,/gene_id "[^"]+"/)) {
        gene=substr($0,RSTART+9,RLENGTH-10)
    }if (gene!="") {
        key=$1 SUBSEP gene 
        if (!(key in min)) {
            min[key]=$4
            max[key]=$5
            chr[key]=$1
        } else {
            if ($4 < min[key]) min[key]=$4
            if ($5 > max[key]) max[key]=$5
        }
    }
}
END {
    for (key in min) {
        split(key,a,SUBSEP)
        print chr[key], min[key]-1, max[key], a[2]
    }
}' GCA_033216535.1_TgRH_pasteur.augustus.gtf > genes.bed

* 6.2. Contar genes por ventana  (La cuarta columna es la cantidad de genes)

### bedtools intersect -a windows_50kb.bed  -b genes.bed -c > windows_genes.bed

* 6.3. Calcular la densidad genica (La ultima columna es la genes/MB)

awk 'BEGIN{OFS="\t"} {
    print $1,$2,$3,$4,$4/0.05
}' windows_genes.bed > gene_density_50kb.tsv


7. Crear la tabla con toda la informacion chrom	Cromosoma
- start	Inicio de ventana
- end	Final de ventana
- GC_percent	Porcentaje G+C
- N_count	Número de bases ambiguas
- N_percent	Porcentaje de bases ambiguas
- genes	Número de genes
- genes_per_Mb	Densidad génica

awk 'BEGIN{OFS="\t"}
NR==FNR {
    if(FNR>1)
        genes[$1 FS $2 FS $3]=$4
    next
}
FNR==1 {
    print "chrom","start","end","GC_percent","N_count","N_percent","genes","genes_per_Mb"
    next
}
{
    key=$1 FS $2 FS $3
    if(key in genes)
        print $1,$2,$3,$5*100,$10,($10/$12)*100,genes[key],genes[key]/0.05
}' gene_density_50kb.tsv nuc_50kb.txt > tabla_50kb.tsv

## 7. Graficas de %GC y densidad de genes por  cada cromosoma. 
se utilizo el scrip de R [Grafica_GC_desnidad_genica_unacolumna.R](Grafica_GC_desnidad_genica_unacolumna.R) El cual guarda las graficas en una carpeta tanto indivial como en conjunto, la salida final es una imagen en pdf que contiene todas las graficas [ALL_chromosomes_GC_gene_density.pdf](ALL_chromosomes_GC_gene_density.pdf)

## 8. Grafias de %GC entre cromosomas 

Para realizar la comparacion de %CG entre cromosomas con respecto al promedio global se utilizo el siguiente script [Grafica_comparacion_CG_cromosomas.R](Grafica_comparacion_CG_cromosomas.R)

## 9.PCA segun frecuencia de aa











