import os
from Bio import Entrez, SeqIO
import time

# 1. Configuración
Entrez.email = "tu_correo@ejemplo.com" # Cambia esto
folder_name = "B_rapa_HSFs_Final_Clean"
# Búsqueda amplia en la base de datos de proteínas
search_query = '("Brassica rapa"[Organism]) AND ("heat shock transcription factor" OR "heat shock factor")'

if not os.path.exists(folder_name):
    os.makedirs(folder_name)

def download_and_clean_combined():
    print("1. Buscando proteínas HSF en NCBI para obtener máxima cobertura...")
    search_handle = Entrez.esearch(db="protein", term=search_query, retmax=200)
    search_record = Entrez.read(search_handle)
    search_handle.close()
    
    protein_ids = search_record["IdList"]
    print(f"Se localizaron {len(protein_ids)} registros. Iniciando mapeo a nucleótidos y limpieza...")

    # Diccionario para evitar duplicados: { 'ID_GEN': SeqRecord }
    curated_genes = {}

    for p_id in protein_ids:
        try:
            # 2. Vincular Proteína -> Nucleótido
            link_handle = Entrez.elink(dbfrom="protein", db="nucleotide", id=p_id, linkname="protein_nucleotide_mgc")
            link_record = Entrez.read(link_handle)
            link_handle.close()

            if not link_record[0]["LinkSetDb"]:
                # Reintento con link genérico si el específico falla
                link_handle = Entrez.elink(dbfrom="protein", db="nucleotide", id=p_id)
                link_record = Entrez.read(link_handle)
                link_handle.close()

            if link_record[0]["LinkSetDb"]:
                nuc_id = link_record[0]["LinkSetDb"][0]["Link"][0]["Id"]
                
                # 3. Descargar el GenBank para extraer metadatos
                fetch_handle = Entrez.efetch(db="nucleotide", id=nuc_id, rettype="gb", retmode="text")
                seq_record = SeqIO.read(fetch_handle, "genbank")
                fetch_handle.close()

                # 4. Identificar el Locus Tag (ID único del gen Bra...)
                locus_tag = None
                for feature in seq_record.features:
                    if feature.type == "gene":
                        if "locus_tag" in feature.qualifiers:
                            locus_tag = feature.qualifiers["locus_tag"][0]
                            break
                
                # Si no tiene Locus Tag, usamos el ID de la secuencia
                unique_key = locus_tag if locus_tag else seq_record.id
                
                # 5. Lógica de limpieza: quedarnos con la versión más larga (mejor anotada)
                if unique_key not in curated_genes:
                    curated_genes[unique_key] = seq_record
                    print(f"Nuevo gen detectado: {unique_key}")
                else:
                    if len(seq_record.seq) > len(curated_genes[unique_key].seq):
                        curated_genes[unique_key] = seq_record
                        print(f"Actualizado {unique_key} por una versión más larga.")
                
            time.sleep(0.2) # Pausa rápida

        except Exception as e:
            continue

    # 6. Guardar los resultados únicos
    print(f"\n--- Guardando resultados finales ---")
    for final_id, record in curated_genes.items():
        file_name = "".join([c for c in final_id if c.isalnum() or c in ('_', '.')])
        file_path = os.path.join(folder_name, f"{file_name}.fasta")
        SeqIO.write(record, file_path, "fasta")

    print(f"\nProceso terminado. De {len(protein_ids)} registros iniciales,")
    print(f"hemos obtenido {len(curated_genes)} genes HSF únicos.")

if __name__ == "__main__":
    download_and_clean_combined()