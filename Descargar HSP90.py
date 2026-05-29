import os
import time
import re
from Bio import Entrez
from Bio import SeqIO

def descargar_hsp90_completo_b_rapa():
    # --- CONFIGURACIÓN ---
    Entrez.email = "tu_correo@ejemplo.com" 
    especie = "Brassica rapa"
    familia_nombre = "HSP90"
    
    # Términos de búsqueda exhaustivos para HSP90 en plantas
    # HSP81, HSP82 y HSP83 son nombres alternativos comunes para HSP90
    terminos = [
        'HSP90', 'heat shock protein 90', 'HSP81', 'HSP82', 'HSP83',
        'HSP90-1', 'HSP90-2', 'HSP90-3', 'HSP90-4', 'HSP90.5', 'HSP90.7'
    ]
    
    nombre_carpeta = f"{familia_nombre}_{especie.replace(' ', '_')}_Oficial"
    if not os.path.exists(nombre_carpeta):
        os.makedirs(nombre_carpeta)

    print(f"--- INICIANDO BÚSQUEDA EXHAUSTIVA DE {familia_nombre} EN {especie} ---")
    
    ids_genes = set()

    # PASO 1: Identificar todos los Genes en la base de datos 'gene'
    for t in terminos:
        query = f'("{t}"[Gene/Protein Name]) AND "{especie}"[Organism]'
        try:
            handle = Entrez.esearch(db="gene", term=query, retmax=100)
            record = Entrez.read(handle)
            handle.close()
            
            nuevos_ids = record["IdList"]
            if nuevos_ids:
                print(f"Buscando '{t}': Encontrados {len(nuevos_ids)} genes potenciales.")
                ids_genes.update(nuevos_ids)
            time.sleep(0.3)
        except Exception as e:
            print(f"Error buscando {t}: {e}")

    print(f"\nTotal de genes únicos detectados: {len(ids_genes)}")
    
    if not ids_genes:
        print("No se encontraron genes. Verifica los términos o la especie.")
        return

    # PASO 2: Obtener las secuencias de ARNm (RefSeq) asociadas a esos Genes
    secuencias_descargadas = []
    vistos_accession = set()

    print("Accediendo a las secuencias de nucleótidos oficiales...")
    
    for i, g_id in enumerate(ids_genes, 1):
        try:
            # Vinculamos el ID de gen con sus ARNm (XM_ o NM_)
            link_handle = Entrez.elink(dbfrom="gene", db="nuccore", id=g_id, term="srcdb_refseq[prop] AND mRNA[filter]")
            link_record = Entrez.read(link_handle)
            link_handle.close()

            nuccore_ids = []
            if link_record[0]["LinkSetDb"]:
                for link in link_record[0]["LinkSetDb"][0]["Link"]:
                    nuccore_ids.append(link["Id"])

            for n_id in nuccore_ids:
                handle = Entrez.efetch(db="nuccore", id=n_id, rettype="gb", retmode="text")
                seq_record = SeqIO.read(handle, "genbank")
                handle.close()

                if seq_record.id in vistos_accession:
                    continue
                
                # Filtrar solo por secuencias de calidad (RefSeq: XM_ o NM_)
                if seq_record.id.startswith(("XM_", "NM_")):
                    vistos_accession.add(seq_record.id)
                    
                    # Intentar obtener el nombre del gen para el archivo
                    gene_symbol = "HSP90"
                    for feat in seq_record.features:
                        if feat.type == "gene":
                            gene_symbol = feat.qualifiers.get("gene", feat.qualifiers.get("locus_tag", ["HSP90"]))[0]
                            break
                    
                    nombre_archivo = f"{gene_symbol}_{seq_record.id}.fasta".replace(".", "_")
                    nombre_archivo = re.sub(r'[\\/*?:"<>| ]', "_", nombre_archivo)
                    
                    ruta_salida = os.path.join(nombre_carpeta, nombre_archivo)
                    with open(ruta_salida, "w") as f:
                        SeqIO.write(seq_record, f, "fasta")
                    
                    secuencias_descargadas.append(seq_record)
                    print(f"[{i}/{len(ids_genes)}] Guardado: {nombre_archivo}")
            
            time.sleep(0.4) # Respetar el límite de NCBI

        except Exception as e:
            print(f"Error procesando Gene ID {g_id}: {e}")

    # PASO 3: Guardar el conjunto total
    if secuencias_descargadas:
        archivo_final = os.path.join(nombre_carpeta, f"{familia_nombre}_{especie.replace(' ', '_')}_CONJUNTO.fasta")
        with open(archivo_final, "w") as f_all:
            SeqIO.write(secuencias_descargadas, f_all, "fasta")
        
        print(f"\n--- PROCESO COMPLETADO ---")
        print(f"Se han descargado {len(secuencias_descargadas)} secuencias de HSP90 de alta calidad.")
        print(f"Carpeta: {os.path.abspath(nombre_carpeta)}")
    else:
        print("No se encontraron secuencias de ARNm que cumplan con los filtros de calidad.")

if __name__ == "__main__":
    descargar_hsp90_completo_b_rapa()   
   