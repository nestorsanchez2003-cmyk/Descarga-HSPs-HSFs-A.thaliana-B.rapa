import os
import time
import re
from Bio import Entrez
from Bio import SeqIO

def descargar_hsp40_completo_b_rapa():
    # --- CONFIGURACIÓN ---
    Entrez.email = "tu_correo@ejemplo.com" 
    especie = "Brassica rapa"
    familia_nombre = "HSP40"
    
    # Términos de búsqueda exhaustivos para HSP40/DnaJ
    # Esta familia es gigantesca en plantas.
    terminos = [
        'HSP40', 'heat shock protein 40', 'DnaJ', 'DnaJ-like', 
        'J-domain protein', 'J-protein', 'DnaJA', 'DnaJB', 'DnaJC'
    ]
    
    nombre_carpeta = f"{familia_nombre}_{especie.replace(' ', '_')}_Oficial"
    if not os.path.exists(nombre_carpeta):
        os.makedirs(nombre_carpeta)

    print(f"--- INICIANDO BÚSQUEDA MASIVA DE {familia_nombre}/DnaJ EN {especie} ---")
    print("Nota: Esta familia es muy numerosa, el proceso puede tardar varios minutos.")
    
    ids_genes = set()

    # PASO 1: Identificar todos los Genes en la base de datos 'gene'
    for t in terminos:
        query = f'("{t}"[Gene/Protein Name]) AND "{especie}"[Organism]'
        try:
            # Aumentamos retmax porque HSP40 tiene muchísimos miembros
            handle = Entrez.esearch(db="gene", term=query, retmax=500)
            record = Entrez.read(handle)
            handle.close()
            
            nuevos_ids = record["IdList"]
            if nuevos_ids:
                print(f"Buscando '{t}': Encontrados {len(nuevos_ids)} genes.")
                ids_genes.update(nuevos_ids)
            time.sleep(0.3)
        except Exception as e:
            print(f"  -> Error buscando {t}: {e}")

    print(f"\nTotal de genes únicos detectados: {len(ids_genes)}")
    
    if not ids_genes:
        print("No se encontraron genes.")
        return

    # PASO 2: Obtener las secuencias de ARNm (RefSeq) asociadas
    secuencias_descargadas = []
    vistos_accession = set()

    print("Accediendo a las secuencias de nucleótidos (ARNm)...")
    
    # Ordenamos los IDs para tener un proceso estable
    lista_ids = sorted(list(ids_genes))

    for i, g_id in enumerate(lista_ids, 1):
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
                
                # Filtrado de calidad estricto
                if seq_record.id.startswith(("XM_", "NM_")):
                    vistos_accession.add(seq_record.id)
                    
                    # Intentar extraer nombre del gen o locus tag
                    gene_symbol = "DnaJ"
                    for feat in seq_record.features:
                        if feat.type == "gene":
                            gene_symbol = feat.qualifiers.get("gene", feat.qualifiers.get("locus_tag", ["DnaJ"]))[0]
                            break
                    
                    nombre_archivo = f"{gene_symbol}_{seq_record.id}.fasta".replace(".", "_")
                    nombre_archivo = re.sub(r'[\\/*?:"<>| ]', "_", nombre_archivo)
                    
                    ruta_salida = os.path.join(nombre_carpeta, nombre_archivo)
                    with open(ruta_salida, "w") as f:
                        SeqIO.write(seq_record, f, "fasta")
                    
                    secuencias_descargadas.append(seq_record)
                    
                    if i % 10 == 0: # Imprimir progreso cada 10 genes para no saturar la consola
                        print(f"Progreso: {i}/{len(ids_genes)} genes procesados...")
            
            time.sleep(0.35) 

        except Exception as e:
            # Errores puntuales de conexión son normales en volúmenes grandes
            pass

    # PASO 3: Guardar el conjunto total
    if secuencias_descargadas:
        archivo_final = os.path.join(nombre_carpeta, f"{familia_nombre}_CONJUNTO_TOTAL.fasta")
        with open(archivo_final, "w") as f_all:
            SeqIO.write(secuencias_descargadas, f_all, "fasta")
        
        print(f"\n--- PROCESO COMPLETADO ---")
        print(f"Se han descargado {len(secuencias_descargadas)} secuencias de HSP40/DnaJ.")
        print(f"Carpeta: {os.path.abspath(nombre_carpeta)}")
    else:
        print("No se encontraron secuencias válidas.")

if __name__ == "__main__":
    descargar_hsp40_completo_b_rapa()