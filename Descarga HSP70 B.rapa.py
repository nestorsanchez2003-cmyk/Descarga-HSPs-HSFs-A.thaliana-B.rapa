import os
import time
import re
from Bio import Entrez
from Bio import SeqIO

def descargar_hsp70_final_science():
    # --- CONFIGURACIÓN ---
    Entrez.email = "tu_correo@ejemplo.com" 
    familia_genica = "HSP70"
    especie = "Brassica rapa"
    
    nombre_carpeta = f"{familia_genica}_{especie.replace(' ', '_')}"
    if not os.path.exists(nombre_carpeta):
        os.makedirs(nombre_carpeta)

    # --- ESTRATEGIA DE BÚSQUEDA REFINADA ---
    # Añadimos "biomol rna" y "protein coding" para evitar los XR_ (no codificantes)
    query = (
        f'("{familia_genica}"[Title] OR "heat shock protein 70"[Title] OR "heat shock 70 kDa protein"[Title]) '
        f'AND "{especie}"[Organism] '
        f'AND "biomol rna"[Properties] '
        f'AND 1000:5000[SLEN]'
    )

    print(f"Buscando secuencias codificantes de HSP70 en {especie}...")

    try:
        handle = Entrez.esearch(db="nucleotide", term=query, retmax=100)
        record = Entrez.read(handle)
        handle.close()

        id_list = record["IdList"]
        print(f"Se han identificado {len(id_list)} secuencias codificantes potenciales.")

        lista_records = []
        
        for i, seq_id in enumerate(id_list, 1):
            try:
                # Descarga en GenBank para metadatos
                handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
                seq_record = SeqIO.read(handle, "genbank")
                handle.close()

                # --- FILTRO CIENTÍFICO: Solo queremos mRNAs (XM_ o NM_) ---
                # Esto descarta XR_ (ncRNAs) y fragmentos de DNA genómico
                if not (seq_record.id.startswith("XM_") or seq_record.id.startswith("NM_")):
                    print(f"[{i}/{len(id_list)}] Saltando {seq_record.id} (No es mRNA codificante)")
                    continue

                # --- EXTRACCIÓN DE NOMBRE ---
                gene_symbol = ""
                for feature in seq_record.features:
                    if feature.type == "gene":
                        gene_symbol = feature.qualifiers.get("gene", [""])[0]
                        if not gene_symbol:
                            gene_symbol = feature.qualifiers.get("locus_tag", [""])[0]
                        break
                
                if not gene_symbol:
                    gene_symbol = "HSP70_like"

                # Limpiar nombre: Ej. LOC103836222_XM_033277924
                nombre_base = re.sub(r'[\\/*?:"<>| ]', "_", gene_symbol)
                nombre_archivo = f"{nombre_base}_{seq_record.id}.fasta"
                
                # --- GUARDAR INDIVIDUAL ---
                ruta_fasta = os.path.join(nombre_carpeta, nombre_archivo)
                with open(ruta_fasta, "w") as f_out:
                    SeqIO.write(seq_record, f_out, "fasta")
                
                lista_records.append(seq_record)
                print(f"[{i}/{len(id_list)}] Guardado: {nombre_archivo}")
                
                time.sleep(0.4) 

            except Exception as e:
                print(f"Error con ID {seq_id}: {e}")

        # --- GUARDAR CONJUNTO ---
        if lista_records:
            archivo_total = os.path.join(nombre_carpeta, f"{familia_genica}_{especie.replace(' ', '_')}_CONJUNTO.fasta")
            with open(archivo_total, "w") as f_all:
                SeqIO.write(lista_records, f_all, "fasta")
            
            print(f"\nProceso terminado. Total de secuencias codificantes: {len(lista_records)}")
            print(f"Carpeta: {nombre_carpeta}")
        else:
            print("No se pudieron procesar secuencias.")

    except Exception as e:
        print(f"Error crítico: {e}")

if __name__ == "__main__":
    descargar_hsp70_final_science()