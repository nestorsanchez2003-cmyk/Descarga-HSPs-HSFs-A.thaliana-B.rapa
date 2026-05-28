import os
from Bio import Entrez, SeqIO
import time

# 1. Configuración de parámetros
Entrez.email = "tu_correo@ejemplo.com"  # RECUERDA CAMBIAR ESTO
folder_name = "A.thaliana HSFs sec"

# Lista oficial de los 21 HSFs de Arabidopsis (Códigos AGI)
agi_hsfs = [
    "AT4G17750", "AT5G16820", "AT1G32330", "AT3G02990",  # Clase A1 (a, b, d, e)
    "AT2G26150", "AT5G03720", "AT4G18880", "AT5G45710",  # A2, A3, A4a, A4c
    "AT4G13960", "AT5G43840", "AT3G22830", "AT3G51910",  # A5, A6a, A6b, A7a
    "AT3G63350", "AT1G67970", "AT5G54070",              # A7b, A8, A9
    "AT4G36990", "AT5G62020", "AT4G11660", "AT4G17140",  # Clase B1, B2a, B2b, B3
    "AT1G46264", "AT3G24520"                             # B4, Clase C1
]

# 2. Crear la carpeta
if not os.path.exists(folder_name):
    os.makedirs(folder_name)
    print(f"Directorio creado: {folder_name}")

def download_by_agi():
    print(f"Iniciando descarga de {len(agi_hsfs)} genes HSF...")
    
    for agi in agi_hsfs:
        try:
            # Buscamos el código AGI específicamente en Arabidopsis y filtramos por mRNA
            search_term = f"{agi}[All Fields] AND Arabidopsis thaliana[Organism] AND mRNA[Filter]"
            
            # Buscamos el ID de NCBI (GI/Accession) para ese código AGI
            search_handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=1)
            search_record = Entrez.read(search_handle)
            search_handle.close()
            
            if search_record["IdList"]:
                ncbi_id = search_record["IdList"][0]
                
                # Descargamos la secuencia en FASTA
                fetch_handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="fasta", retmode="text")
                seq_record = SeqIO.read(fetch_handle, "fasta")
                fetch_handle.close()
                
                # Guardamos con el nombre del código AGI para mayor orden
                file_path = os.path.join(folder_name, f"{agi}.fasta")
                SeqIO.write(seq_record, file_path, "fasta")
                print(f"Éxito: {agi} guardado.")
            else:
                print(f"Aviso: No se encontró secuencia para {agi}")
            
            # Pequeña pausa para no saturar el servidor del NCBI
            time.sleep(0.5)
            
        except Exception as e:
            print(f"Error con {agi}: {e}")

if __name__ == "__main__":
    download_by_agi()
    print("\n--- Proceso completado ---")