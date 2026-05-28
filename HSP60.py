from Bio import Entrez
import os

# Tu correo (obligatorio para NCBI)
Entrez.email = "tu_correo@ejemplo.com"

# Carpeta de salida
output_folder = "HSP60 A.t"

# Crear carpeta si no existe
os.makedirs(output_folder, exist_ok=True)

# HSP60 / chaperonina CPN60 (GroEL-like en plantas)
hsp60_genes = {
    "CPN60A1": "AT2G28000",
    "CPN60A2": "AT5G20720",
    "CPN60B1": "AT1G55490",
    "CPN60B2": "AT3G23990",
    "CPN60G1": "AT1G14980",
    "CPN60G2": "AT3G13470"
}

for gene, agi in hsp60_genes.items():

    print(f"Descargando {gene}...")

    search = Entrez.esearch(
        db="protein",
        term=f"{agi}[Gene Name] AND Arabidopsis thaliana[Organism]"
    )

    record = Entrez.read(search)

    if not record["IdList"]:
        print(f"No encontrado: {gene}")
        continue

    protein_id = record["IdList"][0]

    fetch = Entrez.efetch(
        db="protein",
        id=protein_id,
        rettype="fasta",
        retmode="text"
    )

    fasta = fetch.read()

    output_path = os.path.join(output_folder, f"{gene}.fasta")

    with open(output_path, "w") as f:
        f.write(fasta)

    print(f"{gene} descargado correctamente")

print("Descarga de HSP60 completada")