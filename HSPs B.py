from Bio import Entrez
import os

# Tu email obligatorio para NCBI
Entrez.email = "nsanrom@correo.ugr.es"

# Carpeta de salida
output_folder = "HSP100_Brassica_rapa"
os.makedirs(output_folder, exist_ok=True)

# HSP100 / ClpB en Brassica rapa (homólogos putativos)
hsp100_genes = {
    "ClpB-cyt": "LOC103849120",
    "ClpB-m": "LOC103872340",
    "ClpB-c": "LOC103855670",
    "ClpB-like1": "LOC103861210"
}

for gene, locus in hsp100_genes.items():

    print(f"Descargando {gene}...")

    search = Entrez.esearch(
        db="protein",
        term=f"{locus}[Gene Name] AND Brassica rapa[Organism]"
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

print("Descarga HSP100 completada en Brassica rapa")