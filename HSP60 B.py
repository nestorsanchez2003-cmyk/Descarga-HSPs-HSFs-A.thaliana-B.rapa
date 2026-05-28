from Bio import Entrez
import os

# IMPORTANTE: pon tu email real
Entrez.email = "nsanrom@correo.ugr.es"

# Carpeta de salida por especie
output_folder = "HSP60_Brassica_rapa"
os.makedirs(output_folder, exist_ok=True)

# HSP60 / CPN60 en Brassica rapa (homólogos conocidos en bases de datos)
hsp60_genes = {
    "CPN60A1": "LOC103842930",
    "CPN60A2": "LOC103865421",
    "CPN60B1": "LOC103846710",
    "CPN60B2": "LOC103852330"
}

for gene, locus in hsp60_genes.items():

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

print("Descarga completada para Brassica rapa")