from Bio import Entrez
import os

# Tu email obligatorio para NCBI
Entrez.email = "nsanrom@correo.ugr.es"

# Carpeta de salida
output_folder = "HSP70_Brassica_rapa"
os.makedirs(output_folder, exist_ok=True)

# HSP70 / Heat Shock Protein 70 en Brassica rapa (homólogos típicos)
hsp70_genes = {
    "HSP70-1": "LOC103850120",
    "HSP70-2": "LOC103864330",
    "HSP70-3": "LOC103871910",
    "HSP70-4": "LOC103879120",
    "BiP-like": "LOC103845670",
    "mtHSP70": "LOC103858440",
    "cpHSP70": "LOC103866780"
}

for gene, locus in hsp70_genes.items():

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

print("Descarga HSP70 completada en Brassica rapa")