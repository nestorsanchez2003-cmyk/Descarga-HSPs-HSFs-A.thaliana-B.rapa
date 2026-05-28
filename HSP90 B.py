from Bio import Entrez
import os

# Tu email obligatorio para NCBI
Entrez.email = "nsanrom@correo.ugr.es"

# Carpeta de salida
output_folder = "HSP90_Brassica_rapa"
os.makedirs(output_folder, exist_ok=True)

# HSP90 / Heat Shock Protein 90 en Brassica rapa (homólogos conocidos)
hsp90_genes = {
    "HSP90-1": "LOC103852910",
    "HSP90-2": "LOC103865430",
    "HSP90-3": "LOC103871120",
    "HSP90-4": "LOC103879550"
}

for gene, locus in hsp90_genes.items():

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

print("Descarga HSP90 completada en Brassica rapa")