from Bio import Entrez
import os

# Tu email obligatorio para NCBI
Entrez.email = "tu_correo@ejemplo.com"

# Carpeta de salida
output_folder = "sHSP_Brassica_rapa"
os.makedirs(output_folder, exist_ok=True)

# sHSP / small Heat Shock Proteins en Brassica rapa (representativo)
shsp_genes = {
    "HSP17.6": "LOC103845910",
    "HSP17.8": "LOC103852140",
    "HSP18.1": "LOC103861300",
    "HSP18.2": "LOC103867450",
    "HSP21": "LOC103874120",
    "HSP22": "LOC103879910",
    "HSP23.5": "LOC103883210",
    "HSP23.6": "LOC103888450",
    "HSP26.5": "LOC103891120"
}

for gene, locus in shsp_genes.items():

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

print("Descarga sHSP completada en Brassica rapa")