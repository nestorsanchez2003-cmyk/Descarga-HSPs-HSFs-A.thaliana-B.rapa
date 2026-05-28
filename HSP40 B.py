from Bio import Entrez
import os

# Tu email obligatorio para NCBI
Entrez.email = "nsanromorreo.ugr.es"

# Carpeta de salida
output_folder = "HSP40_Brassica_rapa"
os.makedirs(output_folder, exist_ok=True)

# HSP40 / DNAJ (co-chaperonas J-domain proteins) en Brassica rapa
hsp40_genes = {
    "DNAJ1": "LOC103846210",
    "DNAJ2": "LOC103853440",
    "DNAJ3": "LOC103861980",
    "DNAJ4": "LOC103872150",
    "DNAJ5": "LOC103879300",
    "DNAJ6": "LOC103884120"
}

for gene, locus in hsp40_genes.items():

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

print("Descarga HSP40 completada en Brassica rapa")