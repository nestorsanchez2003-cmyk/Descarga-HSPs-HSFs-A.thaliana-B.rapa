from Bio import Entrez

Entrez.email = "nsanrom@correo.ugr.es"

hsp100_genes = {
     "HSP90-1": "AT5G52640",
    "HSP90-2": "AT5G56030",
    "HSP90-3": "AT5G56010",
    "HSP90-4": "AT5G56000",
    "HSP90-5": "AT4G24190",
    "HSP90-6": "AT2G04030",
    "HSP90-7": "AT4G12400"
}

for gene, agi in hsp100_genes.items():

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

    with open(f"{gene}.fasta", "w") as f:
        f.write(fasta)

    print(f"{gene} descargado correctamente")

print("Descarga completada")