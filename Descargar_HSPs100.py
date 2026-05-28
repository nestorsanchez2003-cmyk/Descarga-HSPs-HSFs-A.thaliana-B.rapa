from Bio import Entrez

Entrez.email = "nsanrom@correo.ugr.es"

hsp100_genes = {
    "HSP101": "AT1G74310",
    "ClpB2": "AT5G15450",
    "ClpB3": "AT1G74360",
    "ClpB4": "AT2G25140",
    "ClpB5": "AT5G22830"
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