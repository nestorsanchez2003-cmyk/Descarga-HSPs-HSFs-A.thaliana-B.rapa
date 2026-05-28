from Bio import Entrez

Entrez.email = "nsanrom@correo.ugr.es"

hsp100_genes = {
      "HSP70-1": "AT5G02500",
    "HSP70-2": "AT5G02490",
    "HSP70-3": "AT3G09440",
    "HSP70-4": "AT3G12580",
    "HSP70-5": "AT1G16030",
    "HSP70-6": "AT1G56410",
    "HSP70-7": "AT5G49910",
    "BiP1": "AT5G28540",
    "BiP2": "AT5G42020",
    "BiP3": "AT1G09080",
    "cpHSC70-1": "AT4G24280",
    "cpHSC70-2": "AT5G49910",
    "mtHSC70-1": "AT4G37910",
    "mtHSC70-2": "AT5G09590"
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