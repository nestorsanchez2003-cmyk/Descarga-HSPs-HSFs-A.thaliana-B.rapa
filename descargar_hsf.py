from Bio import Entrez

# Pon aquí tu correo
Entrez.email = "nsanrom@correo.ugr.es"

hsf_genes = {
    "HSFA1A": "AT4G17750",
    "HSFA1B": "AT5G16820",
    "HSFA1D": "AT1G32330",
    "HSFA1E": "AT3G02990",
    "HSFA2": "AT2G26150",
    "HSFA3": "AT5G03720",
    "HSFA4A": "AT4G18880",
    "HSFA4C": "AT5G45710",
    "HSFA5": "AT4G13980",
    "HSFA6A": "AT5G43840",
    "HSFA6B": "AT3G22830",
    "HSFA7A": "AT3G51910",
    "HSFA7B": "AT3G63350",
    "HSFA8": "AT1G67970",
    "HSFA9": "AT5G54070",
    "HSFB1": "AT4G36990",
    "HSFB2A": "AT5G62020",
    "HSFB2B": "AT4G11660",
    "HSFB3": "AT2G25140",
    "HSFB4": "AT1G46264",
    "HSFC1": "AT3G24520"
}

for gene, agi in hsf_genes.items():

    print(f"Descargando {gene}...")

    search = Entrez.esearch(
        db="gene",
        term=f"{agi}[Gene Name] AND Arabidopsis thaliana[Organism]"
    )

    record = Entrez.read(search)

    if not record["IdList"]:
        print(f"No encontrado: {gene}")
        continue

    gene_id = record["IdList"][0]

    fetch = Entrez.efetch(
        db="gene",
        id=gene_id,
        rettype="fasta",
        retmode="text"
    )

    fasta = fetch.read()

    with open(f"{gene}.fasta", "w") as f:
        f.write(fasta)

    print(f"{gene} descargado correctamente")

print("Todas las secuencias fueron descargadas")