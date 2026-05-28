from Bio import Entrez
from Bio import SeqIO
from io import StringIO

# ==============================
# CONFIGURACIÓN
# ==============================

Entrez.email = "nsanrom@correo.ugr.es"

# ==============================
# FAMILIA HSP40 / DNAJ COMPLETA
# Arabidopsis thaliana
# SIN DUPLICADOS
# ==============================

hsp40_genes = {

    "DJA1": "AT1G28210",
    "DJA2": "AT5G22060",
    "DJA3": "AT3G44110",
    "DJA4": "AT1G72440",
    "DJA5": "AT1G76730",
    "DJA6": "AT5G06910",
    "DJA7": "AT5G14140",
    "DJA8": "AT1G80920",
    "DJA9": "AT3G07370",
    "DJA10": "AT5G25550",

    "DJB1": "AT5G05200",
    "DJB2": "AT4G13830",
    "DJB3": "AT2G20560",
    "DJB4": "AT3G13310",
    "DJB5": "AT1G79940",
    "DJB6": "AT3G47660",
    "DJB7": "AT1G03020",
    "DJB8": "AT5G61210",

    "DJC1": "AT4G10060",
    "DJC2": "AT3G47940",
    "DJC3": "AT5G17810",
    "DJC4": "AT5G23040",
    "DJC5": "AT1G74470",
    "DJC6": "AT3G62600",
    "DJC7": "AT4G21130",
    "DJC8": "AT1G80950",
    "DJC9": "AT5G06900",
    "DJC10": "AT4G36040",
    "DJC11": "AT3G08970",
    "DJC12": "AT1G21080",
    "DJC13": "AT3G28740",
    "DJC14": "AT1G22360",
    "DJC15": "AT2G22360",
    "DJC16": "AT5G48030",
    "DJC17": "AT1G80030",
    "DJC18": "AT5G17710",
    "DJC19": "AT4G35780",
    "DJC20": "AT1G73540",
    "DJC21": "AT3G13445",
    "DJC22": "AT5G13870",
    "DJC23": "AT2G33110",
    "DJC24": "AT1G15730",
    "DJC25": "AT4G09340",
    "DJC26": "AT5G03160",
    "DJC27": "AT3G18100",
    "DJC28": "AT5G03030",
    "DJC29": "AT1G55490",
    "DJC30": "AT3G44110",

    "AtJ1": "AT1G80920",
    "AtJ2": "AT4G36040",
    "AtJ3": "AT3G44110",
    "AtJ4": "AT1G28210",
    "AtJ5": "AT5G22060",
    "AtJ6": "AT5G06910",
    "AtJ7": "AT5G14140",
    "AtJ8": "AT1G76730",
    "AtJ9": "AT2G20560",
    "AtJ10": "AT3G07370"
}

# ==============================
# ELIMINAR DUPLICADOS
# ==============================

unique_genes = {}

for gene, agi in hsp40_genes.items():
    if agi not in unique_genes.values():
        unique_genes[gene] = agi

print(f"Genes únicos encontrados: {len(unique_genes)}")

# ==============================
# DESCARGA FASTA PROTEICA
# ==============================

combined_fasta = []

for gene, agi in unique_genes.items():

    print(f"Descargando {gene} ({agi})...")

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

    fasta_data = fetch.read()

    # Guardar FASTA individual
    with open(f"{gene}.fasta", "w") as f:
        f.write(fasta_data)

    # Añadir al FASTA combinado
    combined_fasta.append(fasta_data)

    print(f"{gene} descargado correctamente")

# ==============================
# FASTA COMBINADO
# ==============================

with open("HSP40_Arabidopsis_COMPLETO.fasta", "w") as out_fasta:
    for fasta in combined_fasta:
        out_fasta.write(fasta)

print("===================================")
print("DESCARGA COMPLETADA")
print("FASTA combinado generado:")
print("HSP40_Arabidopsis_COMPLETO.fasta")
print("===================================")