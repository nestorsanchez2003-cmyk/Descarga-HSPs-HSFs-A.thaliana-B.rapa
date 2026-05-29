
from Bio import Entrez
import os
import time

# =========================================
# CONFIGURACIÓN GENERAL
# =========================================

Entrez.email = "nsanrom@correo.ugr.es"

# =========================================
# FUNCIÓN ROBUSTA DE DESCARGA
# =========================================


def descargar_familia_hsp(familia_nombre, genes_dict):

    output_dir = f"{familia_nombre}_Arabidopsis_thaliana"
    os.makedirs(output_dir, exist_ok=True)

    # Eliminar duplicados
    unique_genes = {}

    for gene, agi in genes_dict.items():
        if agi not in unique_genes.values():
            unique_genes[gene] = agi

    print(f"\n====================================")
    print(f"Procesando familia: {familia_nombre}")
    print(f"Genes únicos: {len(unique_genes)}")
    print(f"====================================")

    combined_fasta = []
    not_found = []

    for gene, agi in unique_genes.items():

        print(f"Descargando {gene} ({agi})...")

        try:

            # Búsqueda robusta priorizando RefSeq y mRNA completos
            search = Entrez.esearch(
                db="nucleotide",
                term=(
                    f"{agi}[Gene] AND "
                    f"Arabidopsis thaliana[Organism] AND "
                    f"biomol_mrna[PROP]"
                ),
                retmax=10
            )

            record = Entrez.read(search)

            if not record["IdList"]:
                print(f"No encontrado: {gene}")
                not_found.append(gene)
                continue

            selected_id = None

            # Priorización robusta
            for seq_id in record["IdList"]:

                summary_handle = Entrez.esummary(db="nucleotide", id=seq_id)
                summary = Entrez.read(summary_handle)

                title = summary[0].get("Title", "")
                accession = summary[0].get("Caption", "")

                # Prioridad máxima: RefSeq curado
                if accession.startswith("NM_"):
                    selected_id = seq_id
                    break

                # Segunda prioridad: RefSeq computacional
                elif accession.startswith("XM_"):
                    selected_id = seq_id

                # Tercera prioridad: mRNA completo
                elif "mRNA" in title and "partial" not in title.lower():
                    if selected_id is None:
                        selected_id = seq_id

            # Fallback
            if selected_id is None:
                selected_id = record["IdList"][0]

            fetch = Entrez.efetch(
                db="nucleotide",
                id=selected_id,
                rettype="fasta",
                retmode="text"
            )

            fasta_data = fetch.read()

            # Guardar FASTA individual
            with open(os.path.join(output_dir, f"{gene}.fasta"), "w") as f:
                f.write(fasta_data)

            combined_fasta.append(fasta_data)

            print(f"{gene} descargado correctamente")

            # Evitar saturar NCBI
            time.sleep(0.4)

        except Exception as e:
            print(f"Error con {gene}: {e}")
            not_found.append(gene)

    # FASTA combinado
    combined_path = os.path.join(
        output_dir,
        f"{familia_nombre}_Arabidopsis_thaliana_COMBINADO.fasta"
    )

    with open(combined_path, "w") as out_fasta:
        for fasta in combined_fasta:
            out_fasta.write(fasta)

    # Log de genes no encontrados
    if not_found:

        log_path = os.path.join(output_dir, "genes_no_encontrados.txt")

        with open(log_path, "w") as log:
            for gene in not_found:
                log.write(gene + "\n")

    print(f"\nFamilia {familia_nombre} completada")
    print(f"FASTA combinado generado:\n{combined_path}")


# =========================================
# HSP40 / DNAJ
# =========================================

HSP40 = {

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
    "DJC30": "AT3G44110"
}

# =========================================
# HSP100
# =========================================

HSP100 = {

    "HSP101": "AT1G74310",
    "HSP98.7": "AT5G15450",
    "HSP93V": "AT5G50920",
    "HSP93III": "AT5G51070",
    "CLPB3": "AT2G25140",
    "CLPP2": "AT5G23140"
}

# =========================================
# HSP90
# =========================================

HSP90 = {

    "HSP90.1": "AT5G52640",
    "HSP90.2": "AT5G56030",
    "HSP90.3": "AT5G56010",
    "HSP90.4": "AT5G56000",
    "HSP90.5": "AT2G04030",
    "HSP90.6": "AT3G07770",
    "HSP90.7": "AT4G24190"
}

# =========================================
# HSP70
# =========================================

HSP70 = {

    "HSP70-1": "AT5G02500",
    "HSP70-2": "AT5G02490",
    "HSP70-3": "AT3G12580",
    "HSP70-4": "AT3G09440",
    "HSP70-5": "AT1G16030",
    "HSP70-6": "AT1G56410",
    "HSP70-7": "AT5G42020",
    "HSP70-8": "AT1G79930",
    "BiP1": "AT5G28540",
    "BiP2": "AT5G42020",
    "BiP3": "AT1G09080"
}

# =========================================
# HSP60
# =========================================

HSP60 = {

    "CPN60A1": "AT2G28000",
    "CPN60A2": "AT1G55490",
    "CPN60B1": "AT1G55490",
    "CPN60B2": "AT3G13470",
    "CPN60B3": "AT5G56500",
    "CPN60B4": "AT1G26230"
}

# =========================================
# sHSPs
# =========================================

sHSP = {

    "HSP17.4": "AT3G46230",
    "HSP17.6A": "AT5G12020",
    "HSP17.6II": "AT5G12030",
    "HSP17.6C": "AT1G53540",
    "HSP17.7": "AT5G12020",
    "HSP18.1": "AT5G59720",
    "HSP18.2": "AT5G59730",
    "HSP21": "AT4G27670",
    "HSP22": "AT4G10250",
    "HSP23.5": "AT5G51440",
    "HSP23.6": "AT4G25200",
    "HSP26.5": "AT1G52560"
}

# =========================================
# EJECUCIÓN AUTOMÁTICA
# =========================================

familias = {
    "HSP40": HSP40,
    "HSP100": HSP100,
    "HSP90": HSP90,
    "HSP70": HSP70,
    "HSP60": HSP60,
    "sHSP": sHSP
}

for nombre, genes in familias.items():
    descargar_familia_hsp(nombre, genes)

print("\n====================================")
print("DESCARGA TOTAL COMPLETADA")
print("====================================")
