"""
Gene download functions
"""
import xml.etree.ElementTree as ET
from Bio import Entrez, SeqIO


def get_gene_id(organism, gene_name):
    Entrez.email = "hans@orikami.nl"  # Always tell NCBI who you are

    handle = Entrez.esearch(db="gene", term=f'"{organism}"[Organism] AND {gene_name}[gene]')
    record = Entrez.read(handle)
    gi_list = record["IdList"]
    # return the first hit (probably the best hit)
    return gi_list[0]


def get_gene_record(gene_id):
    Entrez.email = "hans@orikami.nl"  # Always tell NCBI who you are

    handle = Entrez.efetch(db="gene", id=[gene_id], retmode='xml')
    # records = Entrez.parse(handle)
    root = ET.fromstring(handle.read())
    #
    seq_id_gis = list()
    froms = list()
    tos = list()
    strands = list()
    for gen in root.find('Entrezgene').find('Entrezgene_locus').findall('Gene-commentary'):
        if gen.find('Gene-commentary_type').get('value') == 'genomic':
            froms.append([val.text for val in gen.find('Gene-commentary_seqs').iter('Seq-interval_from')])
            tos.append([val.text for val in gen.find('Gene-commentary_seqs').iter('Seq-interval_to')])
            seq_id_gis.append([val.text for val in gen.find('Gene-commentary_seqs').iter('Seq-id_gi')])
            strands.append([val.attrib['value'] for val in gen.find('Gene-commentary_seqs').iter('Na-strand')])

    # strand is given in plus or minus. Rework that into 1 and -1
    strands = [-1 if s[0].lower() == "minus" else 1 for s in strands]

    index = 0
    handle = Entrez.efetch(db="nucleotide",
                           id=seq_id_gis[index][0],
                           rettype="gb",
                           strand=strands[index],
                           seq_start=int(froms[index][0]) + 1,
                           seq_stop=int(tos[index][0]) + 1)
    record = SeqIO.read(handle, "gb")

    return record
