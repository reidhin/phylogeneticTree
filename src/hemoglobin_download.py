# libraries
from Bio import Entrez
import xml.etree.ElementTree as ET
import os
from Bio import SeqIO

Entrez.email = "hans@orikami.nl"  # Always tell NCBI who you are

# check dit
# https://www.ncbi.nlm.nih.gov/gene/?term=(%22Ape%22%5BOrganism%5D)+AND+HBB%5BGene%5D

handle = Entrez.esearch(db="gene", term='"Homo sapiens"[Organism] AND HBA1[gene]')
record = Entrez.read(handle)
gi_list = record["IdList"]
handle = Entrez.efetch(db="gene", id=gi_list, retmode='xml')
#records = Entrez.parse(handle)
root = ET.fromstring(handle.read())
#
seq_id_gis = list()
froms = list()
tos = list()
for gen in root.find('Entrezgene').find('Entrezgene_locus').findall('Gene-commentary'):
    if gen.find('Gene-commentary_type').get('value') == 'genomic':
        froms.append([val.text for val in gen.find('Gene-commentary_seqs').iter('Seq-interval_from')])
        tos.append([val.text for val in gen.find('Gene-commentary_seqs').iter('Seq-interval_to')])
        seq_id_gis.append([val.text for val in gen.find('Gene-commentary_seqs').iter('Seq-id_gi')])

index = 1
handle = Entrez.efetch(db="nucleotide",
                       id=seq_id_gis[index][0],
                       rettype="gb",
                       strand=1,
                       seq_start=int(froms[index][0])+1,
                       seq_stop=int(tos[index][0])+1)
record = SeqIO.read(handle, "gb")

cds_seq = ''
for feat in record.features:
    if feat.type.lower()=='cds':
        cds_seq=feat
        break
# get cds sequence
cds_seq.extract(record.seq).translate()

for i in range(len(record.seq)):
    if i in cds_seq:
        print(record.seq[i].upper(), end="")
    else:
        print(record.seq[i].lower(), end="")

(record.seq[37:132] + record.seq[249:454]+ record.seq[603:732]).translate()

with open(os.path.join("..", "data", "Output.xml"), "w") as text_file:
    text_file.write(ET.tostring(root, encoding='utf8').decode('utf8'))
pass