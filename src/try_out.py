# libraries
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "hans@orikami.nl"  # Always tell NCBI who you are
handle = Entrez.einfo()
result = handle.read()
'''
handle = Entrez.esearch(db="nucleotide", term="Cypripedioideae[Orgn] AND matK[Gene]", idtype="acc")
record = Entrez.read(handle)
print(record["Count"])
print(record["IdList"])
'''
handle = Entrez.esearch(db="gene", term='"Homo sapiens"[Organism] AND HBA1[gene]')
record = Entrez.read(handle)
gi_list = record["IdList"]
handle = Entrez.efetch(db="gene", id=gi_list, retmode='xml')
#records = Entrez.parse(handle)
text = handle.read()
#print(text)
records = SeqIO.parse(handle, "seqxml")
for record in records:
    print("%s, length %i, with %i features"
          % (record.name, len(record), len(record.features)))


handle.close()

print('finished')
# load relevant genetic data
# from homo sapiens
# from pan
# from mouse?


# build tree - https://biopython.org/wiki/Phylo

# visualize
