# libraries
from Bio import Entrez
from Bio import SeqIO

Entrez.email = 'hans@orikami.nl'
h = Entrez.einfo()
record = Entrez.read(h)
h = Entrez.esearch(
    db='nucleotide',
#    term='ebola AND complete genome[title]'
    term="NC_002549.1"
)
result = Entrez.read(h)

with Entrez.efetch(
    db="nucleotide", rettype="fasta", retmode="text", id=result['IdList'][0]
) as handle:
    seq_record = SeqIO.read(handle, "fasta")


h.close()
print('finished')