# libraries
from Bio import Entrez
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO

aln = AlignIO.read('../data/msa.phy', 'phylip')
print(aln)
calculator = DistanceCalculator('blastn')
dm = calculator.get_distance(aln)
print(dm)

Entrez.email = "hans@orikami.nl"  # Always tell NCBI who you are
handle = Entrez.einfo()
result = handle.read()
handle = Entrez.esearch(db="nucleotide", term="Cypripedioideae[Orgn] AND matK[Gene]", idtype="acc")
record = Entrez.read(handle)
print(record["Count"])
print(record["IdList"])
handle.close()
print('finished')
# load relevant genetic data
# from homo sapiens
# from pan
# from mouse?


# build tree - https://biopython.org/wiki/Phylo

# visualize
