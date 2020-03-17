# libraries
from Bio import Entrez
from Bio import SeqIO

# https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?HostLineage_ss=Homo%20sapiens%20(human),%20taxid:9606&SeqType_s=Nucleotide&Flags_ss=refseq&VirusLineage_ss=Middle%20East%20respiratory%20syndrome-related%20coronavirus%20(MERS-CoV),%20taxid:1335626
# https://www.cdc.gov/flu/about/viruses/types.htm

Entrez.email = 'hans@orikami.nl'
h = Entrez.einfo()
record = Entrez.read(h)
h = Entrez.esearch(
    db='nucleotide',
    term='ebola AND refseq[type]'
#    term="NC_002549.1"
)
result = Entrez.read(h)

with Entrez.efetch(
    db="nucleotide", rettype="fasta", retmode="text", id=result['IdList'][0]
) as handle:
    seq_record = SeqIO.read(handle, "fasta")


h.close()
print('finished')