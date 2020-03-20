# libraries
import os
from Bio import Entrez
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline
from io import StringIO
from Bio import AlignIO


# https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?HostLineage_ss=Homo%20sapiens%20(human),%20taxid:9606&SeqType_s=Nucleotide&Flags_ss=refseq&VirusLineage_ss=Middle%20East%20respiratory%20syndrome-related%20coronavirus%20(MERS-CoV),%20taxid:1335626
# https://www.cdc.gov/flu/about/viruses/types.htm

# please cite us: https://www.ncbi.nlm.nih.gov/pubmed/27899678


def download_data(genomeAccessionNumbers: list):
    Entrez.email = 'hans@orikami.nl'  # Always tell NCBI who you are
    search = " ".join(genomeAccessionNumbers)
    result = Entrez.read(Entrez.esearch(db="nucleotide", term=search, retmode="xml"))
    seq_records = []
    for id in result['IdList']:
        handle = Entrez.efetch(db="nucleotide", id=id, rettype="fasta", retmode="text")
        seq_records.append(SeqIO.read(handle, "fasta"))
    return seq_records

'''
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
'''

if __name__ == '__main__':
    print(os.getcwd())
    accession_numbers = {
        'MERS': 'NC_019843',
        'Corona': 'NC_045512',
#        'H1N1': 'NC_026436'
    }
    seqs = download_data(list(accession_numbers.values()))
    SeqIO.write(seqs, os.path.join('..', 'data', "example.fasta"), "fasta")
    muscle_exe = os.path.join('..', 'muscle.exe')
    muscle_cline = MuscleCommandline(muscle_exe, input=os.path.join('..', 'data', "example.fasta"))
    stdout, stderr = muscle_cline()
    align = AlignIO.read(StringIO(stdout), "fasta")
    print(align)
    print('finished')

