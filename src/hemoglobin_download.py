# libraries
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import os
from io import StringIO

from Bio import Entrez, SeqIO, SeqRecord, SeqFeature, AlignIO
from Bio.SeqUtils import GC
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline


# check dit
# https://www.ncbi.nlm.nih.gov/gene/?term=(%22Ape%22%5BOrganism%5D)+AND+HBB%5BGene%5D

# https://www.ensembl.org/Homo_sapiens/Gene/Compara_Tree?db=core;g=ENSG00000244734;r=11:5225464-5229395

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


def plot_basics(recs, labels):
    fig, axs = plt.subplots(1, 2)
    axs[0].bar(range(len(recs)), [len(rec) for rec in recs])
    axs[0].set_ylabel('Length')
    axs[1].bar(range(len(recs)), [GC(rec.seq) for rec in recs])
    axs[1].set_ylabel('GC content')
    for ax in axs:
        ax.set_xticks(range(len(recs)))
        ax.set_xticklabels(
            labels,
            rotation=45,
            ha='right'
        )
    plt.tight_layout()
    return fig


def get_feature(record: SeqRecord, feature:str) -> SeqFeature:
    """
    Returns feature from record
    :param record: Bio.SeqRecord record to return feature from
    :param feature: string of feature to return
    :return: Bio.SeqFeature
    """
    feat_to_return = None
    for feat in record.features:
        if feat.type.lower() == feature:
            feat_to_return = feat
            break
    return feat_to_return


def print_record(record):
    # get features
    feat_source = get_feature(record, 'source')
    feat_gene = get_feature(record, 'gene')
    feat_cds = get_feature(record, 'cds')
    # flip if the strands are different for correct printing
    if feat_source.strand == feat_cds.strand:
        feat_cds_for_printing = feat_cds
    else:
        feat_cds_for_printing = feat_cds._flip(len(record))
    aminos = feat_cds.extract(record.seq).translate()
    list_matches = list("  ".join(["|"]*len(aminos)) + "  ")
    list_aminos_sep = list("==".join([a for a in aminos]) + "==")
    line_length = int(100)
    line = int(0)
    while line*line_length < len(record.seq):
        r = range(line*line_length, (line+1)*line_length)
        # print gene
        print(feat_gene.extract(record.seq)[r.start:r.stop])
        # print matches and aminos
        for temp_list in [list_matches, list_aminos_sep]:
            for i in r:
                if i in feat_cds_for_printing:
                    print(temp_list.pop(0), end ="")
                else:
                    print(" ", end="")
            print()
        print()
        line = line + 1


'''
cds_seq = ''
for feat in record.features:
    if feat.type.lower() == 'cds':
        cds_seq = feat
        break
# get cds sequence
cds_seq.extract(record.seq).translate()

for i in range(len(record.seq)):
    if i in cds_seq:
        print(record.seq[i].upper(), end="")
    else:
        print(record.seq[i].lower(), end="")

(record.seq[37:132] + record.seq[249:454] + record.seq[603:732]).translate()

with open(os.path.join("..", "data", "Output.xml"), "w") as text_file:
    text_file.write(ET.tostring(root, encoding='utf8').decode('utf8'))
pass
'''

def align_sequences(input_file: str, output_file: str = "alignment.fasta") -> MultipleSeqAlignment:
    """
    Aligns the sequences using the muscle algorithm
    :param input_file: fasta-file with the input sequences
    :param output_file: save as aligned fasta-file
    :return: MultipleSeqAlignment with the alignment result
    """
    # run   muscle to align all sequences
    # can also be ran online:
    # https://www.ebi.ac.uk/Tools/services/web/toolresult.ebi?jobId=muscle-I20200329-210908-0063-87869209-p2m
    '''
    /nfs/public/ro/es/appbin/linux-x86_64/muscle-3.8.31/muscle -in muscle-I20200329-210908-0063-87869209-p2m.upfile -verbose -log muscle-I20200329-210908-0063-87869209-p2m.output -quiet -fasta -out muscle-I20200329-210908-0063-87869209-p2m.fasta -tree2 muscle-I20200329-210908-0063-87869209-p2m.dnd
    '''

    # specify where the muscle.exe is located
    muscle_exe = os.path.join('..', 'muscle3.8.31_i86linux64')

    # define the command line for muscle
    muscle_cline = MuscleCommandline(muscle_exe, input=os.path.join('..', 'data', input_file))

    # use 2 iterations; when sequences are far apart, the attempt to reach a more finer alignment leads to an error
    muscle_cline.maxiters = 2

    # report the final command line
    print(muscle_cline)

    # execute the command
    stdout, stderr = muscle_cline()

    # save for later faster processing or testing
    with open(os.path.join('..', 'data', output_file), "w") as alignment_file:
        alignment_file.write(stdout)

    # return the aligned sequences
    return AlignIO.read(StringIO(stdout), "fasta")


def view_alignment(aln, labels):
    #Todo: letters that are not the same should have another color
    """matplotlib sequence alignment view"""

    # make sequence and id lists from the aln object
    mat = [
        [0 if nucleotide == '-' else 1 for nucleotide in rec.seq] for rec in aln
    ]
    cmap = ListedColormap(['w', 'r'])
    fig, ax = plt.subplots(1, 1)
    ax.matshow(mat, cmap=cmap)
    ax.set_aspect(ax.get_xlim()[1] / ax.get_ylim()[0] / 3)
    ax.set_yticks(range(len(aln)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.tick_params(axis="x", bottom=True, top=False, labelbottom=True, labeltop=False)
    plt.tight_layout()
    return fig


if __name__ == '__main__':
    organism_list = [
        "Homo sapiens",
        "Pan troglodytes"
    ]

    gene_name = 'HBB'

    '''
    # get gene ids
    gene_ids = [get_gene_id(organism, gene_name) for organism in organism_list]

    # get gene records
    gene_records = [get_gene_record(gene_id) for gene_id in gene_ids]

    # write them in file for later upload
    SeqIO.write(gene_records, os.path.join('..', 'data', "trog_downloads.gb"), "genbank")
    '''

    gene_records = list(SeqIO.parse(os.path.join('..', 'data', "trog_downloads.gb"), "genbank"))

    # print one gene record
    print_record(gene_records[0])

    fig = plot_basics(gene_records, organism_list)
    #fig.savefig(os.path.join('..', 'figures', 'basic.png'))

    # save as fasta-file for alignment
    SeqIO.write(gene_records, os.path.join('..', 'data', "trog_before_alignment.fasta"), "fasta")

    # align the sequences
    align = align_sequences("trog_before_alignment.fasta", output_file="trog_after_alignment.fasta")
#    align = AlignIO.read(os.path.join("..", "data", "alignment.fasta"), "fasta")
    print(align)
    fig = view_alignment(align, organism_list)
    fig.savefig(os.path.join('..', 'figures', 'trog_alignment.png'), bbox_inches='tight', dpi=300)


    print('finished!')
    # download the sequences
    # recs = download_data(list(accession_numbers.keys()))
