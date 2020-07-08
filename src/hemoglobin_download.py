# libraries
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import os
import pandas as pd
import seaborn as sns
from io import StringIO

from Bio import Entrez, SeqIO, SeqRecord, SeqFeature, AlignIO, Phylo, Align
from Bio.SeqUtils import GC
from Bio.Align import MultipleSeqAlignment, substitution_matrices
from Bio.Align.Applications import MuscleCommandline

from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


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


def print_record(record, line_length=int(100)):
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
    line = int(0)
    out = ""
    while line*line_length < len(record.seq):
        r = range(line*line_length, (line+1)*line_length)
        # print gene
        out = out + feat_gene.extract(record.seq)[r.start:r.stop] + "\n"
        # print(feat_gene.extract(record.seq)[r.start:r.stop])
        # print matches and aminos
        for temp_list in [list_matches, list_aminos_sep]:
            for i in r:
                if i in feat_cds_for_printing:
                    out = out + temp_list.pop(0)
                    # print(temp_list.pop(0), end ="")
                else:
                    out = out + " "
                    # print(" ", end="")
            out = out + "\n"
            # print()
        out = out + "\n"
        # print()
        line = line + 1
    return out

def align_sequences(input_file: str, output_file: str = "virus_alignment.fasta") -> MultipleSeqAlignment:
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


def plot_phylo_tree(align: MultipleSeqAlignment, accession_numbers: dict, title=""):
    """
    Plots a phylogenetic tree
    :param align: MultipleSeqAlignment with the alignment result to be plotted
    :param accession_numbers: dict of accession numbers and their translation to human-understandable names
    :param title: string of title
    :return: figure-handle of the plotted phylogenetic tree
    """
    # calculate distance - https://biopython.org/wiki/Phylo
    calculator = DistanceCalculator('trans')
    dm = calculator.get_distance(align)

    # construct a tree
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)

    # order terminals
    tree.ladderize()

    # remove the names for the non-terminals for better visual appeal
    for non_terminal in tree.get_nonterminals():
        non_terminal.name = ''

    # change accession numbers into human more understandable names
    for terminal in tree.get_terminals():
        terminal.name = accession_numbers[terminal.name]

    print(Phylo.draw_ascii(tree))

    # plot the tree
    fig, ax = plt.subplots(1, 1)
    # draw the resulting tree
    Phylo.draw(tree, show_confidence=False, axes=ax, do_show=False)
    ax.set_xlim(right=1.3*ax.get_xlim()[1])
    ax.set_title(title)
    return fig


def create_records(records=list, feature='cds') -> list:
    """
    :param records:
    :param feature:
    :return:
    """
    out_records = list()
    if feature == 'cds':
        # get cds sequence
        for rec in records:
            feat = get_feature(rec, 'cds')
            out_records.append(
                SeqRecord.SeqRecord(
                    seq=feat.extract(rec.seq),
                    id=rec.id,
                    description=rec.description
                )
            )
    elif feature == 'amino':
        for rec in gene_records:
            feat = get_feature(rec, 'cds')
            out_records.append(
                SeqRecord.SeqRecord(
                    seq=feat.extract(rec.seq).translate(to_stop=True),
                    id=rec.id,
                    description=rec.description
                )
            )
    elif feature == 'cds_and_introns':
        for rec in records:
            feat = get_feature(rec, 'cds')
            feat = SeqFeature.SeqFeature(
                SeqFeature.FeatureLocation(
                    feat.location.start,
                    feat.location.end,
                    feat.strand
                )
            )
            out_records.append(
                SeqRecord.SeqRecord(
                    seq=feat.extract(rec.seq),
                    id=rec.id,
                    description=rec.description
                )
            )
    else:
        raise ValueError("feature should be cds, amino, or cds_and_introns")
    return out_records


def get_distance_dataframe(align: Align.MultipleSeqAlignment) -> pd.DataFrame:
    """
    :param align:
    :return:
    """
    # calculate distance - https://biopython.org/wiki/Phylo
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(align)

    df = pd.DataFrame(
        dm.matrix,
        index=[trans_dict[ac] for ac in dm.names],
        columns=[trans_dict[ac] for ac in dm.names]
    )

    return df


if __name__ == '__main__':
    organism_list = [
        "Homo sapiens",
        "Pan troglodytes",
        "Gorilla gorilla",
        "Pongo abelii",
        "Pan paniscus"
#        "Hylobates moloch"
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

    # translate dict
    trans_dict = dict(
        zip(
            [gr.id for gr in gene_records], organism_list
        )
    )

    # print one gene record
    fig = plt.figure(figsize=(8, 8.5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.text(0.05, 0.95, print_record(gene_records[0], int(100)),
            horizontalalignment='left',
            verticalalignment='top',
            family='monospace',
            fontsize=8,
            transform=ax.transAxes)
    ax.set_axis_off()
    fig.savefig(os.path.join("..", "figures", "human_gene_DNA.png"))

    # get amino-acid sequence
    amino_records = create_records(gene_records, 'amino')
    #fig = plot_basics(amino_records, organism_list)

    # save as fasta-file for alignment
    SeqIO.write(amino_records, os.path.join('..', 'data', "trog_before_alignment.fasta"), "fasta")

    # align the sequences
    align_amino = align_sequences("trog_before_alignment.fasta", output_file="trog_after_alignment.fasta")
    print(align_amino)

    aligner = Align.PairwiseAligner()
    alignments = aligner.align(amino_records[0].seq, amino_records[1].seq)
    for alignment in sorted(alignments):
        print("Score = %.1f:" % alignment.score)
        print(alignment)

    # get cds records
    cds_records = create_records(gene_records, 'cds')
    #fig = plot_basics(cds_records, organism_list)

    # save as fasta-file for alignment
    SeqIO.write(cds_records, os.path.join('..', 'data', "trog_before_alignment.fasta"), "fasta")

    # align the sequences
    align_cds = align_sequences("trog_before_alignment.fasta", output_file="trog_after_alignment.fasta")
    print(align_cds)

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    alignments = aligner.align(cds_records[0].seq, cds_records[1].seq)
    for alignment in sorted(alignments):
        print("Score = %.1f:" % alignment.score)
        print(alignment)

    # get introns sequence
    full_records = create_records(gene_records, 'cds_and_introns')
    #fig = plot_basics(full_records, organism_list)

    # save as fasta-file for alignment
    SeqIO.write(full_records, os.path.join('..', 'data', "trog_before_alignment.fasta"), "fasta")

    # align the sequences
    align_full = align_sequences("trog_before_alignment.fasta", output_file="trog_after_alignment.fasta")
    print(align_full)

    for al, name in zip([align_amino, align_cds, align_full], ["amino-acids", "exons", "exons and introns"]):
        df = get_distance_dataframe(al)
        plt.figure()
        sns.heatmap(
            df*100,
            fmt='3.2f',
            annot=True,
            linewidths=0.5,
            cmap=sns.light_palette("navy"),
            cbar=False,
            square=True
        )
        plt.title(f'Percent difference based on {name}')
        plt.tight_layout()
        plt.savefig(os.path.join('..', 'figures', f'distance_{name}_{gene_name}.png'))

    # plot the resulting tree
    fig = plot_phylo_tree(align_cds, trans_dict, title="Phylogenetic tree based on the CDS of the HBB gene")
    fig.savefig(os.path.join('..', 'figures', f'tree_{gene_name}.png'))

    df_cds = get_distance_dataframe(align_cds)
    df_full = get_distance_dataframe(align_full)
    df_cds = df_cds.fillna(0) + df_cds.transpose().fillna(0)
    df_full = df_full.fillna(0) + df_full.transpose().fillna(0)
    df_intron = (df_full * align_full.get_alignment_length() - df_cds * align_cds.get_alignment_length()) / \
                (align_full.get_alignment_length() - align_cds.get_alignment_length())
    sns.heatmap(
        df_intron*100,
        fmt='3.2f',
        annot=True,
        linewidths=0.5,
        cmap=sns.light_palette("navy"),
        cbar=False,
        square=True
    )
    plt.title('Percent difference')
    plt.tight_layout()
    plt.show()

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # save as fasta-file for alignment
    out_records = list()
    for i in range(2):
        cds = cds_records[i]
        cds.id = f"cds_{i}"
        full = full_records[i]
        full.id = f"full_{i}"
        out_records = out_records + [cds, full]

    SeqIO.write(
        out_records,
        os.path.join('..', 'data', "trog_before_alignment.fasta"),
        "fasta"
    )

    # align the sequences
    align_print = align_sequences("trog_before_alignment.fasta", output_file="trog_after_alignment.fasta")
    print(align_print)

    # todo: check if nucleotide is part of the cds by position
    d = dict(zip([a.id for a in align_print], [i for i in range(len(align_print))]))
    comb_0 = "".join([f.lower() if c == '-' else f for f, c in zip(align_print[d['full_0']], align_print[d['cds_0']])])
    comb_1 = "".join([f.lower() if c == '-' else f for f, c in zip(align_print[d['full_1']], align_print[d['cds_1']])])
    match = "".join(["|" if c1 == c2 else " " for c1, c2 in zip(comb_0, comb_1)])

    line_length = 100
    line = int(0)
    out = ""
    while line * line_length < align_print.get_alignment_length():
        r = range(line * line_length, (line + 1) * line_length)
        # print first sequence
        out = out + comb_0[r.start:r.stop] + "\n"
        out = out + match[r.start:r.stop] + "\n"
        out = out + comb_1[r.start:r.stop] + "\n"
        out = out + "\n"
        # print()
        line = line + 1

    # print one gene record
    fig = plt.figure(figsize=(8, 8.5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.text(0.05, 0.95, out,
            horizontalalignment='left',
            verticalalignment='top',
            family='monospace',
            fontsize=8,
            transform=ax.transAxes)
    ax.set_axis_off()
    fig.savefig(os.path.join("..", "figures", "human_to_chimp.png"))
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    print('finished!')
