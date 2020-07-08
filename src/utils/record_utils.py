"""
Functions for plotting
"""
from Bio import SeqRecord, SeqFeature


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
        for rec in records:
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