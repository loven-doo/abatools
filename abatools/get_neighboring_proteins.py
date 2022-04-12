import os
import zlib
import json
import argparse
import traceback
from io import StringIO
from collections.abc import Iterable

from bioservices import ENA
from bioservices.uniprot import UniProt
from Bio import SeqIO, SeqFeature, SeqRecord


WINDOW = 10000


def get_neighboring_proteins_cmd():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i',
                        '--in_file',
                        help='Path to an input file with list of Uniprot AC',
                        required=True,)
    parser.add_argument('-o',
                        '--out_file',
                        help='Path to JSON output',
                        required=True,)
    parser.add_argument('-w',
                        '--window',
                        type=int,                        
                        help='Search area for neighboring genes: 2*window + gene_length; default 10000',
                        required=False,
                        default=WINDOW,)

    cmd_args = parser.parse_args()
    uniprot_ac_list = _get_input_list(path_to_input=cmd_args.in_file)
    with open(cmd_args.out_file, 'w') as output_file:
        json.dump(get_neighboring_proteins(uniprot_ac_list=uniprot_ac_list,
                                           window=cmd_args.window),
                  output_file,
                  indent=2)


def _get_input_list(path_to_input):
    with open(path_to_input, 'r') as input_file:
        if os.path.splitext(path_to_input)[1] == '.json':
            return json.load(input_file)
        else:
            return [x.strip() for x in filter(lambda y: y.strip(), input_file.read().split())]


def get_neighboring_proteins(uniprot_ac_list: Iterable, window=WINDOW):
    output_list = []
    for uniprot_ac in uniprot_ac_list:
        print('Preparing AC "%s"...' % uniprot_ac)
        try:
            res = get_protein_neighbors(uniprot_ac=uniprot_ac, window=window)
            output_list.append(res)
        except (TypeError, ValueError):
            print(traceback.format_exc())
    return output_list


def get_protein_neighbors(uniprot_ac: str, window: int):
    embl_id = get_embl_id(uniprot_ac=uniprot_ac)
    print('Found EMBL ID: %s' % embl_id)
    embl_data_str = get_embl_data_str(embl_id=embl_id)
    print('Downloaded EMBL data for ID %s' % embl_id)
    our_feature, our_contig = find_initial_records(uniprot_ac, embl_data_str, embl_id)
    if our_feature is None or our_contig is None:
        raise ValueError('No initial protein in EMBL data')
    protein_start, protein_end, protein_strand = get_feature_location(our_feature)
    protein_neighbors = get_closest_features(window=window,
                                             our_feature=our_feature,
                                             our_contig=our_contig)
                                             
    return {
        "AC": uniprot_ac,
        "contig EMBL ID": embl_id,
        "start": protein_start,
        "end": protein_end,
        "strand": protein_strand,
        "neighbors": protein_neighbors,
    }
    

def get_uniprot_ac(feature: SeqFeature):
    db_xref = feature.qualifiers.get('db_xref')
    if db_xref:
        for xref in db_xref:
            if 'UniProt' in xref:
                ac = xref.split(':')[-1]
                return ac


def feature_to_record(feature: SeqFeature):
    ac = ''
    start = ''
    end = ''
    strand = ''
    product = ''
    ac = get_uniprot_ac(feature)
    feature_start, feature_end, feature_strand = get_feature_location(feature=feature)
    product_qual = feature.qualifiers.get('product')
    if product_qual:
        product = ';'.join(product_qual)
        
    return {
        "AC": ac, 
        "product": product, 
        "start": str(feature_start), 
        "end": str(feature_end), 
        "strand": feature_strand
    }


def get_embl_id(uniprot_ac: str):
    u = UniProt(verbose=False)
    uniprot_record = u.retrieve(uniprot_ac, frmt="txt")
    for record in SeqIO.parse(StringIO(uniprot_record), "swiss"):
        for xref in record.dbxrefs:
            if 'EMBL' in xref:
                embl_id = xref.split(':')[-1].strip()
    return embl_id


def get_embl_data_str(embl_id):
    s = ENA(verbose=False)
    embl_data = s.get_data(embl_id, frmt='text', expanded=True)

    if isinstance(embl_data, bytes):
        return zlib.decompress(embl_data, 15+32).decode('utf-8')
    elif isinstance(embl_data, str):
        return embl_data
    else:
        raise TypeError('Wrong type of EMBL data: "%s"' % type(embl_data))


def find_initial_records(uniprot_ac: str, embl_data_str: str, embl_id: str):
    our_feature = None
    our_contig = None
    for record in SeqIO.parse(StringIO(embl_data_str), "embl"):
        if embl_id in record.id:
            our_contig = record
            break

    if our_contig is not None:
        for feature in our_contig.features:
            if get_uniprot_ac(feature) == uniprot_ac:
                our_feature = feature
                break

    return our_feature, our_contig


def get_feature_location(feature: SeqFeature):
    feature_start = int(feature.location.start)
    feature_end = int(feature.location.end)
    feature_strand = str(feature.strand)
    return feature_start, feature_end, feature_strand


def get_closest_features(window, our_feature: SeqFeature, our_contig: SeqRecord):
    our_start, our_end, our_strand = get_feature_location(our_feature)
    list_of_features = []
    for feature in our_contig.features:
        if feature.type == 'CDS':
            feature_start, feature_end, feature_strand = get_feature_location(feature=feature)
            if feature_end > our_start - window:
                if feature_start < our_end + window:
                    if not (feature_start == our_start and feature_end == our_end):
                        list_of_features.append(feature_to_record(feature))
    return list_of_features


if __name__ == "__main__":
    get_neighboring_proteins_cmd()

