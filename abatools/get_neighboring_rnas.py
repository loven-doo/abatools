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


def get_neighboring_rnas_cmd():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i',
                        '--in_file',
                        help='Path to the input file with list of Uniprot AC: JSON (.json extension) or text file with ACs splitted with spaces, tabs or enters)',
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
    parser.add_argument('-tsv',
                        '--tsv-out',
                        help='Flag to add tsv output',
                        required=False,
                        action='store_true')    

    cmd_args = parser.parse_args()
    uniprot_ac_list = _get_input_list(path_to_input=cmd_args.in_file)
    result = get_neighboring_rnas(uniprot_ac_list=uniprot_ac_list, window=cmd_args.window)
    with open(cmd_args.out_file, 'w') as output_file:
        json.dump(result, output_file, indent=2)
    if cmd_args.tsv_out:
        store_tsv(result, tsv_path=cmd_args.out_file+".tsv")
                

def store_tsv(result, tsv_path):
    column_names = ('RE_AC', 'contig', 'start', 'end', 'strand', 
                    'neighbor_locus', 'n_type', 'n_product', 'n_start', 'n_end', 'n_strand')
    with open(tsv_path, 'w') as tsv_output_file:
        tsv_output_file.write("\t".join(column_names)+"\n")
        for prot in result:
            prot_info = [prot["AC"], prot["contig EMBL ID"], prot["start"], prot["end"], prot["strand"]]
            for neigh in prot["neighbors"]:
                neigh_info = [neigh["locus"], neigh["type"], neigh["product"], neigh["start"], neigh["end"], neigh["strand"]]
                tsv_output_file.write("\t".join(map(str, prot_info+neigh_info))+"\n")


def _get_input_list(path_to_input):
    with open(path_to_input, 'r') as input_file:
        if os.path.splitext(path_to_input)[1] == '.json':
            return json.load(input_file)
        else:
            return [x.strip() for x in filter(lambda y: y.strip(), input_file.read().split())]


def get_neighboring_rnas(uniprot_ac_list: Iterable, window=WINDOW):
    output_list = []
    for uniprot_ac in uniprot_ac_list:
        print('Preparing AC "%s"...' % uniprot_ac)
        try:
            res = get_protein_neighbors(uniprot_ac=uniprot_ac, window=window)
            output_list.append(res)
        except (TypeError, ValueError):
            print(traceback.format_exc())
        except:
            print(traceback.format_exc())
            break
    return output_list


def get_protein_neighbors(uniprot_ac: str, window: int):
    embl_id = get_embl_id(uniprot_ac=uniprot_ac)
    if embl_id is None:
        raise ValueError('No EMBL record ID found')
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
        "contig EMBL ID": our_contig.id,
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


def rna_feature2record(feature: SeqFeature):
    locus = ''
    start = ''
    end = ''
    strand = ''
    product = ''
    rna_type = feature.type
    locus = feature.qualifiers.get('locus_tag')
    feature_start, feature_end, feature_strand = get_feature_location(feature=feature)
    product_qual = feature.qualifiers.get('product')
    note_qual = feature.qualifiers.get('note')
    if product_qual:
        product = ';'.join(product_qual)
    elif note_qual:
        product = ';'.join(note_qual)
        
    return {
        "locus": locus,
        "type": rna_type,
        "product": product,
        "start": str(feature_start), 
        "end": str(feature_end), 
        "strand": feature_strand
    }


def get_embl_id(uniprot_ac: str):
    embl_id = None
    u = UniProt(verbose=False)
    uniprot_record = u.retrieve(uniprot_ac, frmt="txt")
    for record in SeqIO.parse(StringIO(uniprot_record), "swiss"):
        for xref in record.dbxrefs:
            if 'EMBL' in xref:
                embl_id = xref.split(':')[-1].strip()
    return embl_id


def get_embl_data_str(embl_id, frmt='text'):
    s = ENA(verbose=False)
    embl_data = s.get_data(embl_id, frmt=frmt, expanded=True)

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
        for feature in record.features:
            if get_uniprot_ac(feature) == uniprot_ac:
                our_feature = feature
                our_contig = record
                break
        if our_feature is not None and our_contig is not None:
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
        if "RNA" in feature.type:
            feature_start, feature_end, feature_strand = get_feature_location(feature=feature)
            if feature_end > our_start - window:
                if feature_start < our_end + window:
                    if not (feature_start == our_start and feature_end == our_end):
                        list_of_features.append(rna_feature2record(feature))
    return list_of_features


if __name__ == "__main__":
    get_neighboring_rnas_cmd()
