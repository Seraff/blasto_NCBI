#!/usr/bin/env python3
from Bio import SeqIO
from pprint import pprint
import uuid
import csv

from natsort import os_sorted

GENES_FILENAME = 'genes.fasta'
FUNC_ANNOT_FILENAME = 'func_annot.tsv'
STRUCT_ANNOT_FILENAME = 'struct_annot.gff'
TR_RNA_TABLE = 't_r_rnas.tsv'
OUTPUT_PATH = 'final/result.gff'
# OUTPUT_PATH = 'final/LWC14.tbl'

LOCUS_PREFIX = 'LWC14'
DEFAULT_FUNCTION = 'hypothetical protein'


def make_feautre_header(contig_id):
    return f">Feature\t{contig_id}\t\t\t\n"


def make_table_entry(gene_obj):
    result = ''

    coords = gene_obj['coords']

    if gene_obj['source'] != 'RNA':
        if gene_obj['direction'] == '-':
            coords.reverse()
            coords[1] = f">{coords[1]}"

        elif gene_obj['direction'] == '+':
            coords[1] = f">{coords[1]}"
        else:
            raise(Exception(f'Weird direction in {gene_obj}'))

    # dirty workaround
    if gene_obj['gene_name'] == '62d6d129cb58c511855043acab2823dd':
        coords[0] = f"<{coords[0]}"

    coords = [str(e) for e in coords]

    prefix = '\t'.join(coords)
    result += f"{prefix}\tgene\t\t\n"
    result += f"\t\t\tlocus_tag\t{gene_obj['locus_tag']}\n"

    prefix = '\t'.join(coords)
    result += f"{prefix}\t{'CDS' if gene_obj['type'] == 'gene' else gene_obj['type']}\t\t\n"

    gene_func = gene_obj['gene_func'] if gene_obj['gene_func'] else DEFAULT_FUNCTION
    result += f"\t\t\tproduct\t{gene_func}\n"

    return result


def dict_to_gff_description(d):
    return ';'.join([f"{k}={v}" for k, v in d.items()])


def join_gff_list(l):
    return '\t'.join([str(e) for e in l])

def make_gff_entry(gene_obj):
    # {
    # 'contig_id': 'Ctg44_length_309375',
    # 'contig_length': 309375,
    # 'coords': [118692, 120407],
    # 'direction': '+',
    # 'frame': 2,
    # 'gene_func': 'Rieske-like [2Fe-2S] domain containing protein, putative',
    # 'gene_name': '47b2ceaaed138458aded4cc06e4178bc',
    # 'locus_tag': 'LWC14_00001',
    # 'source': 'annotation',
    # 'type': 'gene',
    # 'zoi_coords': [118434, 121751]
    # }

    rows = []

    source = gene_obj['source']
    product = gene_obj['gene_func'] if gene_obj['gene_func'] else DEFAULT_FUNCTION

    if gene_obj['scheme'] == 'annotation':

        # 1. gene (ID=Pelo_1;locus_tag=Pelo_1)
        record = [gene_obj['contig_id'], source, 'gene']
        record.extend(gene_obj['zoi_coords'])
        record.extend(['.', gene_obj['direction'], gene_obj['frame']])

        gene_id = gene_obj['locus_tag']
        description = {'ID': gene_id, 'locus_tag': gene_id}
        record.append(dict_to_gff_description(description))
        rows.append(join_gff_list(record))


        # 2. mRNA
        record = [gene_obj['contig_id'], source, 'mRNA']
        record.extend(gene_obj['zoi_coords'])
        record.extend(['.', gene_obj['direction'], gene_obj['frame']])

        mrna_id = f"{gene_id}_1"
        description = {'ID': mrna_id, 'Parent': gene_id, 'product': product}
        record.append(dict_to_gff_description(description))
        rows.append(join_gff_list(record))


        # 3. exon
        record = [gene_obj['contig_id'], source, 'exon']
        record.extend(gene_obj['coords'])
        record.extend(['.', gene_obj['direction'], gene_obj['frame']])

        exon_id = f"{gene_id}.exon1"
        description = {'ID': exon_id, 'Parent': mrna_id, 'product': product}
        record.append(dict_to_gff_description(description))
        rows.append(join_gff_list(record))

        # 4. CDS
        record = [gene_obj['contig_id'], source, 'CDS']
        record.extend(gene_obj['coords'])
        record.extend(['.', gene_obj['direction'], gene_obj['frame']])

        cds_id = f"cds.{gene_id}"
        description = {'ID': cds_id, 'Parent': mrna_id, 'product': product}
        record.append(dict_to_gff_description(description))
        rows.append(join_gff_list(record))

    elif gene_obj['scheme'] == 'trRNA':
        # 1. gene (ID=Pelo_1;locus_tag=Pelo_1)
        record = [gene_obj['contig_id'], source, 'gene']
        record.extend(gene_obj['coords'])
        record.extend(['.', gene_obj['direction'], gene_obj['frame']])

        gene_id = gene_obj['locus_tag']
        description = {'ID': gene_id, 'locus_tag': gene_id}
        record.append(dict_to_gff_description(description))
        rows.append(join_gff_list(record))

        # 2. t- r- RNA
        record = [gene_obj['contig_id'], source, gene_obj['type']]
        record.extend(gene_obj['coords'])
        record.extend(['.', gene_obj['direction'], gene_obj['frame']])

        mrna_id = f"{gene_id}_1"
        description = {'ID': mrna_id, 'Parent': gene_id, 'product': product}
        record.append(dict_to_gff_description(description))
        rows.append(join_gff_list(record))


    # [print(e + "\n") for e in rows]

    result = '\n'.join(rows)
    result += '\n'

    return result


def parse_struct_annotation_description(raw_field):
    description = raw_field.strip().split(';')
    description = [[e.split('=')[0], e.split('=')[1]] for e in description]
    description = dict(description)
    return description


def main():
    annotated_genes = {}
    struct_annot_data = {}
    gene_cnt = 1

    # get data from structural annotation

    with open(STRUCT_ANNOT_FILENAME) as f:
        for line in f:
            line = line.strip()
            if line == '':
                continue

            splitted = line.split('\t')
            key = (int(splitted[3]), int(splitted[4]))
            if key not in struct_annot_data:
                struct_annot_data[key] = {}
            struct_annot_data[key]['direction'] = splitted[6].strip()

            description = parse_struct_annotation_description(splitted[8].strip())

            struct_annot_data[key]['zoi_start'] = description['zoi_start']
            struct_annot_data[key]['zoi_finish'] = description['zoi_finish']

            struct_annot_data[key]['contig_id'] = splitted[0].strip()
            struct_annot_data[key]['frame'] = splitted[7].strip()


    # Preparing predicted genes, fetch data from functional annotation
    with open(FUNC_ANNOT_FILENAME, newline='') as f:
        for csv_rec in csv.DictReader(f, delimiter='\t'):
            contig_id = csv_rec['qseqid'].split(' ')[1].split('|')[0]

            if contig_id not in annotated_genes:
                annotated_genes[contig_id] = []

            gene_obj = {}

            gene_obj['scheme'] = 'annotation'

            seq_id = csv_rec['qseqid'].split(' ')[0]
            gene_obj['gene_name'] = seq_id
            gene_obj['locus_tag'] = f'{LOCUS_PREFIX}_{gene_cnt:05d}'

            gene_obj['coords'] = csv_rec['qseqid'].split('|')[1].split('-')
            gene_obj['coords'] = [int(e) for e in gene_obj['coords']]

            key_coords = tuple([int(e) for e in gene_obj['coords']])
            if key_coords not in struct_annot_data:
                raise(Exception(f'Cannot find coordinates {key_coords} in structural annotation'))
            else:
                gene_obj['direction'] = struct_annot_data[key_coords]['direction']
                zoi_coords = []
                zoi_coords.append(struct_annot_data[key_coords]['zoi_start'])
                zoi_coords.append(struct_annot_data[key_coords]['zoi_finish'])
                gene_obj['zoi_coords'] = [int(e) for e in zoi_coords]
                gene_obj['contig_id'] = struct_annot_data[key_coords]['contig_id']
                gene_obj['frame'] = struct_annot_data[key_coords]['frame']

            contig_length = int(contig_id.split('_')[-1])
            gene_obj['contig_length'] = contig_length

            if 'no hit found' in csv_rec['sseqdef']:
                gene_obj['gene_func'] = None
            else:
                gene_obj['gene_func'] = csv_rec['sseqdef']

            gene_obj['type'] = 'gene'
            gene_obj['source'] = 'nonstop_annotator'

            annotated_genes[contig_id].append(gene_obj)

            gene_cnt += 1


    # Checking if we have func. annotation for all our annotated genes

    for rec in SeqIO.parse(GENES_FILENAME, 'fasta'):
        contig_id = rec.description.split(' ')[1].split('|')[0]

        if rec.id not in [e['gene_name'] for e in annotated_genes[contig_id]]:
            raise(Exception(f'{rec.id} is not in Kika\'s annotation! Meh...'))


    # preparing RNAs
    with open(TR_RNA_TABLE, newline='') as f:
        for csv_rec in csv.DictReader(f, delimiter='\t'):
            contig_id = csv_rec['contig'].strip()

            if contig_id not in annotated_genes:
                annotated_genes[contig_id] = []

            gene_obj = {}

            gene_obj['scheme'] = 'trRNA'

            gene_obj['contig_id'] = csv_rec['contig']
            gene_obj['type'] = csv_rec['type']
            gene_obj['source'] = csv_rec['source'].lower()

            gene_obj['gene_name'] = uuid.uuid4().hex
            gene_obj['locus_tag'] = f'{LOCUS_PREFIX}_{gene_cnt:05d}'
            gene_obj['coords'] = sorted([int(csv_rec['start']), int(csv_rec['end'])])
            gene_obj['direction'] = csv_rec['strain']
            gene_obj['frame'] = '.'
            gene_obj['gene_func'] = csv_rec['product']

            annotated_genes[contig_id].append(gene_obj)
            gene_cnt += 1

    # Making the actual file

    with open(OUTPUT_PATH, 'w') as out_f:
        contig_ids = os_sorted(list(annotated_genes.keys()))

        for contig_id in contig_ids:
            data = annotated_genes[contig_id]
            # >Feature	Ctg44_length_309375
            # out_f.write(make_feautre_header(contig_id))

            data.sort(key=lambda x: x['coords'][0])

            for gene_obj in data:
                # out_f.write(make_table_entry(gene_obj))
                out_f.write(make_gff_entry(gene_obj))


if __name__ == '__main__':
    main()
