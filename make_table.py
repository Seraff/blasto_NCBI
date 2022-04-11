#!/usr/bin/env python3
from Bio import SeqIO
from pprint import pprint
import uuid

GENES_FILENAME = 'genes.fasta'
FUNC_ANNOT_FILENAME = 'func_annot.tsv'
RNA_ANNOT_FILENAME = 't_r_rnas.gff'
LOCUS_PREFIX = 'LWC14'
DEFAULT_FUNCTION = 'hypothetical protein'


def make_feautre_header(contig_id):
    return f">Feature\t{contig_id}\t\t\t\n"


def make_table_entry(gene_obj):
    result = ''

    prefix = '\t'.join(gene_obj['coords'])
    result += f"{prefix}\tgene\t\t\n"
    result += f"\t\t\tlocus_tag\t{gene_obj['locus_tag']}\n"

    prefix = '\t'.join(gene_obj['coords'])
    result += f"{prefix}\t{'CDS' if gene_obj['type'] == 'gene' else gene_obj['type']}\t\t\n"

    gene_func = gene_obj['gene_func'] if gene_obj['gene_func'] else DEFAULT_FUNCTION
    result += f"\t\t\tproduct\t{gene_func}\n"

    return result


def main():
    annot_gff = {}
    rna_gff = {}
    gene_cnt = 1

    # Preparing predicted genes

    with open(FUNC_ANNOT_FILENAME) as f:
        for line in f:
            line = line.strip()
            if line == '' or line.startswith('qseqid'):
                continue

            splitted = line.split('\t')

            contig_id = splitted[0].split(' ')[1].split('|')[0]
            if contig_id not in annot_gff:
                annot_gff[contig_id] = []

            gene_obj = {}

            seq_id = splitted[0].split(' ')[0]
            gene_obj['gene_name'] = seq_id
            gene_obj['locus_tag'] = f'{LOCUS_PREFIX}_{gene_cnt:05d}'

            gene_obj['coords'] = splitted[0].split('|')[1].split('-')

            if 'no hit found' in splitted[-1]:
                gene_obj['gene_func'] = None
            else:
                func_data = [e.strip() for e in splitted[3].split('|')]
                func_data = dict([e.split('=') for e in func_data if e != ''])

                gene_obj['gene_func'] = func_data['gene_product']

            gene_obj['type'] = 'gene'

            annot_gff[contig_id].append(gene_obj)
            gene_cnt += 1

    # preparing RNAs

    with open(RNA_ANNOT_FILENAME) as f:
        cnt = 1

        for line in f:
            line = line.strip()
            splitted = line.split('\t')

            if line == '' or splitted[2] == 'gene':
                continue

            contig_id = splitted[0]
            if contig_id not in annot_gff:
                annot_gff[contig_id] = []

            gene_obj = {}

            gene_obj['gene_name'] = uuid.uuid4().hex
            gene_obj['locus_tag'] = f'{LOCUS_PREFIX}_{gene_cnt:05d}'
            gene_obj['coords'] = [splitted[3], splitted[4]]
            gene_obj['gene_func'] = splitted[8].split(';')[-1].split('=')[-1]
            gene_obj['type'] = splitted[2]

            annot_gff[contig_id].append(gene_obj)
            gene_cnt += 1

    # Checking if everythin is ok

    for rec in SeqIO.parse(GENES_FILENAME, 'fasta'):
        contig_id = rec.description.split(' ')[1].split('|')[0]

        if rec.id not in [e['gene_name'] for e in annot_gff[contig_id]]:
            raise(Exception(f'{rec.id} is not in Kika\'s annotation! Meh...'))

    # Making the actual file

    with open('Tb_func_annot.tbl', 'w') as out_f:
        for contig_id, data in annot_gff.items():
            out_f.write(make_feautre_header(contig_id))
            for gene_obj in data:
                out_f.write(make_table_entry(gene_obj))


if __name__ == '__main__':
    main()
