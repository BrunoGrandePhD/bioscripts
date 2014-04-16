#!/usr/bin/env python


"""
tabulate_mafs.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This script takes in a number of MAF files and tabulates them
in terms of mutated genes. More specifically, for each gene, it
determine if there is a mutation that belongs to one of a few
categories of mutation types

Known Issues
~~~~~~~~~~~~
- Loads all genes into memory before writing out the results

Example Use
~~~~~~~~~~~

"""


import argparse


MAF_FIELDNAMES = [
    'Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build',
    'Chromosome', 'Start_Position', 'End_Position', 'Strand',
    'Variant_Classification', 'Variant_Type', 'Reference_Allele',
    'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS',
    'dbSNP_Val_Status', 'Tumor_Sample_Barcode',
    'Matched_Norm_Sample_Barcode', 'Match_Norm_Seq_Allele1',
    'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1',
    'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1',
    'Match_Norm_Validation_Allele2', 'Verification_Status',
    'Validation_Status', 'Mutation_Status', 'Sequencing_Phase',
    'Sequence_Source', 'Validation_Method', 'Score', 'BAM_File',
    'Sequencer', 'Tumor_Sample_UUID', 'Matched_Norm_Sample_UUID',
    'Annotation_Transcript', 'Transcript_Position', 'cDNA_Change',
    'Protein_Change', 'effect', 'categ']

# Different categories of mutation types, in decreasing order of severity.
# Names of mutation types are according to VEP classification, stored in
# the Variant_Classification column in the MAF files.
CATEGORY_1 = ['Nonsense_Mutation']

CATEGORY_2 = ['Splice_Site']

CATEGORY_3 = []

CATEGORY_4 = ['Missense_Mutation']

CATEGORY_5 = ['Silent']


def main():
    print
    parser = argparse.ArgumentParser(description='Tabulate a number of MAF ' +
                                     'files where the number of times a gene' +
                                     ' is mutated in various ways is ' +
                                     'counted (once per sample).')
    parser.add_argument('-i', '--input', nargs='+', required=True,
                        type=argparse.FileType('r'),
                        help='Specify a list of MAF files whose mutated ' +
                        'genes will be tabulated.')
    parser.add_argument('-o', '--output', nargs=1, required=True,
                        type=argparse.FileType('w'),
                        help='Specify the output file where to store the ' +
                        'results of the MAF files tabulation.')
    args = parser.parse_args()
    input_mafs = args.input
    genes_total = {}
    gene_template = {
        'category_1': 0,
        'category_2': 0,
        'category_3': 0,
        'category_4': 0,
        'category_5': 0,
    }

    # Loop over every MAF file
    for maf in input_mafs:
        genes_per_file = {}
        for row in maf:
            parsed_row = parse_maf_row(row)
            current_gene = parsed_row['Hugo_Symbol']
            current_variant_class = parsed_row['Variant_Classification']
            # Create gene entry if doesn't exist already
            if parsed_row['Hugo_Symbol'] not in genes_per_file:
                genes_per_file[current_gene] = gene_template.copy()
            # Stratify variant
            if current_variant_class in CATEGORY_1:
                genes_per_file[current_gene]['category_1'] += 1
            if current_variant_class in CATEGORY_2:
                genes_per_file[current_gene]['category_2'] += 1
            if current_variant_class in CATEGORY_3:
                genes_per_file[current_gene]['category_3'] += 1
            if current_variant_class in CATEGORY_4:
                genes_per_file[current_gene]['category_4'] += 1
            if current_variant_class in CATEGORY_5:
                genes_per_file[current_gene]['category_5'] += 1
        # Transfer from genes_per_file to genes_total
        for gene, variants in genes_per_file.items():
            # Create gene entry if doesn't exist already
            if gene not in genes_total:
                genes_total[gene] = gene_template.copy()
            # Stratifying again!
            if variants['category_1'] > 0:
                genes_total[gene]['category_1'] += 1
            if (variants['category_1'] > 0 or variants['category_2'] > 0):
                genes_total[gene]['category_2'] += 1
            if (variants['category_1'] > 0 or variants['category_2'] > 0 or
                    variants['category_3'] > 0):
                genes_total[gene]['category_3'] += 1
            if (variants['category_1'] > 0 or variants['category_2'] > 0 or
                    variants['category_3'] > 0 or variants['category_4'] > 0):
                genes_total[gene]['category_4'] += 1
            if (variants['category_1'] > 0 or variants['category_2'] > 0 or
                    variants['category_3'] > 0 or variants['category_4'] > 0 or
                    variants['category_5'] > 0):
                genes_total[gene]['category_5'] += 1
    print print_tabulation(genes_total)


def parse_maf_row(row):
    """Parse MAF row by creating a dictionary, with a key-value
    pair for each column, according to MAF_FIELDNAMES.
    """
    row_items = row.rstrip('\n').split('\t')
    row_dict = dict(zip(MAF_FIELDNAMES, row_items))
    return row_dict


def print_tabulation(genes_total):
    """Generate string for the number of variants in each category
    for each gene.
    """
    tabulation = ''
    for gene, variants in sorted(genes_total.items(), reverse=True,
                                 key=lambda x:x[1]['category_5']):
        if gene is '':
            continue
        tabulation += gene + '\t'
        tabulation += str(variants['category_1']) + '\t'
        tabulation += str(variants['category_2']) + '\t'
        tabulation += str(variants['category_3']) + '\t'
        tabulation += str(variants['category_4']) + '\t'
        tabulation += str(variants['category_5']) + '\n'
    return tabulation


if __name__ == '__main__':
    main()
