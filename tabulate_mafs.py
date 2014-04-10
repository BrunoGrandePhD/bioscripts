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


def main():
    parser = argparse.ArgumentParser(description='Integrate VCF files from ' +
                                     'an exome and genome.')
    parser.add_argument('-i', '--input', nargs='+', required=True,
                        type=argparse.FileType('r'),
                        help='Specify both VCF files, one from the exome ' +
                        'and the other from the genome, separated by a space.')
    parser.add_argument('-o', '--output', nargs=1, required=True,
                        type=argparse.FileType('w'),
                        help='Specify the output file where to store the' +
                        'results of the MAF files tabulation.')
    args = parser.parse_args()
    input_mafs = args.input
    genes = {}
    gene_template = {}

    # Loop over every MAF file
    for maf in input_mafs:
        pass


if __name__ == '__main__':
    main()
