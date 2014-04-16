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
import time
#next three lines are for Python REST API
import httplib2
http = httplib2.Http(".cache")
server = "http://beta.rest.ensembl.org/"

local_database = True

if local_database:
    import cancerGenome
    db = cancerGenome.cancerGenomeDB(
        database_name='lymphoma_meta_hg19',
        database_host='jango.bcgsc.ca',
        database_user='rmorin',
        database_password='rmorin'
    )


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
    parser.add_argument('-s', '--sort', nargs=1, default='5',
                        choices=['1', '2', '3', '4', '5'],
                        help='Specify which category of mutations you want ' +
                        'the genes to be sorted in decreasing order. See ' +
                        'script header for mutation categories.')

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
    sorted_category = 'category_' + args.sort[0]

    # To ensure that the same mutation from a given patient is counted twice
    # key is patient_id (the normal sample name)
    # value is dictionary with key = chr:position:alt_allele and value = 1
    mutations_seen_in_patients = {}

    # Loop over every MAF file
    for maf in input_mafs:
        genes_per_file = {}
        for row in maf:
            parsed_row = parse_maf_row(row)
            mutation_id = (parsed_row['Chromosome'] + ':' +
                           parsed_row['Start_Position'] + ':' +
                           parsed_row['Tumor_Seq_Allele1'])
            patient_id = parsed_row['Matched_Norm_Sample_Barcode']
            if not patient_id in mutations_seen_in_patients:
                mutations_seen_in_patients[patient_id] = {}
            # Check if this mutation has already been seen before in this
            # patient
            if (mutation_id in mutations_seen_in_patients[patient_id]):
                continue
            else:  # If not, add it to the list of mutations seen in this
                   # patient
                mutations_seen_in_patients[patient_id][mutation_id] = 1
            current_gene = parsed_row['Entrez_Gene_Id']
            current_variant_class = parsed_row['Variant_Classification']

            # Create gene entry if doesn't exist already
            if parsed_row['Entrez_Gene_Id'] not in genes_per_file:
                genes_per_file[current_gene] = gene_template.copy()
                genes_per_file[current_gene]['Hugo_Symbol'] = \
                    parsed_row['Hugo_Symbol']
                transcript_length = get_transcript_length(
                    parsed_row['Annotation_Transcript']
                )
                genes_per_file[current_gene]['transcript_length'] = \
                    transcript_length
            else:  # Update the length if a longer transcript is found
                new_transcript_length = get_transcript_length(
                    parsed_row['Annotation_Transcript']
                )
                if (new_transcript_length >
                        genes_per_file[current_gene]['transcript_length']):
                    genes_per_file[current_gene]['transcript_length'] = \
                        new_transcript_length

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
        for gene, results in genes_per_file.items():
            # Create gene entry if doesn't exist already
            if gene not in genes_total:
                genes_total[gene] = gene_template.copy()
                genes_total[gene]['Hugo_Symbol'] = results['Hugo_Symbol']
                genes_total[gene]['transcript_length'] = \
                    results['transcript_length']
            else:
                new_transcript_length = results['transcript_length']
                if (new_transcript_length >
                        genes_total[gene]['transcript_length']):
                    genes_total[gene]['transcript_length'] = \
                        new_transcript_length
            # Stratifying again!
            if results['category_1'] > 0:
                genes_total[gene]['category_1'] += 1
            if (results['category_1'] > 0 or results['category_2'] > 0):
                genes_total[gene]['category_2'] += 1
            if (results['category_1'] > 0 or results['category_2'] > 0 or
                    results['category_3'] > 0):
                genes_total[gene]['category_3'] += 1
            if (results['category_1'] > 0 or results['category_2'] > 0 or
                    results['category_3'] > 0 or results['category_4'] > 0):
                genes_total[gene]['category_4'] += 1
            if (results['category_1'] > 0 or results['category_2'] > 0 or
                    results['category_3'] > 0 or results['category_4'] > 0 or
                    results['category_5'] > 0):
                genes_total[gene]['category_5'] += 1
    args.output[0].write(print_tabulation(genes_total, sorted_category))


def parse_maf_row(row):
    """Parse MAF row by creating a dictionary, with a key-value
    pair for each column, according to MAF_FIELDNAMES.
    """
    row_items = row.rstrip('\n').split('\t')
    row_dict = dict(zip(MAF_FIELDNAMES, row_items))
    return row_dict


def print_tabulation(genes_total, sorted_category):
    """Generate string for the number of variants in each category
    for each gene.
    """
    tabulation = ''
    for gene_id, results in sorted(genes_total.items(), reverse=True,
                                   key=lambda x: x[1][sorted_category]):
        if results['Hugo_Symbol'] is '':
            gene_name = gene_id
        else:
            gene_name = results['Hugo_Symbol']
        transcript_length = float(results['transcript_length'])
        tabulation += gene_name + '\t'
        tabulation += str(int(transcript_length)) + '\t'
        tabulation += str(results['category_1']) + '\t'
        tabulation += str(results['category_2']) + '\t'
        tabulation += str(results['category_3']) + '\t'
        tabulation += str(results['category_4']) + '\t'
        tabulation += str(results['category_5']) + '\t'
        if int(transcript_length) != 0:
            tabulation += str(int(results['category_1'] * 1000000 /
                              transcript_length)) + '\t'
            tabulation += str(int(results['category_2'] * 1000000 /
                              transcript_length)) + '\t'
            tabulation += str(int(results['category_3'] * 1000000 /
                              transcript_length)) + '\t'
            tabulation += str(int(results['category_4'] * 1000000 /
                              transcript_length)) + '\t'
            tabulation += str(int(results['category_5'] * 1000000 /
                              transcript_length)) + '\n'
        else:  #
            tabulation += '\t' + '\t' + '\t' + '\t' + '\n'
    return tabulation


def get_transcript_length(ensembl_transcript):
    if local_database:
        print "getting length from local db"

        trans_obj = cancerGenome.Transcript(
            db.db, ensembl_transcript_id=ensembl_transcript
        )
        try:
            transcript_length = trans_obj.cds_length
        except AttributeError:
            #transcript probably not in database
            transcript_length = 0
    else:
        time.sleep(0.1)  # so we don't get in trouble from Ensembl
        ext = ("sequence/id/" + ensembl_transcript +
               "?content-type=text/plain;type=cds;")
        resp, text_content = http.request(server+ext, method="GET")
        if not resp.status == 200:
            print ("Invalid response for " + ensembl_transcript + ": ",
                   resp.status)
            print ext
            return 0
        print ensembl_transcript + ': ' + str(len(text_content))
        transcript_length = len(text_content)
    return transcript_length


if __name__ == '__main__':
    main()
