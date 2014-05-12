#!/usr/bin/env python

"""
bioscripts.py
~~~~~~~~~~~~~
This Python script serves as a module for other scripts in Bioscripts.
For instance, standard parsers are stored here.
"""


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


def main():
    pass


def parse_maf(row):
    """Parse MAF row by creating a dictionary, with a key-value
    pair for each column, according to MAF_FIELDNAMES.
    """
    row_items = row.rstrip('\n').split('\t')
    row_dict = dict(zip(MAF_FIELDNAMES, row_items))
    row_dict.update({
        'chromosome': row_dict['Chromosome'],
        'position': int(row_dict['Start_Position']),
        'reference_allele': row_dict['Reference_Allele'],
        'alternate_allele': row_dict['Tumor_Seq_Allele1'],
        'row': row
    })
    return row_dict


def parse_vcf(row, format=None):
    columns = row.strip().split('\t')
    row_dict = {
        'chromosome': columns[0],
        'position': int(columns[1]),
        'id': columns[2],
        'reference_allele': columns[3],
        'alternate_allele': columns[4],
        'quality': columns[5],
        'filter': columns[6],
        'info': columns[7],
        'row': row
    }
    if format is None:
        return row_dict
    else:
        if format is 'museq':
            return parse_vcf_museq(row_dict)
        elif format is 'strelka':
            return parse_vcf_strelka(row_dict)
        else:
            raise Exception('Invalid format for parse_vcf.')


def parse_vcf_museq(row_dict):
    info_items = row_dict['info'].split(';')
    row_dict.update({
        'probability': float(info_items[0][3:]),
        'tumour_number_ref': int(info_items[1][3:]),
        'tumour_number_alt': int(info_items[2][3:]),
        'normal_number_ref': int(info_items[3][3:]),
        'normal_number_alt': int(info_items[4][3:])
    })
    return row_dict


def parse_vcf_strelka(row_dict):
    columns = row_dict['row'].split('\t')
    row_dict = {
        'format': columns[8],
        'normal': columns[9],
        'tumour': columns[10]
    }

    info_items = row_dict['info'].split(';')
    format_items = row_dict['format'].split(':')
    normal_items = row_dict['normal'].split(':')
    tumour_items = row_dict['tumour'].split(':')
    normal_dict = dict(zip(format_items, normal_items))
    tumour_dict = dict(zip(format_items, tumour_items))

    tumour_number_ref = tumour_dict[
        row_dict['reference_allele'] + 'U'
    ].split(',')[1]
    normal_number_ref = normal_dict[
        row_dict['reference_allele'] + 'U'
    ].split(',')[1]

    # Sometimes, the alternative allele isn't there (it's a . instead)
    if row_dict['alternate_allele'] is not '.':

        # Sometimes, there is more than one alternate allele
        if len(row_dict['alternate_allele']) == 1:
            tumour_number_alt = tumour_dict[
                row_dict['alternate_allele'] + 'U'
            ].split(',')[1]
            normal_number_alt = normal_dict[
                row_dict['alternate_allele'] + 'U'
            ].split(',')[1]
        else:  # There is more than one alternate allele
            alt_bases = row_dict['alternate_allele'].split(',')
            count = {}
            for base in alt_bases:
                number = tumour_dict[
                    base + 'U'
                ].split(',')[1]
                count[base] = number
            best_alt_allele = max(count.items(), key=lambda x: x[1])[0]
            tumour_number_alt = tumour_dict[
                best_alt_allele + 'U'
            ].split(',')[1]
            normal_number_alt = normal_dict[
                best_alt_allele + 'U'
            ].split(',')[1]

    else:
        tumour_number_alt = 0
        normal_number_alt = 0

    row_dict.update({
        'probability': float(info_items[1][4:]),
        'tumour_number_ref': int(tumour_number_ref),
        'tumour_number_alt': int(tumour_number_alt),
        'normal_number_ref': int(normal_number_ref),
        'normal_number_alt': int(normal_number_alt)
    })

    return row_dict


if __name__ == '__main__':
    main()
