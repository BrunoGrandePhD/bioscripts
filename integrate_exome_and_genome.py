#!/usr/bin/env python


"""
integrate_exome_and_genome.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This script takes as input two MutationSeq 4 VCF files originating from the
exome and genome of the same tumour-normal pair. It uses some heuristics to
stratify SNV calls into three categories:
- Good candidates (good_calls)
- Poor candidates that warrant validation (poor_calls)
- Rejected candidates (bad_calls)

Known Issues
~~~~~~~~~~~~
- Comment lines aren't dealt with.
- Assumes that the positions in the rows of the VCF files match up, i.e., the
  same number of lines in both VCF files.

Example Use
~~~~~~~~~~~
python integrate_exome_and_genome.py -i exome.vcf genome.vcf \
    -o good_candidates.vcf
"""


import argparse
import os


# What is the fractional threshold for deciding that there is enough
# presence of a "SNV" in the germline reads that it is probably not somatic?
# Requirement: Positive float >= 1
GERMLINE_THRESHOLD = 0.1

# What is the fractional threshold for the number of reads supporting a SNV
# in both the exome and genome for it to be considered a good candidate?
# Requirement: Positive float >= 1
PARTIAL_SUPPORT_THRESHOLD = 0.1

# If a SNV isn't present in either the exome or genome, how much coverage
# in the one that doesn't support the SNV would be required to label a
# candidate as poor.
# Requirement: Positive integer
SUPPORT_AGAINST_MINIMUM_COVERAGE = 10


def main():
    parser = argparse.ArgumentParser(description='Integrate VCF files from' +
                                     'an exome and genome.')
    parser.add_argument('-i', '--input', nargs=2, type=argparse.FileType('r'),
                        required=True,
                        help='Specify both VCF files, one from the exome ' +
                        'and the other from the genome, separated by a space.')
    parser.add_argument('-o', '--output', nargs=1, type=argparse.FileType('w'),
                        required=True,
                        help='Specify the output VCF file for the good' +
                        'candidates, and the files for the other two ' +
                        'categories will be created automatically.')
    args = parser.parse_args()
    input_file_1, input_file_2 = args.input
    output_good_calls = args.output[0]
    output_poor_calls, output_bad_calls = \
        create_output_files(output_good_calls)

    for row_dict_1 in parse_museq_vcf(input_file_1):
        row_dict_2 = parse_museq_vcf(input_file_2).next()
        if row_dict_1['museq_PR'] > row_dict_1['museq_PR']:
            best_call = row_dict_1
        else:
            best_call = row_dict_2
        # Germline Test (for bad_calls)
        if test_germline(row_dict_1, row_dict_2):
            # Write row to bad_calls
            output_bad_calls.write(create_vcf_row(best_call))
            continue
        # Partial Support Test (for good_calls)
        if test_partial_support(row_dict_1, row_dict_2):
            # Write row to good_calls
            output_good_calls.write(create_vcf_row(best_call))
            continue
        # Support Against Test (for poor_calls)
        if test_support_against(row_dict_1, row_dict_2):
            # Write row to poor_calls
            output_poor_calls.write(create_vcf_row(best_call))
            continue
        # If the rows reach this point, they're sorted to good_calls
        output_good_calls.write(create_vcf_row(best_call))

    input_file_1.close()
    input_file_2.close()
    output_good_calls.close()
    output_poor_calls.close()
    output_bad_calls.close()


def create_output_files(output_good_calls):
    original_name = output_good_calls.name
    root, ext = os.path.splitext(original_name)
    output_poor_calls_name = root + '.poor_calls' + ext
    output_poor_calls = open(output_poor_calls_name, 'w')
    output_bad_calls_name = root + '.bad_calls' + ext
    output_bad_calls = open(output_bad_calls_name, 'w')
    return output_poor_calls, output_bad_calls


def parse_museq_vcf(opened_vcf_file):
    # Sample VCF row
    # 1  14057566    .   G   A   16.20   INDL    INFO
    for line in opened_vcf_file:
        if line.strip()[0] is '#':
            continue
        columns = line.strip().split('\t')
        # Sample MutationSeq info column
        #PR=0.98;TR=174;TA=27;NR=152;NA=0;TC=CGG;NI=2;ND=1
        info_items = columns[7].split(';')
        row_dict = {
            'chromosome': columns[0],
            'position': int(columns[1]),
            'id': columns[2],
            'reference_allele': columns[3],
            'alternate_allele': columns[4],
            'quality': float(columns[5]),
            'filter': columns[6],
            'info': columns[7],
            'museq_PR': float(info_items[0][3:]),
            'museq_TR': int(info_items[1][3:]),
            'museq_TA': int(info_items[2][3:]),
            'museq_NR': int(info_items[3][3:]),
            'museq_NA': int(info_items[4][3:]),
            'museq_TC': info_items[5][3:],
            'museq_NI': int(info_items[6][3:]),
            'museq_ND': int(info_items[7][3:])
        }
        yield row_dict


def create_vcf_row(row_dict):
    row_list = [row_dict['chromosome'], str(row_dict['position']),
                row_dict['id'], row_dict['reference_allele'],
                row_dict['alternate_allele'], str(row_dict['quality']),
                row_dict['filter'], row_dict['info']]
    return '\t'.join(row_list) + '\n'


def test_germline(row_dict_1, row_dict_2):
    total_germline_1 = row_dict_1['museq_NR'] + row_dict_1['museq_NA']
    total_germline_2 = row_dict_2['museq_NR'] + row_dict_2['museq_NA']
    germline_allele_ratio_1 = float(row_dict_1['museq_NA'])/total_germline_1
    germline_allele_ratio_2 = float(row_dict_2['museq_NA'])/total_germline_2
    return (germline_allele_ratio_1 > GERMLINE_THRESHOLD or
            germline_allele_ratio_2 > GERMLINE_THRESHOLD)


def test_partial_support(row_dict_1, row_dict_2):
    total_tumour_1 = row_dict_1['museq_TR'] + row_dict_1['museq_TA']
    total_tumour_2 = row_dict_2['museq_TR'] + row_dict_2['museq_TA']
    tumour_allele_ratio_1 = float(row_dict_1['museq_TA'])/total_tumour_1
    tumour_allele_ratio_2 = float(row_dict_2['museq_TA'])/total_tumour_2
    if row_dict_1['museq_PR'] > row_dict_2['museq_PR']:
        # File 2 is the one that's less supported
        return tumour_allele_ratio_2 > PARTIAL_SUPPORT_THRESHOLD
    else:
        # File 1 is the one that's less supported
        return tumour_allele_ratio_1 > PARTIAL_SUPPORT_THRESHOLD


def test_support_against(row_dict_1, row_dict_2):
    total_tumour_1 = row_dict_1['museq_TR'] + row_dict_1['museq_TA']
    total_tumour_2 = row_dict_2['museq_TR'] + row_dict_2['museq_TA']
    if row_dict_1['museq_PR'] > row_dict_2['museq_PR']:
        # File 2 is the one that's less supported
        return (total_tumour_2 > SUPPORT_AGAINST_MINIMUM_COVERAGE and
                row_dict_2['museq_TA'] == 0)
    else:
        # File 1 is the one that's less supported
        return (total_tumour_1 > SUPPORT_AGAINST_MINIMUM_COVERAGE and
                row_dict_1['museq_TA'] == 0)


if __name__ == '__main__':
    main()
