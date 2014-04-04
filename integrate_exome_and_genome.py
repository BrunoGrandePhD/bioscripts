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

The general pipeline is the following:
1) Call SNVs using MutationSeq 4 in paired mode with PR > 0.8 on both the
   exome and genome.
2) (Optional) Restrict SNVs to exonic regions.
3) Pool both sets of SNVs and generate a positions (.pos) file for MutationSeq
   and run MutationSeq in position-specific mode with a threshold of 0 on both
   the exome and genome.
4) Ensure that there are the same number of SNV calls in these new VCF files.
   Also, the positions should be in the same order in both VCF files.
5) Run integrate_exome_and_genome.py on these two files.

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
    parser = argparse.ArgumentParser(description='Integrate VCF files from ' +
                                     'an exome and genome.')
    parser.add_argument('-i', '--input', nargs=2, type=argparse.FileType('r'),
                        required=True,
                        help='Specify both VCF files, one from the exome ' +
                        'and the other from the genome, separated by a space.')
    parser.add_argument('-p', '--prefix', nargs=1, required=True,
                        help='Specify the prefix for the output files, where' +
                        'the prefix is the path and the beginning of the ' +
                        'file name, e.g., "/path/to/file" becomes ' +
                        '/path/to/file.best_calls.vcf, etc.')
    parser.add_argument('-t', '--type', nargs=1, choices=['mutationseq',
                        'strelka'], default='mutationseq',
                        help='')
    args = parser.parse_args()
    input_file_1, input_file_2 = args.input
    output_prefix = args.prefix[0]
    output_good_calls, output_poor_calls, output_bad_calls = \
        create_output_files(output_prefix)
    variant_caller = args.type[0]

    # Choose parser
    if variant_caller == 'mutationseq':
        parse_vcf = parse_museq_vcf
    elif variant_caller == 'strelka':
        parse_vcf = parse_strelka_vcf

    # Going through the VCF files line by line
    # Ensuring that lines are compared with one another
    # And lines that don't have a match are automatically considered good
    row_dict_1 = next(parse_vcf(input_file_1), None)
    row_dict_2 = next(parse_vcf(input_file_2), None)

    while row_dict_1 is not None or row_dict_2 is not None:

        # Whenever one of the rows is None, handle it here
        if row_dict_1 is None:
            output_good_calls.write(create_vcf_row(row_dict_2,
                                    variant_caller))
            row_dict_2 = next(parse_vcf(input_file_2), None)
            continue
        if row_dict_2 is None:
            output_good_calls.write(create_vcf_row(row_dict_1,
                                    variant_caller))
            row_dict_1 = next(parse_vcf(input_file_1), None)
            continue

        # Check if the rows match one another
        if (row_dict_1['chromosome'] == row_dict_2['chromosome'] and
                row_dict_1['position'] == row_dict_2['position']):

            # If they do match, store the best row, which will be the one
            # that will be stored in the output files
            if row_dict_1['probability'] > row_dict_2['probability']:
                best_call = row_dict_1
            else:
                best_call = row_dict_2

            # Germline Test (for bad_calls)
            if test_germline(row_dict_1, row_dict_2):
                # Write row to bad_calls
                output_bad_calls.write(create_vcf_row(best_call,
                                       variant_caller))

            # Partial Support Test (for good_calls)
            elif test_partial_support(row_dict_1, row_dict_2):
                # Write row to good_calls
                output_good_calls.write(create_vcf_row(best_call,
                                        variant_caller))

            # Support Against Test (for poor_calls)
            elif test_support_against(row_dict_1, row_dict_2):
                # Write row to poor_calls
                output_poor_calls.write(create_vcf_row(best_call,
                                        variant_caller))

            # If the rows reach this point, they're sorted to good_calls
            else:
                output_good_calls.write(create_vcf_row(best_call,
                                        variant_caller))

            # Increment both VCF files
            row_dict_1 = next(parse_vcf(input_file_1), None)
            row_dict_2 = next(parse_vcf(input_file_2), None)

        else:  # The rows don't match

            # Check if the chromosomes match
            if row_dict_1['chromosome'] == row_dict_2['chromosome']:

                # Now that the chromosomes match, check if the position match
                if row_dict_1['position'] > row_dict_2['position']:
                    output_good_calls.write(create_vcf_row(row_dict_2,
                                            variant_caller))
                    row_dict_2 = next(parse_vcf(input_file_2), None)

                else:  # row_dict_1['position'] < row_dict_2['position']
                    output_good_calls.write(create_vcf_row(row_dict_1,
                                            variant_caller))
                    row_dict_1 = next(parse_vcf(input_file_1), None)

            else:  # The chromosome don't match

                if row_dict_1['chromosome'] > row_dict_2['chromosome']:
                    output_good_calls.write(create_vcf_row(row_dict_2,
                                            variant_caller))
                    row_dict_2 = next(parse_vcf(input_file_2), None)

                else:  # row_dict_1['chromosome'] < row_dict_2['chromosome']:
                    output_good_calls.write(create_vcf_row(row_dict_1,
                                            variant_caller))
                    row_dict_1 = next(parse_vcf(input_file_1), None)

    input_file_1.close()
    input_file_2.close()
    output_good_calls.close()
    output_poor_calls.close()
    output_bad_calls.close()


def create_output_files(output_prefix):
    output_good_calls_name = output_prefix + '.good_calls.vcf'
    output_good_calls = open(output_good_calls_name, 'w')
    output_poor_calls_name = output_prefix + '.poor_calls.vcf'
    output_poor_calls = open(output_poor_calls_name, 'w')
    output_bad_calls_name = output_prefix + '.bad_calls.vcf'
    output_bad_calls = open(output_bad_calls_name, 'w')
    return output_good_calls, output_poor_calls, output_bad_calls


def parse_museq_vcf(opened_vcf_file):
    # Sample VCF row
    # 1  14057566    .   G   A   16.20   INDL    INFO
    for line in opened_vcf_file:
        if line.strip()[0] is '#':
            continue
        columns = line.strip().split('\t')
        row_dict = {
            'chromosome': columns[0],
            'position': int(columns[1]),
            'id': columns[2],
            'reference_allele': columns[3],
            'alternate_allele': columns[4],
            'quality': columns[5],
            'filter': columns[6],
            'info': columns[7],
        }
        # Sample MutationSeq info column
        #PR=0.98;TR=174;TA=27;NR=152;NA=0;TC=CGG;NI=2;ND=1
        info_items = columns[7].split(';')
        row_dict.update({
            'probability': float(info_items[0][3:]),
            'tumour_number_ref': int(info_items[1][3:]),
            'tumour_number_alt': int(info_items[2][3:]),
            'normal_number_ref': int(info_items[3][3:]),
            'normal_number_alt': int(info_items[4][3:])
        })
        yield row_dict


def parse_strelka_vcf(opened_vcf_file):
    for line in opened_vcf_file:
        if line.strip()[0] is '#':
            continue
        columns = line.strip().split('\t')
        row_dict = {
            'chromosome': columns[0],
            'position': int(columns[1]),
            'id': columns[2],
            'reference_allele': columns[3],
            'alternate_allele': columns[4],
            'quality': columns[5],
            'filter': columns[6],
            'info': columns[7],
            'format': columns[8],
            'normal': columns[9],
            'tumour': columns[10]
        }

        # Sample Strelka INFO column
        #NT=ref;QSS=153;QSS_NT=153;SGT=GG->GT;SOMATIC;TQSS=1;TQSS_NT=1
        # Sample Strelka FORMAT column
        #DP:FDP:SDP:SUBDP:AU:CU:GU:TU
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
        yield row_dict


def create_vcf_row(row_dict, variant_caller):
    row_list = [row_dict['chromosome'], str(row_dict['position']),
                row_dict['id'], row_dict['reference_allele'],
                row_dict['alternate_allele'], str(row_dict['quality']),
                row_dict['filter'], row_dict['info']]
    if variant_caller == 'strelka':
        row_list.extend([row_dict['format'], row_dict['normal'],
                        row_dict['tumour']])
    return '\t'.join(row_list) + '\n'


def test_germline(row_dict_1, row_dict_2):
    total_germline_1 = (row_dict_1['normal_number_ref'] +
                        row_dict_1['normal_number_alt'])
    total_germline_2 = (row_dict_2['normal_number_ref'] +
                        row_dict_2['normal_number_alt'])

    # Handle cases where there is no coverage in the germline
    if total_germline_1 != 0:
        germline_allele_ratio_1 = (float(row_dict_1['normal_number_alt']) /
                                   total_germline_1)
    else:
        germline_allele_ratio_1 = 0

    if total_germline_2 != 0:
        germline_allele_ratio_2 = (float(row_dict_2['normal_number_alt']) /
                                   total_germline_2)
    else:
        germline_allele_ratio_2 = 0

    return (germline_allele_ratio_1 > GERMLINE_THRESHOLD or
            germline_allele_ratio_2 > GERMLINE_THRESHOLD)


def test_partial_support(row_dict_1, row_dict_2):
    total_tumour_1 = (row_dict_1['tumour_number_ref'] +
                      row_dict_1['tumour_number_alt'])
    total_tumour_2 = (row_dict_2['tumour_number_ref'] +
                      row_dict_2['tumour_number_alt'])

    # Handle cases where there is no coverage
    if total_tumour_1 != 0:
        tumour_allele_ratio_1 = (float(row_dict_1['tumour_number_alt']) /
                                 total_tumour_1)
    else:
        tumour_allele_ratio_1 = 0

    if total_tumour_2 != 0:
        tumour_allele_ratio_2 = (float(row_dict_2['tumour_number_alt']) /
                                 total_tumour_2)
    else:
        tumour_allele_ratio_2 = 0

    if row_dict_1['probability'] > row_dict_2['probability']:
        # File 2 is the one that's less supported
        return tumour_allele_ratio_2 > PARTIAL_SUPPORT_THRESHOLD
    else:
        # File 1 is the one that's less supported
        return tumour_allele_ratio_1 > PARTIAL_SUPPORT_THRESHOLD


def test_support_against(row_dict_1, row_dict_2):
    total_tumour_1 = (row_dict_1['tumour_number_ref'] +
                      row_dict_1['tumour_number_alt'])
    total_tumour_2 = (row_dict_2['tumour_number_ref'] +
                      row_dict_2['tumour_number_alt'])
    if row_dict_1['probability'] > row_dict_2['probability']:
        # File 2 is the one that's less supported
        return (total_tumour_2 > SUPPORT_AGAINST_MINIMUM_COVERAGE and
                row_dict_2['tumour_number_alt'] == 0)
    else:
        # File 1 is the one that's less supported
        return (total_tumour_1 > SUPPORT_AGAINST_MINIMUM_COVERAGE and
                row_dict_1['tumour_number_alt'] == 0)


if __name__ == '__main__':
    main()
