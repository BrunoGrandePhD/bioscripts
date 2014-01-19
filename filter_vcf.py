#!/usr/bin/env python

"""
filter_vcf.py
~~~~~~~~~~~~~~

This Python script will output
Tested with Python 2.7.

Input:
The expected input are standard VCF file(s).

Output:
A summary file is generated with the SNV calls that match the specified
filter, organized by VCF file in which they're found. This is mostly useful
when you are looking at multiple VCF files.
"""


import argparse
import datetime


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    ## TODO Finish explaining the --filter option.
    parser.add_argument('--filter', '-f', help='<Help to be added.>',
                        action='store', required=True)
    parser.add_argument('--output', '-o', help='Instead of outputting to ' +
                        'the current directory, you can specify where you ' +
                        'want the summary file.', action='store')
    parser.add_argument('files', help='Space-separated list of VCF files.',
                        nargs='+')
    args = parser.parse_args()
    filter_dict = parse_filter(args.filter)
    if args.output:
        output_file_path = args.output
    else:
        date_and_time = datetime.datetime.today().strftime('%Y-%m-%d_%H-%M')
        output_file_path = 'filter_vcf--' + date_and_time + '.vcf'
    files = args.files

    # Generate the output files
    with open(output_file_path, 'w') as opened_output_file:
        for current_input_file in files:
            with open(current_input_file) as opened_input_file:
                apply_filter(opened_input_file, opened_output_file,
                             filter_dict)


def apply_filter(input_file_handle, output_file_handle, filter_dict):
    output_file_handle.write('# ' + input_file_handle.name + '\n')
    for row_dict in parse_vcf_generator(input_file_handle):
        if (row_dict['chromosome'] == filter_dict['chromosome'] and
                row_dict['position'] > filter_dict['min_position'] and
                row_dict['position'] < filter_dict['max_position']):
            output_file_handle.write(createVcfRow(row_dict))
    output_file_handle.write('\n')
    return


def parse_filter(raw_filter):
    filter_dict = {}
    chromosome, positions = raw_filter.split(':')
    min_position, max_position = positions.split('-')
    filter_dict['chromosome'] = chromosome
    filter_dict['min_position'] = int(min_position)
    filter_dict['max_position'] = int(max_position)
    return filter_dict

    # if len(raw_filter_split) == 1 or (len(raw_filter_split) == 2 and
    #                                   (raw_filter_split[1].strip() == '' or
    #                                    raw_filter_split[1].strip() == '-')):
    #     if raw_filter_split[0].find('chr') == :
    #         filter_dict['chromosome'] = raw_filter_split[0]
    #     elif:

    #     filter_dict['min_position'] = 0
    #     filter_dict['max_position'] = float('inf')
    # elif (len(raw_filter_split) == 2 and
    #       len(raw_filter_split[1].split('-')) == 2):
    #     filter_dict['chromosome'] = raw_filter_split[0]
    #     position_split = raw_filter_split[1].split('-')
    # else:
    #     raise Exception('Your filter doesn\'t match the expected format, ' +
    #                     'i.e., chromosome_id:position1-position2. ' +
    #                     'See the --help instructions for more info.')


def parse_vcf_generator(vcf_file_handle):
    for line in vcf_file_handle:
        if line.strip()[0] is '#':
            continue
        columns = line.split('\t')
        row_dict = {
            'chromosome': columns[0],
            'position': int(columns[1]),
            'id': columns[2],
            'reference_allele': columns[3],
            'alternate_allele': columns[4],
            'quality': float(columns[5]),
            'filter': columns[6],
            'info': columns[7]
        }
        yield row_dict


def createVcfRow(row_dict):
    row_list = [row_dict['chromosome'], str(row_dict['position']),
                row_dict['id'], row_dict['reference_allele'],
                row_dict['alternate_allele'], str(row_dict['quality']),
                row_dict['filter'], row_dict['info']]
    return '\t'.join(row_list)


if __name__ == '__main__':
    main()
