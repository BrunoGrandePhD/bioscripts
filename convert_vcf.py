#!/usr/bin/env python

"""
convert_vcf.py
~~~~~~~~~~~~~~

By default, this Python script generates POS position files from standard
VCF files. Optionally, it can output minimal BED files compatible with UCSC.
See below for more information about input and output file formats.
Tested with Python 2.7.

Input:
The expected input VCF file(s) must have the chromosome identifier and
the position as the first two tab-separated columns, respectively.
The script ignores the rest.

Output:
The POS position file format consists of two tab-separated columns,
the first being the chromosome identifier and the second being the position.
Optionally, you can output minimal BED files consisting of three tab-separated
columns: the chromosome identifier according to UCSC (e.g., chr17), the
starting position and the ending position. The starting position will be the
one stated in the input VCF file, while the ending position will simply be
the starting position incremented by 1.
"""


import argparse
import os


# Supported formats in lowercase
SUPPORTED_OUTPUT_FORMATS = ['pos', 'bed', 'museq.pos']


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--format', '-f', help='Specify the desired output ' +
                        'file format. Currently supported formats are: ' +
                        ', '.join(SUPPORTED_OUTPUT_FORMATS), action='store',
                        default='pos')
    parser.add_argument('files', help='Space-separated list of VCF files.',
                        nargs='+')
    args = parser.parse_args()
    format = args.format.lower()
    files = args.files

    # Check if specified format is supported
    if format not in SUPPORTED_OUTPUT_FORMATS:
        raise Exception('Specified format not supported. Currently ' +
                        'supported formats are: ' +
                        ', '.join(SUPPORTED_OUTPUT_FORMATS))

    # Generate the output files
    for current_input_file in files:
        current_output_file = (os.path.splitext(current_input_file)[0]
                               + '.' + format)
        with open(current_input_file) as opened_input_file:
            with open(current_output_file, 'w') as opened_output_file:
                if format == 'pos':
                    generate_pos(opened_input_file, opened_output_file)
                elif format == 'bed':
                    generate_bed(opened_input_file, opened_output_file)
                elif format == 'museq.pos':
                    generate_museq_pos(opened_input_file, opened_output_file)


def generate_pos(vcf_file_handle, output_file_handle):
    for line_items in parse_vcf_generator(vcf_file_handle):
        chromosome, position = line_items[0], line_items[1]
        output_file_handle.write(chromosome + '\t' + position + '\n')


def generate_museq_pos(vcf_file_handle, output_file_handle):
    for line_items in parse_vcf_generator(vcf_file_handle):
        chromosome, position = line_items[0], line_items[1]
        output_file_handle.write(chromosome + ':' + position + '\n')


def generate_bed(vcf_file_handle, output_file_handle):
    for line_items in parse_vcf_generator(vcf_file_handle):
        chromosome, position = line_items[0], line_items[1]
        output_file_handle.write('chr' + chromosome + '\t' + position + '\t' +
                                 str(int(position) + 1) + '\n')


def parse_vcf_generator(vcf_file_handle):
    for line in vcf_file_handle:
        if line.strip()[0] is '#':
            continue
        columns = line.split('\t')
        yield columns


if __name__ == '__main__':
    main()
