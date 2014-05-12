#!/usr/bin/env python

"""
maf_update.py
~~~~~~~~~~~~~
This Python script outputs a MAF file given a filtered down VCF file
(the one which gave rise to the MAF file in question).
The use case is having one master MAF file containing everything and
using the information in the corresponding VCF file to filter down the
variants (e.g., according to allele ratios).
Effectively, you're filtering down (or updating) the MAF file according
to filters applied to the VCF file.
"""


import argparse
import os
import bioscripts


def main():
    parser = argparse.ArgumentParser(description='Updates MAF file' +
                                     'according to a filtered down' +
                                     'VCF file.')
    parser.add_argument('-m', '--maf',
                        type=argparse.FileType('r'), nargs=1, required=True,
                        help='Specify the input MAF file which will be ' +
                        'updated according to the filtered down VCF file. ' +
                        'The MAF file must be sorted.')
    parser.add_argument('-v', '--vcf',
                        type=argparse.FileType('r'), nargs=1, required=True,
                        help='Specify the input VCF that has been filtered ' +
                        'down. The VCF file must be sorted.')
    parser.add_argument('-o', '--output',
                        type=argparse.FileType('w'), nargs=1,
                        help='Specify where to output the updated MAF file.')
    args = parser.parse_args()
    input_maf = args.maf[0]
    input_vcf = args.vcf[0]
    if args.output:
        output_maf = args.output[0]
    else:
        root, ext = os.path.splitext(input_maf.name)
        output_maf_name = root + '.updated.' + ext
        if not os.path.isfile(output_maf_name):
            output_maf = open(output_maf_name, 'w')
        else:  # The default output MAF file already exists
            raise Exception('Default output file already exists. ' +
                            'Please specify an output file with the ' +
                            '-o/--output option.')

    maf_row = next(input_maf, None)
    while maf_row[0] is '#':
        maf_row = next(input_maf, None)
    vcf_row = next(input_vcf, None)
    while vcf_row[0] is '#':
        vcf_row = next(input_vcf, None)
    while maf_row is not None and vcf_row is not None:
        maf_row_dict = bioscripts.parse_maf(maf_row)
        vcf_row_dict = bioscripts.parse_vcf(vcf_row)

        print ('(maf) ' + maf_row_dict['chromosome'] + ':' +
               str(maf_row_dict['position']) + ' vs (vcf) ' +
               vcf_row_dict['chromosome'] + ':' +
               str(vcf_row_dict['position']))

        if (maf_row_dict['chromosome'] == vcf_row_dict['chromosome'] and
                maf_row_dict['position'] == vcf_row_dict['position']):
            # There is a match between the VCF and MAF files and needs
            # to be written to the output MAF file
            print 'match'
            output_maf.write(maf_row_dict['row'])
            maf_row = next(input_maf, None)
            vcf_row = next(input_vcf, None)

        elif (maf_row_dict['chromosome'] > vcf_row_dict['chromosome'] or
              (maf_row_dict['chromosome'] == vcf_row_dict['chromosome'] and
               maf_row_dict['position'] > vcf_row_dict['position'])):
            # The VCF file lags behind the MAF file and needs to catch up
            vcf_row = next(input_vcf, None)

        else:
            # The MAF file is behind and needs to catch up
            maf_row = next(input_maf, None)

    input_maf.close()
    input_vcf.close()
    output_maf.close()


if __name__ == '__main__':
    main()
