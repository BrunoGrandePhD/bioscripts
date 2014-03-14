#!/usr/bin/env python

import sys
import os

directory, filename = os.path.split(sys.argv[1])
filename_without_ext, file_extension = os.path.splitext(filename)

if directory is not '':
    directory += '/'

input_vcf = open(sys.argv[1])
output_vcf = open(directory + filename_without_ext + '.prepared' +
                  file_extension, 'w')

for line in input_vcf:
    # Keep header, but don't process it
    if line[0] == '#':
        output_vcf.write(line)
        continue
    columns = line.split('\t')
    # Skip all non-standard chromosomes
    chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
                   '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
                   '21', '22', 'X', 'Y', 'MT']
    if columns[0] not in chromosomes:
        continue
    info = columns[7]
    info_items = info.split(';')
    # Update NR and NA with numbers instead of N/A
    info_items[3] = 'NR=20'
    info_items[4] = 'NR=0'
    updated_info = ';'.join(info_items)
    updated_columns = columns[:]
    updated_columns[7] = updated_info
    updated_line = '\t'.join(updated_columns)
    output_vcf.write(updated_line)

input_vcf.close()
output_vcf.close()
