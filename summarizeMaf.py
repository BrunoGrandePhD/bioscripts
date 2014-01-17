#!/usr/bin/env python

import sys
import os
import re
import datetime

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


def prepareFile(filepath, outputFilenameAddition=None):
    """ A convenience function to open files for reading and,
    optionally, open files for writing if the outputFilenameAddition
    parameter is specified. If the input file is 'foo.txt' and
    outputFilenameAddition='bar', the output file will be foo.bar.txt
    """
    directory, filename = os.path.split(filepath)
    if directory:
        directory + '/'
    filename_without_ext, file_extension = os.path.splitext(filename)
    input_maf = open(filepath)
    if outputFilenameAddition:
        output_maf = open(directory + filename_without_ext + '.' +
                          outputFilenameAddition + file_extension, 'w')
        return input_maf, filename_without_ext.replace(' ', '_'), output_maf
    return input_maf, filename_without_ext.replace(' ', '_')


def isSnp(dbSNP_RS):
    isSnp = False
    identifiers = dbSNP_RS.split('&')
    for i in identifiers:
        if re.search('^rs.*', i):
            isSnp = True
        elif re.search('^COSM.*', i):
            isSnp = False
            break
    return isSnp


# Open file for results summary
date_and_time = datetime.datetime.today().strftime('%Y-%m-%d_%H-%M')
summary_file = open('MAFsummary_' + date_and_time + '.tsv', 'w')
summary_file.write('# Every called SNV was either a SNP or a true SNV.\n')
summary_file.write('# A true SNV was either silent or nonsilent.\n')
summary_file.write('# filename\ttotal\tnumber_snps\tnumber_snvs\t' +
                   'number_snvs_silent\tnumber_snvs_nonsilent\n')

args = sys.argv[1:]

# Option to create a SNP-less MAF file
isRemoveSnp = False
if '-r' in args:
    isRemoveSnp = True
    args.pop(args.index('-r'))

for maf_filepath in args:
    # try:
    if isRemoveSnp:
        input_maf, file_name, output_maf = prepareFile(maf_filepath,
                                                       'noSNP')
    else:
        input_maf, file_name = prepareFile(maf_filepath)
    # except IOError:
    #     continue
    cumulative_count = {
        'number_snps': 0,
        'number_snvs': 0,
        'number_snvs_silent': 0,
        'number_snvs_nonsilent': 0
    }
    for line in input_maf:
        columns = dict(zip(MAF_FIELDNAMES, line.split('\t')))
        lineIsSnp = isSnp(columns['dbSNP_RS'])
        if lineIsSnp:
            cumulative_count['number_snps'] += 1
        else:
            cumulative_count['number_snvs'] += 1
            if columns['effect'] == 'silent':
                cumulative_count['number_snvs_silent'] += 1
            elif columns['effect'] == 'nonsilent':
                cumulative_count['number_snvs_nonsilent'] += 1
            if isRemoveSnp:
                output_maf.write(line)
    total = cumulative_count['number_snps'] + cumulative_count['number_snvs']
    summary_file.write(
        file_name + '\t' +
        str(total) + '\t' +
        str(cumulative_count['number_snps']) + '\t' +
        str(cumulative_count['number_snvs']) + '\t' +
        str(cumulative_count['number_snvs_silent']) + '\t' +
        str(cumulative_count['number_snvs_nonsilent']) + '\n'
    )
    if isRemoveSnp:
        output_maf.close()

summary_file.close()
