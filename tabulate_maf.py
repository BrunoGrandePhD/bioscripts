#!/usr/bin/env python
import csv
import sys
import os
import re

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
    input_maf = open(sys.argv[1])
    output_file_with_snp = open(
        os.path.splitext(os.path.basename(sys.argv[1]))[0]+'.summary.tsv', 'w'
    )
    output_file_without_snp = open(
        os.path.splitext(os.path.basename(
            sys.argv[1]
        ))[0]+'.noSNP.summary.tsv', 'w'
    )
    maf_reader = csv.DictReader(
        input_maf,
        fieldnames=MAF_FIELDNAMES,
        delimiter='\t')
    maf_writer_with_snp = csv.DictWriter(
        output_file_with_snp,
        fieldnames=[
            'gene',
            'number_of_silent_mutations',
            'number_of_missense_mutations',
            'number_of_other_mutations',
            'number_of_nonsilent_mutations'
        ],  delimiter='\t')
    maf_writer_with_snp.writeheader()
    maf_writer_without_snp = csv.DictWriter(
        output_file_without_snp,
        fieldnames=[
            'gene',
            'number_of_silent_mutations',
            'number_of_missense_mutations',
            'number_of_other_mutations',
            'number_of_nonsilent_mutations'
        ],  delimiter='\t')
    maf_writer_without_snp.writeheader()
    maf_rows = list(maf_reader)
    maf_rows.sort(key=lambda row: row['Entrez_Gene_Id'])

    def tabulate(maf_rows, withoutSnp):
        past_Entrez_Gene_Id = None
        cumulative_results = {
            'gene': None,
            'number_of_silent_mutations': 0,
            'number_of_missense_mutations': 0,
            'number_of_other_mutations': 0,
            'number_of_nonsilent_mutations': 0
        }
        for row in maf_rows:
            current_Entrez_Gene_Id = row['Entrez_Gene_Id']
            if current_Entrez_Gene_Id == "":
                continue
            if withoutSnp:
                skip = False
                identifiers = row['dbSNP_RS'].split('&')
                for i in identifiers:
                    if re.search('^rs.*', i):
                        skip = True
                    elif re.search('COS.*', i):
                        skip = False
                        break
                if skip:
                    continue
            if current_Entrez_Gene_Id == past_Entrez_Gene_Id:
                if row['Variant_Classification'] == 'Silent':
                    cumulative_results['number_of_silent_mutations'] += 1
                elif row['Variant_Classification'] == 'Missense_Mutation':
                    cumulative_results['number_of_missense_mutations'] += 1
                else:
                    cumulative_results['number_of_other_mutations'] += 1
            else:
                # If the next gene is new, write current cumulative_results
                # to file.
                if cumulative_results['gene'] is not None:
                    cumulative_results['number_of_nonsilent_mutations'] = (
                        cumulative_results['number_of_missense_mutations'] +
                        cumulative_results['number_of_other_mutations'])
                    if withoutSnp:
                        maf_writer_without_snp.writerow(cumulative_results)
                    else:
                        maf_writer_with_snp.writerow(cumulative_results)
                # Reset cumulative_results
                cumulative_results['gene'] = current_Entrez_Gene_Id
                cumulative_results['number_of_silent_mutations'] = 0
                cumulative_results['number_of_missense_mutations'] = 0
                cumulative_results['number_of_other_mutations'] = 0
                if row['Variant_Classification'] == 'Silent':
                    cumulative_results['number_of_silent_mutations'] += 1
                elif row['Variant_Classification'] == 'Missense_Mutation':
                    cumulative_results['number_of_missense_mutations'] += 1
                else:
                    cumulative_results['number_of_other_mutations'] += 1
                past_Entrez_Gene_Id = current_Entrez_Gene_Id

    tabulate(maf_rows, False)
    tabulate(maf_rows, True)
    input_maf.close()
    output_file_with_snp.close()
    output_file_without_snp.close()


if __name__ == '__main__':
    main()
