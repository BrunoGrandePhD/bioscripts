#!/usr/bin/env python

# This script takes one argument, a SNVMix2 output file. 
# This script will read each line of the SNVMix2 output
# and reformat each line to abide the input format
# required by TITAN, namely the following 6 columns:
# chr, position, ref base, ref count, non-ref base, 
# non-ref count. No header and as TSV. 

import sys
import os

inputFilePath = sys.argv[1]
inputFileDir, inputFileName = os.path.split(inputFilePath)
inputFile = open(inputFilePath, "r")
# Empty output file
open(inputFileDir + "/" + os.path.splitext(inputFileName)[0] + "-titan.txt", "w").close()
outputFile = open(inputFileDir + "/" + os.path.splitext(inputFileName)[0] + "-titan.txt", "a")

for line in inputFile:
    columns = line.split("\t")
    column3 = columns[3].split(",")
    chromosome, position = columns[0].split(":")
    ref_base, ref_count = column3[0].split(":")
    nonref_base, nonref_count = column3[1].split(":")
    outputFile.write(chromosome + "\t" + position + "\t" + ref_base + "\t" +
        ref_count + "\t" + nonref_base + "\t" + nonref_count + "\n")

inputFile.close()
outputFile.close()
