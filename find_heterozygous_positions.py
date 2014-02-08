#!/usr/bin/env python

# This script takes one argument, a SNVMix2 output file.
# This script will read each line of the SNVMix2 output
# and filter it down to positions that have a >0.99
# probability of being heterozygous.
# The output of this script is a .POS file, which is a
# tab-separated format with two columns, chromosome and
# position. This output file will be saved alongside
# the input SNVMix file.

import sys
import os

inputFilePath = sys.argv[1]
inputFileDir, inputFileName = os.path.split(inputFilePath)
inputFile = open(inputFilePath, "r")
# Empty output file
open(inputFileDir + "/" + os.path.splitext(inputFileName)[0] + ".filt_AB.pos", "w").close()
outputFile = open(inputFileDir + "/" + os.path.splitext(inputFileName)[0] + ".filt_AB.pos", "a")

for line in inputFile:
    columns = line.split("\t")
    if float(columns[3].split(",")[3]) > 0.99:
        chromosome, position = columns[0].split(":")
        outputFile.write(chromosome + ":" + position + "\n")

inputFile.close()
outputFile.close()
