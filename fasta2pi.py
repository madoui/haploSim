#!/usr/bin/python
##################################################
## fasta2pi.py compute the nucleotide diversity (pi) from a fasta file
##################################################
## Author: Amin Madoui
## Copyright: Copyright 2020, Commissariat à l'Energie Atomique (CEA)
## Centre National de Séquençage - Genoscope
## License: CC-BY-NC
## Version: 0.0.1
## Maintainer: Amin Madoui
## Email: amadoui@genoscope.cns.fr
##################################################
## input:
## output:
import os, sys, getopt, os.path
from collections import OrderedDict

def fasta2dict(fastafile):
    # Parses a fasta file into a dictionary (cc) of dictionary (path).
    result = Dict()
    last_name = None
    with open(fastafile) as sequences:
        for line in sequences:
            if line.startswith('>'):
                last_seqname = line[1:-1]
                result[last_seqname] = []
            else:
                result[last_seqname].append(line[:-1])
    for name in result:
        result[name] = ''.join(result[name])
    return result

def fasta2pi (fasta_infile):
    fasta_sequences = parse_sequences(fasta_infile)
    seeds = dict()
    # foreach species
    for fasta_name in fasta_sequences:

def main(argv):
   inputFasta = ''
   outputDir = 'fasta2pi_output'
   usageMessage = """
   fasta2pi.py -i inputFasta -o output_prefix

   -i input fasta file
   -o output prefix

   """
   try:
      opts, args = getopt.getopt(argv,"hi:o:")
   except getopt.GetoptError:
      print (usageMessage)
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print (usageMessage)
         sys.exit()
      elif opt in ("-i"):
         inputFasta = arg
      elif opt in ("-o"):
         outputDir = arg

   if not os.path.exists(inputFasta):
       print (inputFasta," is not a file")
       print (usageMessage)
       sys.exit(2)
