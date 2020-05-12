#!/usr/bin/python
##################################################
## simHaplo.py simulated haplotypes from a genome according
##  to a random binary evolutionary tree
##################################################
## Author: Amin Madoui
## Copyright: Copyright 2020, CEA-Genoscope
## License: CC-BY-NC
## Version: 0.1.2
## Maintainer: Amin Madoui
## Email: amadoui@genoscope.cns.fr
##################################################

import os, sys, getopt, os.path
import numpy as np
import pprint
import random
from collections import OrderedDict
from typing import Dict
import math

def parse_sequences(filename: str, ordered: bool=False) -> Dict[str, str]:
    """
    Parses a text file of genome sequences into a dictionary.
    Arguments:
      filename: str - The name of the file containing the genome info.
      ordered: bool - Set this to True if you want the result to be ordered.
    """
    result = OrderedDict() if ordered else {}
    last_name = None
    with open(filename) as sequences:
        for line in sequences:
            if line.startswith('>'):
                last_name = line[1:-1]
                print (last_name)
                result[last_name] = []
            else:
                result[last_name].append(line[:-1])
    for name in result:
        result[name] = ''.join(result[name])
    return result


def draw(discrete_probdist):
    """
    Draw random value from discrete probability distribution
    represented as a dict: P(x=value) = discrete_probdist[value].
    """
    limit = 0
    r = np.random.random()
    for value in discrete_probdist:
        limit += discrete_probdist[value]
        if r < limit:
            return value


def create_markov_chain():
    markov_chain = {}
    for from_base in 'ATGC':
       slice_points = sorted([0] + [np.random.random()for i in range(3)] + [1])
       transition_probabilities = [slice_points[i+1] - slice_points[i] for i in range(4)]
       markov_chain[from_base] = {base: p for base, p in zip('ATGC', transition_probabilities)}
    return markov_chain

def mutate_via_markov_chain(dna, markov_chain):
    """ Take the DNA sequence and the markov chain and returns a one-point mutated DNA sequence"""
    dna_list = list(dna)
    mutation_site = np.random.randint(0, len(dna_list) - 1)
    from_base = dna[mutation_site]
    to_base = draw(markov_chain[from_base])
    dna_list[mutation_site] = to_base
    return ''.join(dna_list)

def createSeeds (fasta_file_name, vRate, vDisp, nSeedLines):
    """ Take a fasta file, a mutation rate, a mutation dipersion and a number of mutated genomes to create"""
    fasta_sequences = parse_sequences(fasta_file_name)
    seeds = dict()
    # foreach species
    for fasta_name in fasta_sequences:
        new_sequences = []
        name, sequence = fasta_name, fasta_sequences[fasta_name]
        print ("Species: "+name)
        mc = create_markov_chain() # create a species specific nucleotide transition matrix (markov chain)
        n = int(math.log(nSeedLines,2)) # the number of iteration along a complete binary tree to create nSeedLines
        new_sequences.append(sequence)
        # Generate mutated sequences along a complete binary tree
        for i in range(2,n+2):
            print ("Tree depth: "+str(i))
            leaves_sequences = []
            for j in range(2**(i-1),(2**i),2):#the last leaves
                print ("Node: "+str(i))
                for k in range(2):# gives two new leaves (mutated sequences)
                    mutated_sequences = muteSeq (new_sequences[int(j/2)-1], vRate/(n-1), vDisp, mc)
                    new_sequences.append(mutated_sequences)
                    leaves_sequences.append (mutated_sequences)
            if i == n+1:
                seeds[name] = leaves_sequences
    return seeds

def muteSeq (sequence, vRate, vDisp, mc):
     """ Take a DNA seq and return a mutated DNA seq """
     new_vRate = random.gauss(vRate, vDisp) # Generate a mutation rate following a gaussian distrib
     if new_vRate<0:
         new_vRate = 0.001
     nmutations = int(len(sequence)*new_vRate)
     print ("Generate "+str(nmutations)+" mutations")
     new_sequence = sequence
     for j in range(nmutations):
         new_sequence = mutate_via_markov_chain(new_sequence, mc)
     return new_sequence

def hamming_distance(s1, s2):
    """ Calculate the hamming distance between two strings"""
    h_dist = 0
    for ch1, ch2 in zip(s1, s2):
        if (ch1 != ch2):
            h_dist = h_dist+1
    return h_dist

def get_nuc_diverge_rate_dict (outputDir, seeds):
    """ Compute the pairwise nucleotide divergence rate between of all seeds """
    nuc_diverge_rate_dict = dict()
    for sp in seeds:
        f = open(outputDir+"/"+sp+"/pairwise_pi.txt", "a")
        species_nuc_diverge_rate_dict = dict()
        for seed1 in range(len(seeds[sp])):
            for seed2 in range(len(seeds[sp])):
                seq1 = seeds[sp][seed1]
                seq2 = seeds[sp][seed2]
                species_nuc_diverge_rate_dict[str(seed1)+"_"+str(seed2)] = hamming_distance(seq1,seq2)/len(seq1)
                f.write (str(seed1)+"\t"+str(seed2)+"\t"+str(hamming_distance(seq1,seq2)/len(seq1))+"\n")
        nuc_diverge_rate_dict[sp] = species_nuc_diverge_rate_dict
        f.close()
    return nuc_diverge_rate_dict

def createHaplotypes (seeds, recombFreq, haploNumber):
    """ Generate recombinated haplotypes"""
    haplotypes = dict()
    for species in seeds:
        print ("Species: ", species)
        seedSize = len(seeds[species][0])
        nbSeed = len(seeds[species])
        nb_recombination = int(seedSize*recombFreq)
        print ("Nb recombination: ",nb_recombination)
        if nb_recombination <1:
            print ("Number of recombination: ",nb_recombination, " you need to change the recombination rate (-r option)")
            sys.exit(2)
        for haplo in range(1,haploNumber+1):
            print ("Haplo: ", haplo)
            randnums = np.random.randint(1,101,nb_recombination)
            randseed = np.random.randint(1,nbSeed,nb_recombination)
            start = 0
            haploSeq = ''
            print (seedSize)
            for recomb in range(nb_recombination):
                #print ("Recombination number: ",recomb+1)
                haploFragmentSize = randnums[recomb]/100
                #print ("Size: ",str(start),"-",str(int(start+(haploFragmentSize*seedSize)))," - seed:",randseed[recomb])
                start = int(start+(haploFragmentSize*seedSize))+1

def write_seq(outputDir, Seeds):
    # Create output directory
    try:
        os.mkdir(outputDir)
        print("Directory " , outputDir ,  " created ")
    except FileExistsError:
        print("Directory " , outputDir ,  " already exists")
    # write fasta
    for sp in Seeds:
        try:
            os.mkdir(outputDir+"/"+sp)
            print("Directory " , outputDir+"/"+sp ,  " created ")
        except FileExistsError:
            print("Directory " , outputDir+"/"+sp ,  " already exists")
        for seeds in range(len(Seeds[sp])):
            f = open(outputDir+"/"+sp+"/seeds_"+str(seeds)+".fasta", "a")
            f.write(">seq_"+str(seeds)+"\n"+Seeds[sp][seeds]+"\n")
            f.close()

def main(argv):
   inputFasta = ''
   outputDir = ''
   varFreq = 0.01
   varDisp = 0.001
   recombFreq = 0.0001
   strainsNumber = 16
   haploNumber = 10
   readLength = 100
   minAbundance = 10
   usageMessage = """
   simHaplo.py -i inputFasta -o outputDir -k varFreq -r recombFreq -s strainsNumber -n haploNumber <...more options...>

   -i input directory that contains fasta file
   -o output directory
   -k variant frequency (default 0.01)
   -d variant frequency dispersion (default 0.001)
   -r recombination frequency (default 0.0001)
   -s number of initial strains to generate (default 16)
   -n number of haplotypes to generate (default 10)
   -g genome size (default 1000000)
   -l read leangth (default 100)
   -m minimam abundance of an haplotype (default 10x)

   """
   try:
      opts, args = getopt.getopt(argv,"hi:o:k:v:r:s:n:d:g:l:m:")
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
      elif opt in ("-k"):
         varFreq = float(arg)
      elif opt in ("-s"):
        strainsNumber = int(arg)
      elif opt in ("-n"):
        haploNumber = int(arg)
      elif opt in ("-d"):
        varDisp = float(arg)
      elif opt in ("-r"):
        recombFreq = float(arg)
      elif opt in ("-g"):
        genomeSize = int(arg)
      elif opt in ("-l"):
        readLength = int(arg)
      elif opt in ("-m"):
        minAbundance = int(arg)
   if not os.path.exists(inputFasta):
       print (inputFasta," is not a file")
       print (usageMessage)
       sys.exit(2)
   print ('Parameters:')
   print ('Input directory: ', inputFasta)
   print ('Output directory: ', outputDir)
   print ('variant frequency: ', varFreq)
   print ('variant frequency dispersion: ', varDisp)
   print ('Recombination frequency: ', recombFreq)
   print ('Number of initial genomic lines to create: ', strainsNumber)
   print ('Number of final haplotypes: ', haploNumber)

   # Create seed n genomes with mutation at rate T, (n = strainsNumber, T = varFreq)
   Seeds = createSeeds (inputFasta, varFreq, varDisp, strainsNumber)

   # Write fasta
   write_seq (outputDir, Seeds)

   # TO DO: Add recombination between seeds
   # haplotypes = createHaplotypes (Seeds, recombFreq, haploNumber)

   # Compute Pi_i_j
   pi_i_j = get_nuc_diverge_rate_dict ( outputDir, Seeds )

main(sys.argv[1:])
