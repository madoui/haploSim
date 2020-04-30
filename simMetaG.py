#!/usr/bin/python
##################################################
## simMetaG.py simulated metagenomic sample from haplotypes
##################################################
## Author: Amin Madoui
## Copyright: Copyright 2020, CEA-Genoscope
## License: CC-BY-NC
## Version: 0.1.2
## Maintainer: Amin Madoui
## Email: amadoui@genoscope.cns.fr
##################################################

from numpy import loadtxt
import numpy as np
import sys,os,getopt
import glob

def getGenomeSize(filename):
    size = 0
    last_name = None
    with open(filename) as sequences:
        for line in sequences:
            if not line.startswith('>'):
                size = size + len(line[:-1])
    return size

def getHaplo (haploDir):
    pi_i_j = dict()
    haplo = dict()
    for spDir in os.scandir(haploDir):
        sp = spDir.path
        sp_pi = dict()
        # Get pairwise Pi values
        for line in open( sp + "/pairwise_pi.txt" ):
            data = line.split("\t")
            comp = str( data[0] ) + "_" + str( data[1] )
            sp_pi[comp] = data[2]
        pi_i_j[sp] = sp_pi
    return pi_i_j

def generateHaploAbundance( all_pi_i_j, n , haplodir, outputDir, read_length, minAbundance, sampleID):
    """ Generate n random haplotypes abundances and compute pi (nucleotide diversity) """
    for sp in all_pi_i_j:
        spp = sp.split("/")
        species = spp[1]
        readsPath = outputDir+"/"+str(sampleID)+"/"+species
        try:
            os.mkdir(outputDir)
            os.mkdir(outputDir+"/"+str(sampleID))
            os.mkdir(outputDir+"/"+str(sampleID)+"/"+species)
        except FileExistsError:
            print("Some files already exist...")
        pi_i_j = all_pi_i_j[sp]
        abundances = [int(np.random.random()*100)+ minAbundance for i in range(n)]
        hapFreq  = [abundances[i]/sum(abundances) for i in range(len(abundances))]
        genomeSize = getGenomeSize( sp + "/seeds_0.fasta")
        log = open(readsPath+"/log.txt", "a")
        # Generate read sets
        for i in range(n):
            read1 = readsPath + "/seeds_" + str(i) + "_read1.fq"
            read2 = readsPath + "/seeds_" + str(i) + "_read2.fq"
            input =  sp + "/seeds_" + str(i) + ".fasta"
            nb_pairs = int( abundances[i] * genomeSize / (2*read_length ))
            log.write ( "wgsim " + input + " " + read1 + " " + read2 + " -1 " + str(read_length) + " -2 " + str(read_length) + " -N " + str( nb_pairs )+"\t#Abundance: "+str(abundances[i])+"\n" )
            os.system ( "wgsim " + input + " " + read1 + " " + read2 + " -1 " + str(read_length) + " -2 " + str(read_length) + " -N " + str( nb_pairs ) )
        # compute pi
        pi = 0
        for i in range(n-1):
            for j in range(i+1,n):
                pi = pi + float(hapFreq[i]) * float(hapFreq[j]) * float(pi_i_j[str(i)+"_"+str(j)])
        print ("Pi: "+str(pi))
        log.write("\npi\t" + str(pi))
        log.close()
        print ( "Species: " + sp + ", pi: " + str(pi) )

def main(argv):
   inputHaplo = ''
   outputDir = 'simulation'
   haploNumber = 5
   readLength = 100
   minAbundance = 10
   usageMessage = """
   simMetaG.py -i inputHaplo -o outputDir -n haploNumber -m <...more options...>

   -i input directory that contains haplotypes fasta files
   -o output directory
   -n minimum of haplotypes per species (default 5)
   -l read leangth (default 100)
   -m minimum abundance of an haplotype (default 10x)

   """
   try:
      opts, args = getopt.getopt(argv,"hi:o:n:l:m:")
   except getopt.GetoptError:
      print (usageMessage)
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print (usageMessage)
         sys.exit()
      elif opt in ("-i"):
         inputHaplo = arg
      elif opt in ("-o"):
         outputDir = arg
      elif opt in ("-n"):
        haploNumber = int(arg)
      elif opt in ("-l"):
        readLength = int(arg)
      elif opt in ("-m"):
        minAbundance = int(arg)
   if not os.path.exists(inputHaplo):
       print (inputHaplo, " is not a file")
       print (usageMessage)
       sys.exit(2)
   print ('Parameters:')
   print ('Input directory: ', inputHaplo)
   print ('Output directory: ', outputDir)
   print ('Number of final haplotypes per sample: ', haploNumber)

   pi_i_j = getHaplo(inputHaplo)

   sampleID = 1
   generateHaploAbundance(pi_i_j, haploNumber, inputHaplo, outputDir, readLength, minAbundance, sampleID)

main(sys.argv[1:])
