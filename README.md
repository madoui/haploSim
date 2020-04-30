# haploSim
Generate metagenomic data from an admixture of haplotypes.

## Quick example on Coronavirus

```
python3 simHaplo.py -i wuhan-hu1.fasta -o CoronaHaplo -s 16
python3 simMetaG.py -i CoronaHaplo -o  Corona_1_reads -n 10

```

## How to simualte haplotypes?
To simulate haplotypes, use the `simHaplo.py` script. Several parameters can be tuned as follow:
```
simHaplo.py -i inputFasta -o outputDir -k varFreq -r recombFreq -s strainsNumber -n haploNumber <...more options...>

   -i input directory that contains fasta file
   -o output directory
   -k variant frequency (default 0.01)
   -d variant frequency dispersion (default 0.001)
   -r recombination frequency (default 0.0001)
   -s number of initial strains to generate (default 3)
   -n number of haplotypes to generate (default 10)
   -g genome size (default 1,000,000)
   -l read length (default 100)
   -m minimum abundance of a haplotype (default 10x)
 ```

## How to simulate metagenomic reads from a haplotypes admixture?
Once the haplotypes have been created, use `simMetaG` to generate illumina reads on a random admixture of haplotypes. Several parameters can tuned as follow:
```
 simMetaG.py -i inputHaplo -o outputDir -n haploNumber ...more options...>

   -i input directory that contains haplotypes fasta files
   -o output directory
   -n minimum of haplotypes per species (default 5)
   -l read leangth (default 100)
   -m minimum abundance of an haplotype (default 10x)
```

