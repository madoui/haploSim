# haploSim
Generate metagenomic data from an admixture of haplotypes.

## Quick example on Coronavirus

```
python3 simHaplo.py -i data/wuhan-hu1.fasta -o CoronaHaplo -s 16
python3 simMetaG.py -i CoronaHaplo -o  Corona_1_reads -n 10

```
