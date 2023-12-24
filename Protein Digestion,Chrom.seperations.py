
## Input file :Streptococcus pneumoniae fasta file ((strain ATCC BAA-255 / R6))

import Bio
import csv
import pprint
from Bio import SeqIO
from pyteomics import parser
from collections import defaultdict
import matplotlib.pyplot as plt

# Parsing the FASTA file
tryptic_peptides = defaultdict(list)

for records in SeqIO.parse('Strepto.fasta','fasta'): ## Replace fasta file
    sequence = str(records.seq)
    fasta_id = str(records.id)
    for pep in parser.cleave(sequence, 'trypsin'):
        tryptic_peptides[pep].append(fasta_id)
print(tryptic_peptides)
pprint.pprint(tryptic_peptides)  # Prints all digest fragments with ids
     
  ## Creating CSV file output    
header = 'seq,modications'
My_peptides = f'{header}\n'
for peptides in tryptic_peptides.keys():
    # print(peptides)
    My_peptides = My_peptides + f'{str(peptides)}\n'
    
f = open("analysis.csv", "w+")
f.write(My_peptides)
f.close()

# General counts for individual peptide sequences from digest
Dict = {}
for i in tryptic_peptides:
    val= i
    if(val in Dict):
          Dict[val]+= 1
    else:
        Dict[val]= 1
print(Dict)

# Non-proteotypic peptides count
non_proteotypic_peptides = [i for i, count in Dict.items() if count > 1]
print(non_proteotypic_peptides)

# Percentage of non-proteotypic peptides
percent_non_proteotypic = len(non_proteotypic_peptides) / len(tryptic_peptides)
print(percent_non_proteotypic)


# Plot for the non-proteotypic peptides
non_proteotypic_peptide_lengths = [len(i) for val in non_proteotypic_peptides]
plt.hist(non_proteotypic_peptide_lengths, bins=range(min(non_proteotypic_peptide_lengths), max(non_proteotypic_peptide_lengths) + 1))
plt.xlabel('Peptide Length')
plt.ylabel('Count')
plt.title('Peptide Lengths for Non-Proteotypic Peptides')
plt.show()











    
