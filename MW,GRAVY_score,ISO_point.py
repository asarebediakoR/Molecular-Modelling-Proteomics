
## Python packages used: Biopython,matplotlib,pyteomics,numpy
## Input file :Streptococcus pneumoniae fasta file ((strain ATCC BAA-255 / R6))
 
  
## Estimating sequesnce ID, Molecular Weight and GRAVY score of fasta sequence
 
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
listAnalysis = []
for sequence in SeqIO.parse('Strepto.fasta','fasta'): # Replace with your fasta file name
    print(sequence.id) # Print sequence ID
    print(repr(sequence.seq)) #Print sequence 
    print(len(sequence))  #Print total nucleotides
    
for seq_record in SeqIO.parse('Strepto.fasta','fasta'): #Replace file name
    my_seq = (seq_record.seq) 
    analysed_seq = ProteinAnalysis(my_seq) 
    listAnalysis.append(analysed_seq.molecular_weight()) # Sends molecular weight results into empty 'listAnalysis' array
    print(seq_record.id, analysed_seq.gravy())    #Reports sequence id and gravy score
Analysis_results = analysed_seq.gravy()  
      
# Analysis of Isoelectric point

from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA
protein = PA("Strepto")
print("IEP of {}: {:.2f}".format(protein.sequence,
                                 protein.isoelectric_point())) # Print out the pI of sequence
print("Charge at pH 4.53: {:.2f}"
      .format(protein.charge_at_pH(4.53)))


### Plotting Gravy score of protein sequences
from matplotlib import pyplot as plt
import numpy as np

##dataset
a = listAnalysis
 
# Creating  the histogram with dataset
fig, ax = plt.subplots(figsize =(10, 7))
ax.hist(a, bins =int(180/5))
plt.show()   ## Shows plot

### Converting  results into a  CSV file

from Bio.Seq import Seq
from Bio import SeqIO        
header = "Sequence ID, Molecular Weight, Isoelectic Point, Gravy Score\n"
for records in SeqIO.parse("Strepto.fasta", "fasta"): #Replace with file name
    analysed_seq = ProteinAnalysis(records.seq)
    data = f"{records.id}, {analysed_seq.molecular_weight()},{analysed_seq.isoelectric_point()},{analysed_seq.gravy()}\n"
    header  = header + (data)
print(header)

f = open("Strepto.csv", "w")
f.write(header)
f.close()
      
