
## HPLC-MS/MS (MS1,MS2 Measurements)

import os
import pymzml
import plotly.graph_objects as go

def main():
    example_file = os.path.join('Try_stepped.mzML') #Replace mzMLfile
    run = pymzml.run.Reader(example_file)
    retention_times = []
    total_ion_intensities = []
    for n, spec in enumerate(run):
        if spec.ms_level != 1:  # Checks the ms level
         continue
        rt = spec.scan_time_in_minutes()
        retention_times.append(rt)
        
        tic = spec.TIC
        total_ion_intensities.append(tic)
    print(retention_times) 
    print(total_ion_intensities)  


    fig = go.Figure()
    fig.add_trace(
    go.Scatter(
        x=retention_times,
        y=total_ion_intensities,
        mode='lines',
    )
)    
    fig.update_layout(title='Total Ion Chromatogram(TIC)')
    fig.show()
if __name__ == "__main__":
    main()    

import pprint
from Bio import SeqIO
from Bio.Seq import Seq
from pyteomics import  mass,parser
from collections import defaultdict
import matplotlib.pyplot as plt

# Parsing the FASTA file
tryptic_peptides = defaultdict(list)

for records in SeqIO.parse('CSG.fasta','fasta'):
    sequence = str(records.seq)
    fasta_id = str(records.id)
    for pep in parser.cleave(sequence, 'trypsin'):
        tryptic_peptides[pep].append(fasta_id)
pprint.pprint(tryptic_peptides)  # Prints all digest fragments with ids

selected_peptides= [p for p in tryptic_peptides if 5 <= len(p) <= 40] 
pprint.pprint(selected_peptides)    

## Calculating mass of  selected peptides
for count, i in enumerate(selected_peptides):
    Pep_weight = mass.calculate_mass(sequeance=i, charge=2)
    # pprint.pprint('weight',Pep_weight)
    x=[]
    y=[]
    for n,spectrum in enumerate(run):
        if spectrum.ms_level == 1:
          retention_time = spectrum.scan_time_in_minutes()
          peakes_tuple = spectrum.peaks('centroided')
          for mz,intensity in peakes_tuple:
              diff = abs(mz-Pep_weight)
              x.append(retention_time)
              y.append(intensity)
        pprint.pprint(x)
        pprint.pprint(y)
              
from matplotlib import pyplot as plt
import numpy as np

fig, ax = plt.subplots()
ax.scatter(retention_time, intensity)
plt.xlabel('Retention Time') 
plt.ylabel('Intensity') 
plt.title("Spectra of all peptides")
plt.show()
plt.savefig('Trptic digest peptides')
        
              
import os
import pymzml
import pandas as pd

time_dependent_intensities = []
def main():
    
    example_file = os.path.join('Try_stepped.mzML') # Replace  mzMLfile
    run = pymzml.run.Reader(example_file)
    
    My_mz = 526.814190328745

    for spectrum in run:  ##Iterating through all spectra for the specific peptide
        if spectrum.ms_level == 2: # Checking for ms levels
            has_peak_matches = spectrum.has_peak(My_mz)
            if has_peak_matches != []:
                for mz, I in has_peak_matches:
                    time_dependent_intensities.append(
                        [spectrum.scan_time_in_minutes(), I, mz]
                    )
    print("RT   \ti   \tmz")
    for rt, i, mz in time_dependent_intensities:
        print("{0:5.3f}\t{1:13.4f}\t{2:10}".format(rt, i, mz))
    return


if __name__ == "__main__":
    main()
    
from matplotlib import pyplot as plt
import numpy as np

df = pd.DataFrame (time_dependent_intensities, columns = (['Retention Time', 'Intensity', 'M/Z']))
retent_time = df["Retention Time"]
intensity = df['Intensity']

fig, ax = plt.subplots()
ax.scatter(retent_time, intensity)
plt.xlabel('Retention Time') 
plt.ylabel('Intensity') 
plt.title("MS2 for VTAHILSVGR2+ peptide")
plt.show()
plt.savefig('Intensities_of_all_Ions_over_Time')


import copy

def main():
    #  defining general layout attributes
    layout = {
        "xaxis": {
            "title": "<i>m/z</i>",
            "tickmode": "auto",
            "showticklabels": True,
            "ticklen": 5,
            "tickwidth": 1,
            "ticks": "outside",
            "showline": True,
            "showgrid": False,
        },
        "yaxis": {
            "color": "#000000",
            "tickmode": "auto",
            "showticklabels": True,
            "ticklen": 5,
            "tickwidth": 1,
            "ticks": "outside",
            "showline": True,
            "showgrid": False,
        },
    }

    example_file = os.path.join('Try_stepped.mzML' ) # Replace file

    # Precisions for MS2
    run = pymzml.run.Reader(example_file, MS_precisions={2:5e-4})
    p = pymzml.plot.Factory()
    plot_layout = {}


    # Specific MS2 spectrum
    ms2_spectrum = run[8016]

    # The MS_precision for the plotting option label.triangle.MS_precision
    p.new_plot(title="MS2 Spectrum Annotated: VTAHILSVGR",MS_precision=5e-4)
    p.add(
        ms2_spectrum.peaks("centroided"),
        color=(0, 0, 0),
        opacity=1,
        style="sticks",
        name="Peptide_spectrum",
    )

    b_ions = {
        "b<sub>2</sub><sup>+2</sup>": 59.54676576,
        "b<sub>3</sub><sup>+2</sup>": 110.0706049,
        "b<sub>4</sub><sup>+2</sup>": 145.5891618,
        "b<sub>2</sub>": 214.118617,
        "b<sub>5</sub><sup>+2</sup>": 270.660649,
        "b<sub>6</sub><sup>+2</sup>": 327.2026817,
        "b<sub>3</sub>": 370.718696,
        "b<sub>7</sub><sup>+2</sup>": 420.2529029,
        "b<sub>4</sub>": 448.763634816, 
    }

    y_ions = {
        "y<sub>1</sub><sup>+2</sup>": 477.2799833,
        "y<sub>2</sub><sup>+2</sup>": 426.75614413,
        "y<sub>1</sub>": 391.2375872,
        "y<sub>3</sub><sup>+2</sup>": 322.7081313,
        "y<sub>4</sub><sup>+2</sup>": 266.1660993,
        "y<sub>2</sub>": 209.6240673393,
        "y<sub>5</sub><sup>+2</sup>": 166.108053137,
        "y<sub>6</sub><sup>+2</sup>": 116.573846180,
        "y<sub>3</sub>": 88.063114320419,
    }

    
    # Using the has_peak() function
    for ion_list in [b_ions, y_ions]:
        label_list = []
        for fragment in ion_list.keys():
            peak = ms2_spectrum.has_peak(ion_list[fragment])
            if len(peak) != 0:
                label_list.append((ion_list[fragment], peak[0][1], fragment))
        if ion_list == b_ions:
            color = (0, 0, 255)
        else:
            color = (0, 255, 0)
        p.add(
            label_list,
            color=color,
            style="label.triangle.MS_precision",
            name="b,y-ions fragments",
        )

    for axis in layout.keys():
        plot_layout["{0}3".format(axis)] = copy.copy(layout[axis])

    
    filename = "Ions_plot_{0}_annotation.html".format(os.path.basename(example_file))
    p.save(filename=filename, layout=plot_layout)
    print("Plotted file: {0}".format(filename))


if __name__ == "__main__":
    main()