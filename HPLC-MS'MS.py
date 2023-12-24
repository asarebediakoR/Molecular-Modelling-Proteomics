## HPLC-MS/MS

import os
import sys
import pymzml
import matplotlib.pyplot as plt

RT = []
INT = []
def main():
    example_file = os.path.join('Try_stepped.mzML') # Replace mzML file
    run = pymzml.run.Reader(example_file)
   
    for n, spec in enumerate(run):
        rt = spec.scan_time_in_minutes()
        INT.append(spec.TIC)
        RT.append(rt)
        print(rt)
        if rt < 30 and rt > 25: 
            if spec.ms_level != 1:
                continue
            peaks_tuples = spec.peaks("centroided") # peaks function returns tuples of mz and intensity
            import pprint
            pprint.pprint(peaks_tuples)
            total_intensities = 0
            for mz,i in peaks_tuples:
               total_intensities += 1   # summing up the intensities for each spectrum
            print(rt,total_intensities)
            break
       
if __name__ == "__main__":
    main()

# Plot total ion chromatogram
plt.plot(RT,INT)
plt.xlabel('Time (min)')
plt.ylabel('Intensity')
plt.title('Total Ion Chromatogram (TIC)')
plt.show()


import os
import pymzml
def main():
   
    example_file = os.path.join('Try_stepped.mzML') # loading the mzml with the Reader function
    run = pymzml.run.Reader(example_file)
    # iterating over run, this will allow you to go through the different spectra in the file
    p = pymzml.plot.Factory()
    for n, spec in enumerate(run):
        rt = spec.scan_time_in_minutes()
        print(rt)
        if rt < 30 and rt > 25: 
            if spec.ms_level != 1:  # checking the ms level of each spectrum
                continue
            p.new_plot()
            p.add(spec.peaks("centroided"), color=(0, 0, 0), style="sticks", name="peaks")
            filename = "MS1_plot_{0}_{1}.html".format(
            os.path.basename(example_file), spec.ID
            )
            p.save(filename=filename)
            print("Plotted file: {0}".format(filename))
            break


if __name__ == "__main__":
    main()


import pymzml
import os
import pyteomics
from pyteomics import mass

# Calculate m/z for VTAHILSVGR with charge 2
peptide = 'VTAHILSVGR'
charge = 2
Mass = mass.calculate_mass(sequence='VTAHILSVGR', charge = 2)
print('m/z',Mass)










def main():
    
    example_file = os.path.join('Try_stepped.mzML')
    mz_to_find = 526.81419003328745
    run = pymzml.run.Reader(example_file)
    for spectrum in run:
        found_peaks = spectrum.has_peak(mz_to_find)
        if found_peaks != []:
            print("Found peaks: {0} in spectrum {1}".format(found_peaks, spectrum.ID))


if __name__ == "__main__":
    main()
    

import os
import pymzml
import pandas as pd

time_dependent_intensities = []
def main():
    
    example_file = os.path.join('Try_stepped.mzML')
    run = pymzml.run.Reader(example_file)
    
    My_mz = 525.806913861975

    for spectrum in run:  ## Iterating through all spectra for the specific peptide
        if spectrum.ms_level == 1:
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
plt.title("Intensities for VTAHILSVGR2+ peptide over time ")
plt.show()
plt.savefig('Intensities_of_all_Ions_over_Time')