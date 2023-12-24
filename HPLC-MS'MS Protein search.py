#HPLC-MS/MS run,Protein Search
#!/usr/bin/env python3
# encoding: utf-8

import ursgal
import os
import sys

def main():
    """
    Simple example script how to generate a target decoy database using one or
    multiple fasta files as input.
    Note:
        A 'shuffled peptide preserving cleavage sites' database is
        generated.
    usage:
        ./generate_target_decoy_from_fasta.py <fasta_database_1> <fasta_database_2> <fasta_database_n> ...
    """
    params = {
        "enzyme": "trypsin",
        # 'decoy_generation_mode': 'reverse_protein',
    }

    fasta_database_list = sys.argv[1:-1]
    target_db_name = sys.argv[-1]

    uc = ursgal.UController(params=params)

    new_target_decoy_db_name = uc.execute_misc_engine(
        input_file=fasta_database_list,
        engine="generate_target_decoy_1_0_0",
        output_file_name=target_db_name,
    )
    print("Generated target decoy database: {0}".format(new_target_decoy_db_name))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(main.__doc__)
    main()


import ursgal
import sys
import glob
import os


def main(folder=None, profile=None, target_decoy_database=None):
    """
    An example test script to search all mzML files which are present in the
    specified folder. The search is currently performed on 4 search engines
    and 2 validation engines.
    The machine profile has to be specified as well as the target-decoy
    database.
    usage:
        ./do_it_all_folder_wide.py <mzML_folder> <profile> <target_decoy_database>
    Current profiles:
        * 'QExactive+'
        * 'LTQ XL low res'
        * 'LTQ XL high res'
    """
    # define folder with mzML_files as sys.argv[1]
    mzML_files = []
    for mzml in glob.glob(os.path.join("{0}".format(folder), "*.mzML")):
        mzML_files.append(mzml)

    mass_spectrometer = profile

    # We specify all search engines and validation engines that we want to use in a list
    # (version numbers might differ on windows or mac):
    search_engines = [
        "omssa_2_1_9",
        #"xtandem_vengeance",
        #"msgfplus_v2016_09_16",
        # 'msamanda_1_0_0_6300',
        # 'myrimatch_2_1_138',
    ]

    validation_engines = [
        "percolator_3_4_0",
        #"qvality",
    ]

    # Modifications that should be included in the search
    all_mods = [
        "C,fix,any,Carbamidomethyl",
        "M,opt,any,Oxidation",
        # 'N,opt,any,Deamidated',
        # 'Q,opt,any,Deamidated',
        # 'E,opt,any,Methyl',
        # 'K,opt,any,Methyl',
        # 'R,opt,any,Methyl',
        "*,opt,Prot-N-term,Acetyl",
        # 'S,opt,any,Phospho',
        # 'T,opt,any,Phospho',
        # 'N,opt,any,HexNAc'
    ]

    # Initializing the Ursgal UController class with
    # our specified modifications and mass spectrometer
    params = {
        "database": target_decoy_database,
        "modifications": all_mods,
        "csv_filter_rules": [
            ["Is decoy", "equals", "false"],
            ["PEP", "lte", 0.01],
        ],
    }

    uc = ursgal.UController(profile=mass_spectrometer, params=params)

    # complete workflow:
    # every spectrum file is searched with every search engine,
    # results are validated (for each engine seperately),
    # validated results are merged and filtered for targets and PEP <= 0.01.
    # In the end, all filtered results from all spectrum files are merged
    for validation_engine in validation_engines:
        result_files = []
        for spec_file in mzML_files:
            validated_results = []
            for search_engine in search_engines:
                unified_search_results = uc.search(
                    input_file=spec_file,
                    engine=search_engine,
                )
                validated_csv = uc.validate(
                    input_file=unified_search_results,
                    engine=validation_engine,
                )
                validated_results.append(validated_csv)

            validated_results_from_all_engines = uc.execute_misc_engine(
                input_file=validated_results,
                engine="merge_csvs_1_0_0",
            )
            filtered_validated_results = uc.execute_misc_engine(
                input_file=validated_results_from_all_engines,
                engine="filter_csv_1_0_0",
            )
            result_files.append(filtered_validated_results)

        results_all_files = uc.execute_misc_engine(
            input_file=result_files,
            engine="merge_csvs_1_0_0",
        )


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(main.__doc__)
        sys.exit(1)
    main(
        folder=sys.argv[1],
        profile=sys.argv[2],
        target_decoy_database=sys.argv[3],
    )


import csv

filename = 'Try_stepped_xtandem_vengeance.csv'
column_names = ['Protein ID', 'Sequence']  # Names of the columns to count
count_dict = {col: {} for col in column_names}  # A dictionary to store the  various counts

with open(filename, 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        for col in column_names:
            value = row[col]
            if value not in count_dict[col]:
                count_dict[col][value] = 0
            count_dict[col][value] += 1

# Print the counts
for col, count_dict_col in count_dict.items():
    print(col)
    for value, count in count_dict_col.items():
        print(f" - {value}: {count}")


