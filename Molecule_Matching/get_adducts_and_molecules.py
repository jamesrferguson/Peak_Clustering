from __future__  import division
import math
from Adducts_Molecules import *
import csv
import re

INTERVAL = 0.01

stndrd1_molecules_filename = "std1_mols.csv"
stndrd2_molecules_filename = "std2_mols.csv"

stndrd1_adducts_files = ["std1-file1.group.peakml_RUN_JF.txt", "std1-file2.group.peakml_RUN_JF.txt", "std1-file3.group.peakml_RUN_JF.txt",
                 "std1-file4.group.peakml_RUN_JF.txt", "std1-file5.group.peakml_RUN_JF.txt"]

stndrd2_adducts_files = ["std2-file1.group.peakml_RUN_JF.txt", "std2-file2.group.peakml_RUN_JF.txt", "std2-file3.group.peakml_RUN_JF.txt",
                 "std2-file4.group.peakml_RUN_JF.txt", "std2-file5.group.peakml_RUN_JF.txt"]

# get standard 1 molecules
stndrd1_molecules = []

with open(stndrd1_molecules_filename, 'rb') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in reader:
        if (row[0] != ''):
            stndrd1_molecules.append(Molecule(row[0], row[1], row[2], float(row[3]), INTERVAL))

stndrd1_adducts = []

for file_name in stndrd1_adducts_files:
    print(file_name)
    file = open(file_name, 'r')

    line = file.readline()
    # ignore 1st line of each file
    line = file.readline()

    while line != "":
        line = line.rstrip('\n')
        params = re.split('\t', line)
        stndrd1_adducts.append(Adduct(float(params[7]), params[4]))
        line = file.readline()

    file.close()

for molecule in stndrd1_molecules:
    for adduct in stndrd1_adducts:
        molecule.check_adduct(adduct)

output_file_name = "Molecule_Adducts.txt"

output_file = open(output_file_name, 'w')

for molecule in stndrd1_molecules:
    adducts = molecule.get_transforms()

    transform_details = ""

    for adduct in adducts:
        transform_details += adduct+" {0:.2f}% ".format(adducts[adduct]/molecule.get_adduct_count()*100)

    output_file.write("\nstandard: {} name: {} formula: {} mass: {} \n transforms: {}".format(molecule.get_standard(), molecule.get_name(), molecule.get_formula(),
                                           molecule.get_mass(), transform_details))

output_file.close()