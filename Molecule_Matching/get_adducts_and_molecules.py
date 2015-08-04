from __future__  import division
from Adducts_Molecules import *
import csv
import re

class Molecule_Transforms:
    def __init__(self, molecules_file, adducts_files, output_file_name, mass_interval):
        self.molecules_file = molecules_file
        self.adducts_files = adducts_files
        self.output_file_name = output_file_name
        self.mass_interval = mass_interval
        self.molecules = []
        self.adducts = []

    def run(self):
        # read molecule details from csv file
        with open(self.molecules_file, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            for row in reader:
                if (row[0] != ''):
                    self.molecules.append(Molecule(row[0], row[1], row[2], float(row[3]), self.mass_interval))

        csvfile.close()

        # read adduct details from peak cluster model output
        for file_name in self.adducts_files:
            print(file_name)
            file = open(file_name, 'r')

            line = file.readline()
            # ignore 1st line of each file
            line = file.readline()

            while line != "":
                line = line.rstrip('\n')
                params = re.split('\t', line)
                self.adducts.append(Adduct(float(params[7]), params[4]))
                line = file.readline()

            file.close()

        # loop through each molecule and check if its mass is within mass_interval of each adduct mass
        for molecule in self.molecules:
            for adduct in self.adducts:
                molecule.check_adduct(adduct)

        with open(self.output_file_name, 'wb') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(["Standard", "Name", "Formula", "Mass", ' ', "Transforms"])
            for molecule in self.molecules:
                adducts = molecule.get_transforms()

                list = [molecule.get_standard(), molecule.get_name(), molecule.get_formula(),
                                 molecule.get_mass(),' ']

                for transform in adducts:
                    list.append(transform)
                    list.append(adducts[transform])
                    list.append(adducts[transform]/molecule.get_adduct_count())

                writer.writerow(list)

        csvfile.close()

'''
RUN
'''

stndrd1_molecules_filename = "std1_mols.csv"
stndrd2_molecules_filename = "std2_mols.csv"

stndrd1_adducts_files = ["std1-file1.group.peakml_RUN_JF.txt", "std1-file2.group.peakml_RUN_JF.txt", "std1-file3.group.peakml_RUN_JF.txt",
                 "std1-file4.group.peakml_RUN_JF.txt", "std1-file5.group.peakml_RUN_JF.txt"]

stndrd2_adducts_files = ["std2-file1.group.peakml_RUN_JF.txt", "std2-file2.group.peakml_RUN_JF.txt", "std2-file3.group.peakml_RUN_JF.txt",
                 "std2-file4.group.peakml_RUN_JF.txt", "std2-file5.group.peakml_RUN_JF.txt"]

stndrd1 = Molecule_Transforms(stndrd1_molecules_filename, stndrd1_adducts_files, "stndrd1_molecules_with_transforms.csv", 0.01)
stndrd2 = Molecule_Transforms(stndrd2_molecules_filename, stndrd2_adducts_files, "stndrd2_molecules_with_transforms.csv", 0.01)

stndrd1.run()
stndrd2.run()