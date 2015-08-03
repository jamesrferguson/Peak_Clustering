from __future__ import division
import math

class Adduct:
    def __init__(self, mass, transform):
        self.mass = mass
        self.transform = transform

    def get_mass(self):
        return self.mass

    def get_transform(self):
        return self.transform

class Molecule:
    def __init__(self, standard, name, formula, mass, mass_interval):
        self.standard = standard
        self.name = name
        self.formula = formula
        self.mass = mass
        self.mass_interval = mass_interval
        self.transforms = {}
        self.adduct_count = 0

    def get_standard(self):
        return self.standard

    def get_name(self):
        return self.name

    def get_formula(self):
        return self.formula

    def get_mass(self):
        return self.mass

    def check_adduct(self, adduct):
        # check if difference between masses is within mass_interval
        if (abs(self.mass - adduct.get_mass()) <= self.mass_interval):
            # if in mass_interval, adduct to transforms/increment this transform's count by 1

            if adduct.get_transform() in self.transforms:
                self.transforms[adduct.get_transform()] += 1
                self.adduct_count += 1
            else:
                self.transforms[adduct.get_transform()] = 1
                self.adduct_count += 1

    def get_transforms(self):
        return self.transforms

    def get_adduct_count(self):
        return self.adduct_count