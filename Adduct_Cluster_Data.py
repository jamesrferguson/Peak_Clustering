from __future__ import division
from Get_Data import *
from Adduct_Details import *
from Gibbs_Sampling import *
from Variational_Bayes import *

class Adduct_Cluster_Data:
    def __init__(self, data_file_name, transform_file_name):
        self.RT_INTERVAL = 10
        self.MASS_INTERVAL = 1e-6
        self.PPM = 5.0

        self.data_file_name = data_file_name
        self.transform_file_name = transform_file_name
        self.peaks = []
        self.only_one_cluster = []
        self.more_than_one_cluster = []
        self.transforms = {}
        self.clusters = []

        # load adduct data
        self.data = read_data(data_file_name)
        # load adduct transformation details
        self.transformation_data = read_transformations(transform_file_name)

        # populate list of peaks
        print("Getting peaks ...")
        peak_count = 0
        for entry in self.data:
            self.peaks.append(Peak(peak_count, entry[0], entry[1], entry[2]))
            peak_count += 1
        print("Got peaks")

        # populate list of transformations
        print("Getting transforms ...")
        for transform in self.transformation_data:
            new_transform = Transform(transform[0], transform[1:])
            self.transforms[new_transform.get_name()] = new_transform
        print("Got transforms")

        # populate clusters by applying M+H transform to each adduct in adducts
        print("Getting clusters ...")
        cluster_number = 0
        for peak in self.peaks:
            self.clusters.append(Cluster(cluster_number, self.transforms['M+H'].apply(peak),
                                         peak.get_RT(), self.RT_INTERVAL, self.PPM,
                                         peak.get_intensity()))
            cluster_number += 1
        print("Got clusters")

        # get possible clusters for each Peak
        print("Finding possible clusters for peaks ...")
        for peak in self.peaks:
            for cluster in self.clusters:
                if (cluster.peak_in_RT_bin(peak) & cluster.peak_in_intensity_bound(peak)):
                    if (self.transforms['M+H'].apply(peak) != cluster.get_mass_mean()):
                        for transform in self.transforms.values():
                            if (cluster.peak_in_mass_bin(peak, transform)):
                                peak.add_possible_cluster(Possible_Cluster(cluster, transform, transform.apply(peak), peak.get_RT()))
                    else:
                        peak.add_possible_cluster(Possible_Cluster(cluster, self.transforms['M+H'], self.transforms['M+H'].apply(peak), peak.get_RT()))

        print("Finished finding possible clusters. Find peaks with only one possible cluster ....")

        for peak in self.peaks:
            possible_clusters = peak.get_possible_clusters()
            if (len(possible_clusters) == 1):
                self.only_one_cluster.append(peak)
            else:
                self.more_than_one_cluster.append(peak)

        print("Finished!")

    def get_peaks(self):
        return self.peaks

    def get_clusters(self):
        return self.clusters

    def get_only_one_cluster(self):
        return self.only_one_cluster

    def get_more_than_one_cluster(self):
        return self.more_than_one_cluster

    def get_transforms(self):
        return self.transforms

'''
RUN DATA
'''

alpha = 1
iterations = 1000

extension_data = "Data/"
extension_output = "Data/Output/"
'''
files = ['std1-file1.group.peakml', 'std1-file2.group.peakml','std1-file3.group.peakml','std1-file4.group.peakml',
         'std1-file5.group.peakml','std2-file1.group.peakml','std2-file2.group.peakml','std2-file3.group.peakml',
         'std2-file4.group.peakml','std2-file5.group.peakml']
'''
files = ['testtxt']
'''
for file in files:

    print("Running {}".format(extension_data+file))

    add_cluster_data  = Adduct_Cluster_Data(extension_data+file+".txt", "mulsub2.txt")

    gibbs = Bivariate_Gibbs(add_cluster_data, alpha, iterations, extension_output+file)

    gibbs.run()

print("DONE")
'''
for file in files:
    print("Running {}".format(extension_data+file))

    add_cluster_data  = Adduct_Cluster_Data(extension_data+file+".txt", "mulsub2.txt")

    vb = Variational_Bayes(add_cluster_data, alpha, iterations, extension_output+file)

    vb.run()
