from __future__ import division
'''
Contains classes Cluster, Transform and Adduct used for MS cluster modelling
'''

class Cluster:
    def __init__(self, number, mass_mean, RT_mean, RT_interval, ppm, max_intensity):
        self.number = number
        self.mass_mean = mass_mean
        self.RT_mean = RT_mean
        # mass bin width for cluster model
        self.ppm = ppm
        self.sigma = [self.mass_mean*self.ppm/3e6, 1]
        self.sigma_0 = [self.mass_mean*self.ppm/3e6, 1]
        # RT bin width for cluster model
        self.RT_bin_width = RT_interval
        self.max_intensity = max_intensity
        self.cluster_count = 0
        self.cluster_mass_sum = 0.0
        self.cluster_RT_sum = 0.0

        # for VB only
        self.z_sum = 0.0
        self.z_sum_mass = 0.0
        self.z_sum_RT = 0.0
        self.updated_alpha = 0.0
        self.updated_mean = [self.mass_mean, self.RT_mean]
        self.updated_sd = self.sigma

    def get_number(self):
        return self.number

    def get_mass_mean(self):
        return self.mass_mean

    def get_RT_mean(self):
        return self.RT_mean

    def get_sigma(self):
        return self.sigma

    def get_sigma_0(self):
        return self.sigma_0

    def peak_in_RT_bin(self, peak):
        if (peak.get_RT() > self.RT_mean-self.RT_bin_width) & (peak.get_RT()
        < self.RT_mean+self.RT_bin_width):
            return True
        else:
            return False

    def peak_in_intensity_bound(self, peak):
        if (peak.get_intensity() <= self.max_intensity):
            return True
        else:
            return False

    def peak_in_mass_bin(self, peak, transform):
        if (abs((transform.apply(peak) - self.mass_mean)/self.mass_mean)) <= self.ppm * 1e-6:
            return True
        else:
            return False

    def get_cluster_count(self):
        return self.cluster_count

    def get_cluster_mass_sum(self):
        return self.cluster_mass_sum

    def get_cluster_RT_sum(self):
        return self.cluster_RT_sum

    def set_cluster_count(self, amount):
        self.cluster_count = amount

    def set_cluster_mass_sum(self, amount):
        self.cluster_mass_sum = amount

    def set_cluster_RT_sum(self, amount):
        self.cluster_RT_sum = amount

    def get_z_sum(self):
        return self.z_sum

    def set_z_sum(self, amount):
        self.z_sum = amount

    def get_z_sum_mass(self):
        return self.z_sum_mass

    def set_z_sum_mass(self, amount):
        self.z_sum_mass = amount

    def get_z_sum_RT(self):
        return self.z_sum_RT

    def set_z_sum_RT(self, amount):
        self.z_sum_RT = amount

    def get_updated_alpha(self):
        return self.updated_alpha

    def set_updated_alpha(self, amount):
        self.updated_alpha = amount

    def get_updated_mean(self):
        return self.updated_mean

    def set_updated_mean(self, amount):
        self.updated_mean = amount

    def get_updated_sd(self):
        return self.updated_sd

    def set_updated_sd(self, amount):
        self.updated_sd = amount

class Possible_Cluster:
    def __init__(self, cluster, transform, transformed_mass, RT_value):
        self.cluster = cluster
        self.cluster_number = cluster.get_number()
        self.transform = transform
        self.transformed_mass = transformed_mass
        self.RT_value = RT_value
        self.times_peak_in_cluster = 0
        self.prob_peak_in_cluster = 0.0

        # VB only
        self.expected_z_value = 0.0

    def get_cluster(self):
        return self.cluster

    def get_cluster_number(self):
        return self.cluster_number

    def get_transform(self):
        return self.transform

    def get_transformed_mass(self):
        return self.transformed_mass

    def get_times_peak_in_cluster(self):
        return self.times_peak_in_cluster

    def get_prob_peak_in_cluster(self):
        return self.prob_peak_in_cluster

    def set_times_peak_in_cluster(self, amount):
        self.times_peak_in_cluster = amount

    def set_prob_peak_in_cluster(self, amount):
        self.prob_peak_in_cluster = amount

    def get_RT_value(self):
        return self.RT_value

    def get_expected_z_value(self):
        return self.expected_z_value

    def set_expected_z_value(self, amount):
        self.expected_z_value = amount

class Transform:
    def __init__(self, name, params):
        self.name = name
        self.params = params

    def get_name(self):
        return self.name

    def apply(self, peak):
        return (peak.get_mass()-self.params[0])/self.params[1] + self.params[2]

class Peak:
    def __init__(self, number, adduct_mass, adduct_RT, adduct_intensity):
        self.number = number
        self.adduct_mass = adduct_mass
        self.adduct_RT = adduct_RT
        self.adduct_intensity = adduct_intensity
        self.possible_clusters = []
        self.current_possible_cluster = None

        # VB only
        self.current_z = {}

    def get_number(self):
        return self.number

    def get_mass(self):
        return self.adduct_mass

    def get_RT(self):
        return self.adduct_RT

    def get_intensity(self):
        return self.adduct_intensity

    def get_possible_clusters(self):
        return self.possible_clusters

    def add_possible_cluster(self, possible_cluster):
        self.possible_clusters.append(possible_cluster)

    def set_current_possible_cluster(self, current):
        self.current_possible_cluster = current

    def get_current_possible_cluster(self):
        return self.current_possible_cluster

    def get_current_z(self):
        return self.current_z

    def set_current_z(self, cluster, z_value):
        self.current_z[cluster] = z_value