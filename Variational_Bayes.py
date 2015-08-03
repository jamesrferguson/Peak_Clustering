from __future__ import division
import math
from scipy.special import psi

class Variational_Bayes:

    def __init__(self, adduct_cluster_data, alpha, iterations):
        self.adduct_cluster_data = adduct_cluster_data
        self.peaks = self.adduct_cluster_data.get_peaks()
        self.only_one_cluster = self.adduct_cluster_data.get_only_one_cluster()
        self.more_than_one_cluster = self.adduct_cluster_data.get_more_than_one_cluster()
        self.K = len(adduct_cluster_data.get_clusters())
        self.alpha = alpha/self.K
        self.iterations = iterations
        self.output_file_name = 'Adducts_Output_VB.txt'

    def run(self):
        # set parameter values for peaks with only one possible cluster
        print("Getting allocations for peaks with only one cluster ...")
        for peak in self.only_one_cluster:
            peak.set_current_possible_cluster(peak.get_possible_clusters()[0])

            possible_clusters = peak.get_possible_clusters()
            for possible_cluster in possible_clusters:
                possible_cluster.set_expected_z_value(1)

                cluster = possible_cluster.get_cluster()
                cluster.set_z_sum(1.0)
                cluster.set_z_sum_mass(possible_cluster.get_transformed_mass())
                cluster.set_z_sum_RT(possible_cluster.get_RT_value())
                cluster.set_updated_alpha(self.alpha)

        print("Obtained allocations for peaks with only one cluster.")

        for i in range(self.iterations):
            for peak in self.more_than_one_cluster:
                possible_clusters = peak.get_possible_clusters()
                for possible_cluster in possible_clusters:
                    # 1. Get alpha values for each cluster
                    cluster = possible_cluster.get_cluster()
                    cluster.set_updated_alpha(self.alpha+cluster.get_z_sum())

                    # 2. Get updated mean parameters for each cluster
                    mu_0_mass = cluster.get_mass_mean()
                    mu_0_RT = cluster.get_RT_mean()
                    sigma = cluster.get_sigma()
                    sigma_0 = cluster.get_sigma_0()
                    z_sum = cluster.get_z_sum()
                    z_sum_mass = cluster.get_z_sum_mass()
                    z_sum_RT = cluster.get_z_sum_RT()

                    b_x = (math.pow(sigma[0]*sigma_0[0], 2))/(math.pow(sigma_0[0], 2)*z_sum+math.pow(sigma[0], 2))
                    b_y = (math.pow(sigma[1]*sigma_0[1], 2))/(math.pow(sigma_0[1], 2)*z_sum+math.pow(sigma[1], 2))

                    a_x = b_x*((z_sum_mass)/(math.pow(sigma[0], 2)) + mu_0_mass/(math.pow(sigma_0[0], 2)))
                    a_y = b_y*((z_sum_RT)/(math.pow(sigma[1], 2)) + mu_0_RT/(math.pow(sigma_0[1], 2)))

                    cluster.set_updated_mean([a_x, a_y])
                    cluster.set_updated_sd([b_x, b_y])

            # update <z_nk> values for each cluster
            for peak in self.more_than_one_cluster:
                possible_clusters = peak.get_possible_clusters()
                alpha_sum = 0.0
                for possible_cluster in possible_clusters:
                    cluster = possible_cluster.get_cluster()
                    alpha_sum += cluster.get_updated_alpha()

                gamma_values = {}

                for possible_cluster in possible_clusters:
                    cluster = possible_cluster.get_cluster()
                    sigma = cluster.get_sigma()
                    sigma_0 = cluster.get_sigma_0()
                    cluster_mass_mean = cluster.get_updated_mean()[0]
                    cluster_RT_mean = cluster.get_updated_mean()[1]
                    cluster_mass_sd = cluster.get_updated_sd()[0]
                    cluster_RT_sd = cluster.get_updated_sd()[1]
                    alpha = cluster.get_updated_alpha()
                    mass = possible_cluster.get_transformed_mass()
                    RT = possible_cluster.get_RT_value()

                    gamma_values[possible_cluster] = psi(alpha)+psi(alpha_sum) \
                                                     - (math.pow(mass, 2) - 2*mass*cluster_mass_mean + math.pow(cluster_mass_sd, 2) + math.pow(cluster_mass_mean, 2))/(2*math.pow(sigma[0], 2)) \
                                                     - (math.pow(RT, 2) - 2*RT*cluster_RT_mean+math.pow(cluster_RT_sd, 2) + math.pow(cluster_RT_mean, 2))/(2*math.pow(sigma[1], 2)) \
                                                     - math.log(2*math.pi*sigma[0]*sigma[1])

                max_value = max(gamma_values.values())

                # store each log value minus prob_max
                gamma_values_minus_max = {}

                for cluster in possible_clusters:
                    gamma_values_minus_max[cluster] = math.exp(gamma_values[cluster]-max_value)

                cluster_z = {}

                for possible_cluster in possible_clusters:
                    cluster = possible_cluster.get_cluster()
                    cluster_z[possible_cluster] = math.exp(gamma_values[possible_cluster]-max_value)/sum(gamma_values_minus_max.values())
                    peak.set_current_z(possible_cluster, cluster_z[possible_cluster])

                    cluster.set_z_sum(cluster.get_z_sum()+cluster_z[possible_cluster])
                    cluster.set_z_sum_mass(cluster.get_z_sum_mass()+cluster_z[possible_cluster]*possible_cluster.get_transformed_mass())
                    cluster.set_z_sum_RT(cluster.get_z_sum_RT()+cluster_z[possible_cluster]*possible_cluster.get_RT_value())

                peak.set_current_possible_cluster(max(cluster_z, key=cluster_z.get))

        print("Obtained allocations for peaks with more than one cluster.")
        print("Writing to output")
        output_file = open(self.output_file_name, 'w')

        for peak in self.peaks:
            current_cluster = peak.get_current_possible_cluster()
            possible_clusters = [(cluster.get_cluster_number(), cluster.get_transform().get_name()) for cluster in peak.get_possible_clusters()]
            output_file.write("\n Adduct {}, un-transformed mass and RT: ({},{}), allocated to cluster {}, transformed mass {}, cluster mass {}, \n possible clusters {}, "
                              "transformation {} \n".format(
                peak.get_number(), peak.get_mass(), peak.get_RT(), current_cluster.get_cluster_number(),
                current_cluster.get_transformed_mass(), current_cluster.get_cluster().get_mass_mean(), possible_clusters,
                current_cluster.get_transform().get_name()))