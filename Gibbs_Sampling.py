from __future__ import division
import random
import math

class Bivariate_Gibbs:

    def __init__(self, adduct_cluster_data, alpha, iterations, file_name):
        self.adduct_cluster_data = adduct_cluster_data
        self.peaks = self.adduct_cluster_data.get_peaks()
        self.only_one_cluster = self.adduct_cluster_data.get_only_one_cluster()
        self.more_than_one_cluster = self.adduct_cluster_data.get_more_than_one_cluster()
        self.K = len(adduct_cluster_data.get_clusters())
        self.alpha = alpha/self.K
        self.iterations = iterations
        #self.output_file_name = 'Adducts_Output.txt'
        self.output_file_name = file_name+'_RUN_JF.txt'

    '''
    Function to update the cluster count and cluster sum variables for a given cluster.
    @:param The Possible_Cluster object corresponding to the Cluster whose attributes are to be updated.
    Boolean variable add_to_cluster - True if peak is being added and False if being removed from cluster.
    '''
    def update_cluster_count(self, current_possible_cluster, add_to_cluster):
        # mod = +/- 1 depending on whether peak is being added/removed from current_possible_cluster
        if add_to_cluster:
            mod = 1
        else:
            mod = -1

        current_cluster = current_possible_cluster.get_cluster()
        transformed_mass = current_possible_cluster.get_transformed_mass()
        RT_value = current_possible_cluster.get_RT_value()

        # updated cluster count and sum values
        current_cluster.set_cluster_count(current_cluster.get_cluster_count()+mod)
        current_cluster.set_cluster_mass_sum(current_cluster.get_cluster_mass_sum()+mod*transformed_mass)
        current_cluster.set_cluster_RT_sum(current_cluster.get_cluster_RT_sum()+mod*RT_value)

    '''
    Function to updated the current_possible_cluster parameter for a given peak.
    @:param The peak being considered, a list of its Possible_Cluster objects, the updated
    probabilities that the peak is in each Possible_Cluster.
    '''
    def multinomial_sample(self, peak, possible_clusters, probs):
        # get new current cluster
        u = random.random()
        cum_prob = 0

        for cluster in possible_clusters:
            cum_prob += probs[cluster]
            if u <= cum_prob:
                peak.current_possible_cluster = cluster
                break
    '''
    Function to write the Gibbs sampler output to a text file.
    '''
    def write_to_output(self):
        print("Obtained allocations for peaks with more than one cluster.")
        print("Writing to output")
        output_file = open(self.output_file_name, 'w')
        output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format("Peak Mass", "RT", "Intensity", "Cluster", "Transform", "Freq.", "Prob.", "M"))
        for peak in self.peaks:
            current_cluster = peak.get_current_possible_cluster()
            possible_clusters = [(cluster.get_cluster_number(), cluster.get_transform().get_name()) for cluster in peak.get_possible_clusters()]
            '''
            output_file.write("\n Adduct {}, un-transformed mass and RT: ({},{}), allocated to cluster {}, transformed mass {}, cluster mass {}, \n possible clusters {}, "
                              "transformation {}, times in cluster {}, prob in cluster {} \n".format(
                peak.get_number(), peak.get_mass(), peak.get_RT(), current_cluster.get_cluster_number(),
                current_cluster.get_transformed_mass(), current_cluster.get_cluster().get_mass_mean(), possible_clusters,
                current_cluster.get_transform().get_name(), current_cluster.get_times_peak_in_cluster(), current_cluster.get_prob_peak_in_cluster()))
            '''

            output_file.write("\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(peak.get_mass(), peak.get_RT(), peak.get_intensity(), current_cluster.get_cluster_number(),
                                                            current_cluster.get_transform().get_name(), current_cluster.get_times_peak_in_cluster(), current_cluster.get_prob_peak_in_cluster(),
                                                                        current_cluster.get_transformed_mass()))

        output_file.close()

    '''
    Main function which run the Gibbs sampler for all of the peaks.

    The Gibbs sampler algorithm
    is then run on peaks with more than one possible cluster.

    Writes output to text file once finished.
    '''
    def run(self):
        # set parameter values for peaks with only one possible cluster
        print("Getting allocations for peaks with only one cluster ...")
        for peak in self.only_one_cluster:
            # set current_possible_cluster for peaks with only one cluster
            peak.set_current_possible_cluster(peak.get_possible_clusters()[0])

            current_possible_cluster = peak.get_current_possible_cluster()

            self.update_cluster_count(current_possible_cluster, True)

            current_possible_cluster.set_times_peak_in_cluster(self.iterations/2)
            current_possible_cluster.set_prob_peak_in_cluster(1.0)

        print("Obtained allocations for peaks with only one cluster.")

        print("Getting allocations for peaks with more than one cluster ...")
        for i in range(self.iterations):
            print("iteration {}".format(i))
            for peak in self.more_than_one_cluster:
                possible_clusters = peak.get_possible_clusters()
                current_possible_cluster = peak.get_current_possible_cluster()

                if current_possible_cluster != None:
                    self.update_cluster_count(current_possible_cluster, False)

                total_cluster_sum_minus_peak = 0

                for possible_cluster in possible_clusters:
                    total_cluster_sum_minus_peak += possible_cluster.get_cluster().get_cluster_count()

                log_probs = {}

                # calculate updated log probabilities for each possible_cluster for peak
                for possible_cluster in possible_clusters:
                    cluster = possible_cluster.get_cluster()
                    transformed_mass = possible_cluster.get_transformed_mass()
                    RT_value = possible_cluster.get_RT_value()
                    mu_0_mass = cluster.get_mass_mean()
                    mu_0_RT = cluster.get_RT_mean()
                    sigma = cluster.get_sigma()
                    sigma_0 = cluster.get_sigma_0()
                    cluster_count = cluster.get_cluster_count()
                    cluster_mass_sum = cluster.get_cluster_mass_sum()
                    cluster_RT_sum = cluster.get_cluster_RT_sum()

                    b_x = math.sqrt((math.pow(sigma[0]*sigma_0[0], 2)
                                       /(math.pow(sigma_0[0], 2)*cluster_count+math.pow(sigma[0], 2))))
                    b_y = math.sqrt((math.pow(sigma[1]*sigma_0[1], 2)
                                           /(math.pow(sigma_0[1], 2)*cluster_count+math.pow(sigma[1], 2))))

                    a_x = math.pow(b_x, 2)*(mu_0_mass/math.pow(sigma_0[0], 2)
                                                        +(cluster_mass_sum)/math.pow(sigma[0], 2))
                    a_y = math.pow(b_y, 2)*(mu_0_RT/math.pow(sigma_0[1], 2)
                                                        +(cluster_RT_sum)/math.pow(sigma[1], 2))

                    alpha_plus_cluster_count = (self.alpha + cluster_count)/(self.alpha*self.K + total_cluster_sum_minus_peak)

                    data_term = math.pow((transformed_mass-a_x), 2)/(2*(math.pow(sigma[0],2)+math.pow(b_x, 2)))\
                                +math.pow((RT_value-a_y), 2)/(2*(math.pow(sigma[1],2)+math.pow(b_y, 2)))

                    pi_term = 2*math.pi*math.sqrt((math.pow(sigma[0],2)+math.pow(b_x, 2))
                                                  *(math.pow(sigma[1], 2)+math.pow(b_y, 2)))

                    log_probs[possible_cluster] = math.log(alpha_plus_cluster_count)-math.log(pi_term)-data_term

                # get max log probability for current peak over all possible clusters
                prob_max = max(log_probs.values())

                # store each log value minus prob_max
                log_probs_minus_max = {}

                for cluster in possible_clusters:
                    log_probs_minus_max[cluster] = math.exp(log_probs[cluster]-prob_max)

                probs = {}

                # calculate probability current peak is in each cluster
                for cluster in possible_clusters:
                    number = cluster.get_cluster_number()
                    probs[cluster] = (math.exp(log_probs[cluster]-prob_max))/sum(log_probs_minus_max.values())

                # get new current_possible_cluster based on updated probabilties for each possible_cluster
                self.multinomial_sample(peak, possible_clusters, probs)

                # update cluster count and sum variables for peak's updated current cluster
                current_possible_cluster = peak.get_current_possible_cluster()

                self.update_cluster_count(current_possible_cluster, True)

                # start recording number of times peak is in current_possible_cluster after burn-in period of iterations/2
                if (i+1) > self.iterations/2:
                    current_possible_cluster.set_times_peak_in_cluster(current_possible_cluster.get_times_peak_in_cluster()+1)
                    current_possible_cluster.set_prob_peak_in_cluster(current_possible_cluster.get_times_peak_in_cluster()/(self.iterations/2))

        # write output to text file Adducts_Output.txt
        self.write_to_output()

        print("DONE!")