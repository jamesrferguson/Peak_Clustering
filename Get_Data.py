import re

'''
Reads in the data for each adduct transformation.
@:param The name of the file containing the transformation details
@:return A dictionary with key value pairs - ("transform name", list containing tranform parameters)
'''
def read_transformations(file_name):

    transform_details = []

    try:
        transform_file = open(file_name, 'r')

        line = transform_file.readline()

        while line != "":
            line = line.rstrip('\n')
            params = re.split(',', line)
            params[1:] = list(map(float, params[1:]))
            transform_details.append(params)
            line = transform_file.readline()

        transform_file.close()

    except:
        print("Error reading transformations file")

    return transform_details

'''
Reads in mass and RT data from file
@:param The file name
@:return A list containing (mass, RT) pairs
'''
def read_data(file_name):
    data = []

    try:
        data_file = open(file_name, 'r')

        line = data_file.readline()

        while line != "":
            line = line.rstrip('\n').split('\t')[0:3]
            data.append(line)
            line = data_file.readline()

        del data[0]

        for line in data:
            line[0], line[1], line[2] = float(line[0]), float(line[1]), float(line[2])

        data_file.close()

    except:
        print("Error reading data file")

    return data

'''
Function to carry out the adduct transformations
@:param A list containing the mass values, the name of the transformation (adduct name), dictionary containing the
transformation parameters
'''
def transform_data(data, transform_name, transform_details):
    transform = transform_name.strip()

    params = transform_details[transform]

    return (data-params[0])/params[1] + params[2]

print(read_transformations("mulsub2.txt"))