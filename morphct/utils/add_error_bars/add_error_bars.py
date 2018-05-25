from morphct.code import helper_functions as hf
import numpy as np
import glob
import argparse
import os
import csv
import matplotlib.pyplot as plt

def split_argument_into_dictionary(argument):
    """
    Uses the string parsing on the argument to
    split the string into a dictionary.
    Requires:
        argument - string in the from of a dictionary
    Reutrns:
        combine_list - dictionary
    """
    #Create empty dictionary
    combine_list = {}
    #Remove curly brackets
    argument = argument[1:-1]
    #Remove whitespace
    argument = "".join([char for char in argument if char != " "])
    #Identify and split the lists based off of the '],' feature
    argument = argument.split("],")
    #Add the ']' back in at places it was just stripped in splitting
    argument = [item+"]" for item in argument[:-1]]+[argument[-1]]
    #Iterate through the key:pair strings
    for item in argument:
        #split based on ':' and take the zeroth element as the key name
        name = item.split(':')[0].split('\'')[1]
        #assume the 1th element is the list
        items = item.split(':')[1]
        #Remove the square brackets around string
        items = items[1:-1]
        #Create empty sublist
        sublist = []
        #Split the list based off commas
        for subitems in items.split(','):
            #Run a glob on the items in the list
            runs = glob.glob(subitems)
            #Add the items in the glob to the sublist
            for run in runs:
                sublist.append(run)
        #Make the key:pair combination of the keys and sublist
        combine_list[name] = sublist
    return combine_list

def extract_mobility(filename, runs_data):
    """
    Iterate through the csv file and create
    a dictionary in which the keys are fields in
    the csv files and the pairs are lists of
    the values from these fields.
    Requires:
        filename - string, name of the csv file
        runs_data - dictionary of lists
    """
    #Open the results.csv file
    data = open(filename, 'r')
    data = csv.reader(data)
    #Turn each row into a list
    data = [row for row in data]
    #Iterate through the rows
    for row in data:
        #Use 'works' flag to determine
        #if the field is blank or a number
        works = False
        #Try to convert the field to a number
        #If possible change 'works' flag to true
        try:
            value = float(row[1])
            works = True
        #If value can't be converted, skip it
        except:
            pass
        #If it can be converted, write data to dictionary
        if works:
            #Create the field if it doesn't exist
            if row[0] not in list(runs_data.keys()):
                runs_data[row[0]] = []
            runs_data[row[0]].append(float(row[1]))
    return runs_data

def iterate_through_runs_to_combine(runs):
    runs_data = {}
    #Iterate through each run over which to average.
    for directory in runs:
        result_file = directory+'/results.csv'
        #Make sure the results.csv file exists.
        if os.path.exists(result_file):
            runs_data = extract_mobility(result_file, runs_data)
    return runs_data

def plot_data(data, title, xlabel = "Order-Semi-Disorder"):
    """
    Creates a simple plot of the data.
    Requires:
        data - array
        title - string from the property calculated
    Optional:
        xlabel - string
    """
    plt.close()
    #Place plot in 'output' directory.

    if not os.path.exists('output'):
        os.makedirs('output')
    fig, ax = plt.subplots()

    #Plot the errors
    plt.errorbar(data[:,0], data[:,1], yerr=data[:,2])
    plt.ylabel(r"$\mu_0$ cm$^2$/Vs")
    plt.xlabel(xlabel)
    plt.title(title)
    plt.savefig("output/{}.pdf".format(title))

def save_mean_data_to_csv(data):
    """
    Write the average and standard error data
    to a file.
    Requires:
        data - list
    """
    #Puts data in directory called 'output'
    #Creates one if not present.
    if not os.path.exists('output'):
        os.makedirs('output')

    text = ""
    for line in data:
        text +="{},{},{}\n".format(line[0], line[1], line[2])

    #Will over write data.csv each time.
    with open("output/data.csv", 'w') as f:
        f.writelines(text)

def calc_mean_and_dev(combine_list, properties = ['hole_mobility']):
    """
    Iterate through the dictionary of lists to get the
    information from the appropriate csv files, then calculate 
    the average and standard error of the data.
    Plot and write this data to files.
    Requires:
        combine_list - dictionary of lists where each key
            has a list of runs over which to be averaged.
    Optional:
        properties - list of calculated properties to plot.
            defaults to the hole_mobility
    """

    #Make a list out of the dictionary keys
    #(Used for assigning numbers to strings)
    keys_list = list(combine_list.keys())

    #Create a big dictionary that will have a dictionary of
    #list for each output of the csv file.
    total_data = {}

    #Iterate through each of run set to combine
    for key, pair in combine_list.items():
        total_data[key] = iterate_through_runs_to_combine(pair)

    #Turn data into lists and arrays for writing and plotting
    #the desired properties.
    data_list = []
    #Go through the properties first for plotting each property separately.
    for prop in properties:
        property_list = []
        #Get the average and deviation from each run for the desired property.
        for key, pair in total_data.items():
            name = "{}_{}".format(key, prop)
            p_data = total_data[key][prop]
            p_mean = np.mean(p_data)
            p_error = np.std(p_data)/np.sqrt(len(p_data))
            data_list.append([name, p_mean, p_error])
            #If the keys are numbers write to list
            try:
                property_list.append([float(key), p_mean, p_error])
            #Else convert to number based on its index in the list
            except:
                property_list.append([keys_list.index(key), p_mean, p_error])
        #Plot the data
        plot_data(np.array(property_list), prop)
    #Write all the data into a single file
    save_mean_data_to_csv(data_list)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--combine", default=None, required=True,
                        type = split_argument_into_dictionary,
                        help=('''Combine mobility calculations from list of runs to create error bars.'
                            'Takes a dictionary where keys are strings and pairs are lists of directories.'
                            'Glob is applied to items within lists so that wildcards will 
                            expand all items fitting the glob into the list.
                            'For example: "{'crystal':[cryst-*, large_p3ht_ordered],}"
                            will produce: {'crystal':['cryst-frame-1', 'cryst-frame-2',
                            ...'cryst-frame-N', 'large_p3ht_ordered']
                            'Function will look for result.csv files for mobility if they exist.'
                            '''))
    args, directory_list = parser.parse_known_args()
    calc_mean_and_dev(args.combine)
