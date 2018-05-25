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
    data = open(filename, 'r')
    data = csv.reader(data)
    data = [row for row in data]
    for row in data:
        works = False
        try:
            value = float(row[1])
            works = True
        except:
            pass
        if works:
            if row[0] not in list(runs_data.keys()):
                runs_data[row[0]] = []
            runs_data[row[0]].append(float(row[1]))
    return runs_data

def iterate_through_runs_to_combine(runs):
    runs_data = {}
    for directory in runs:
        result_file = directory+'/results.csv'
        if os.path.exists(result_file):
            runs_data = extract_mobility(result_file, runs_data)
    return runs_data

def plot_data(data, title):
    plt.close()
    if not os.path.exists('output'):
        os.makedirs('output')
    fig, ax = plt.subplots()
    plt.errorbar(data[:,0], data[:,1], yerr=data[:,2])
    plt.ylabel(r"$\mu_0$ cm$^2$/Vs")
    plt.xlabel("Order-Semi-Disorder")
    plt.title(title)
    plt.savefig("output/{}.pdf".format(title))

def save_mean_data_to_csv(data):
    if not os.path.exists('output'):
        os.makedirs('output')

    text = ""
    for line in data:
        text +="{},{},{}\n".format(line[0], line[1], line[2])
    with open("output/data.csv", 'w') as f:
        f.writelines(text)
    print(text)

def calc_mean_and_dev(combine_list, properties = ['hole_mobility']):

    keys_list = list(combine_list.keys())

    total_data = {}
    for key, pair in combine_list.items():
        total_data[key] = iterate_through_runs_to_combine(pair)

    data_list = []
    for prop in properties:
        property_list = []
        for key, pair in total_data.items():
            name = "{}_{}".format(key, prop)
            p_data = total_data[key][prop]
            p_mean = np.mean(p_data)
            p_error = np.std(p_data)/np.sqrt(len(p_data))
            data_list.append([name, p_mean, p_error])
            try:
                property_list.append([float(key), p_mean, p_error])
            except:
                property_list.append([keys_list.index(key), p_mean, p_error])
        plot_data(np.array(property_list), prop)
    save_mean_data_to_csv(data_list)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--combine", default=None, required=False,
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
