from morphct.code import helper_functions as hf
import numpy as np
import glob
import argparse

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--combine", default=None, required=False,
                        type = split_argument_into_dictionary,
                        help=('''Combine mobility calculations from list of runs to create error bars.'
                            'Takes a dictionary where keys are strings and pairs are lists.'
                            'Glob is applied to items within lists so that wildcards will 
                            expand all items fitting the glob into the list.
                            'For example: "{'crystal':[cryst-*, large_p3ht_ordered],}"
                            will produce: {'crystal':['cryst-frame-1', 'cryst-frame-2',
                            ...'cryst-frame-N', 'large_p3ht_ordered']'''))
    args, directory_list = parser.parse_known_args()
    print(args.combine)
    
