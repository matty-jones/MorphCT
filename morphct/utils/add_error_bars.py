from morphct.code import helper_functions as hf
from collections import OrderedDict
import numpy as np
import glob
import argparse
import os
import csv


plt = None


def split_argument_into_dictionary(argument):
    """
    Uses the string parsing on the argument to
    split the string into a dictionary.
    Requires:
        argument - string in the from of a dictionary
    Reutrns:
        combine_list - dictionary
    """
    # Create empty dictionary
    combine_list = OrderedDict()
    # Remove curly brackets
    argument = argument[1:-1]
    # Remove whitespace
    argument = "".join([char for char in argument if char != " "])
    # Identify and split the lists based off of the '],' feature
    argument = argument.split("],")
    # Add the ']' back in at places it was just stripped in splitting
    argument = ["".join([item, "]"]) for item in argument[:-1]] + [argument[-1]]
    # Iterate through the key:pair strings
    for item in argument:
        # split based on ':' and take the zeroth element as the key name
        name = item.split(":")[0].split("'")[1]
        # assume the 1th element is the list
        items = item.split(":")[1]
        # Remove the square brackets around string
        items = items[1:-1]
        # Create empty sublist
        sublist = []
        # Split the list based off commas
        for subitems in items.split(","):
            # Run a glob on the items in the list (remove quotes from around string)
            if ((subitems[0] == "'") and (subitems[-1] == "'")) or (
                (subitems[0] == '"') and (subitems[-1] == '"')
            ):
                subitems = subitems[1:-1]
            runs = glob.glob(subitems)
            # Add the items in the glob to the sublist
            for run in runs:
                sublist.append(run)
        # Make the key:pair combination of the keys and sublist
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
    # Open the results.csv file
    data = open(filename, "r")
    data = csv.reader(data)
    # Turn each row into a list
    data = [row for row in data]
    # Iterate through the rows
    for row in data:
        # Use 'works' flag to determine
        # if the field is blank or a number
        works = False
        # Try to convert the field to a number
        # If possible change 'works' flag to true
        try:
            value = float(row[1])
            works = True
        # If value can't be converted, skip it
        except:
            pass
        # If it can be converted, write data to dictionary
        if works:
            # Create the field if it doesn't exist
            if row[0] not in list(runs_data.keys()):
                runs_data[row[0]] = []
            runs_data[row[0]].append(float(row[1]))
    return runs_data


def iterate_through_runs_to_combine(runs):
    runs_data = {}
    # Iterate through each run over which to average.
    for directory in runs:
        result_file = os.path.join(directory, "results.csv")
        # Make sure the results.csv file exists.
        if os.path.exists(result_file):
            runs_data = extract_mobility(result_file, runs_data)
    return runs_data


def plot_data(data, title, output_file, xlabel="Order-Semi-Disorder"):
    """
    Creates a simple plot of the data.
    Requires:
        data - array
        title - string from the property calculated
    Optional:
        xlabel - string
    """
    plt.close()
    # Place plot in 'output' directory.

    if not os.path.exists("output"):
        os.makedirs("output")
    fig, ax = plt.subplots()

    x_data, y_data, y_error = zip(*sorted(zip(data[:, 0], data[:, 1], data[:, 2])))
    # Plot the errors
    plt.errorbar(x_data, y_data, yerr=y_error, c="r")#, fmt='o')
    plt.yscale("log")
    # plt.yticks([2E-2, 3E-2, 4E-2, 5E-2, 6E-2, 7E-2, 8E-2, 9E-2, 1E-1], [r'$2\times10^{-2}$','',r'$4\times10^{-2}$','','','','',
    #    '',r'$1\times10^{-1}$'])
    plt.ylabel(r"Mobility (cm$^{2}$V$^{-1}$s$^{-1}$)")
    # plt.ylabel(r"Anisotropy (Arb. U.)")
    plt.xlabel(xlabel)
    # plt.title(title)
    (save_dir, save_file) = os.path.split(output_file)
    os.makedirs(save_dir, exist_ok=True)
    plt.savefig(os.path.join(save_dir, save_file))
    print("Figure saved as", output_file)


def save_mean_data_to_csv(data, prop):
    """
    Write the average and standard error data
    to a file.
    Requires:
        data - list
    """
    # Puts data in directory called 'output'
    # Creates one if not present.
    if not os.path.exists("output"):
        os.makedirs("output")

    text = ""
    for line in data:
        text = "".join(
            [
                text,
                "{0:s},{1:s},{2:s}\n".format(
                    *tuple(map(str, [line[0], line[1], line[2]]))
                ),
            ]
        )

    # Will over write data.csv each time.
    with open(os.path.join("output", "".join([prop, ".csv"])), "w+") as f:
        f.writelines(text)


def calc_mean_and_dev(
    combine_list, sequence, x_label, output_file, prop="hole_mobility",
    cutoff_prop=None, cutoff_val=0
):
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

    # Make a list out of the dictionary keys
    # (Used for assigning numbers to strings)
    keys_list = list(combine_list.keys())

    # Create a big dictionary that will have a dictionary of
    # list for each output of the csv file.
    total_data = {}

    # Iterate through each of run set to combine
    for key, pair in combine_list.items():
        total_data[key] = iterate_through_runs_to_combine(pair)
    # Turn data into lists and arrays for writing and plotting
    # the desired properties.
    data_list = []
    property_list = []
    # Get the average and deviation from each run for the desired property.
    for key, pair in total_data.items():
        name = "{}_{}".format(key, prop)
        # Check the cutoff in each case
        if cutoff_prop is not None:
            p_data = []
            for val_index, cutoff_check in enumerate(total_data[key][cutoff_prop]):
                if float(cutoff_check) >= float(cutoff_val):
                    p_data.append(total_data[key][prop][val_index])
                else:
                    print(
                        "Skipping", prop, "result from", key, val_index, "as",
                        cutoff_prop, "value =", cutoff_check, "which is <", cutoff_val
                    )
        else:
            p_data = total_data[key][prop]
        p_mean = np.mean(p_data)
        p_error = np.std(p_data) / np.sqrt(len(p_data))
        data_list.append([name, p_mean, p_error])
        # If the keys are numbers write to list
        try:
            property_list.append([float(key), p_mean, p_error])
        # Else convert to number based on its index in the list
        except:
            property_list.append([sequence[keys_list.index(key)], p_mean, p_error])
    # Plot the data
    plot_data(np.array(property_list), prop, output_file, xlabel=x_label)
    # Write all the data into a single file
    save_mean_data_to_csv(data_list, prop)


def main():
    global plt
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--combine",
        required=True,
        type=split_argument_into_dictionary,
        help=(
            """Combine mobility calculations from list of runs to create error bars.
                              Takes a dictionary where keys are strings and pairs are lists of directories.
                              Glob is applied to items within lists so that wildcards will
                              expand all items fitting the glob into the list.
                              For example: "{'crystal':[cryst-*, large_p3ht_ordered],}
                              will produce: {'crystal':['cryst-frame-1', 'cryst-frame-2',
                              ...'cryst-frame-N', 'large_p3ht_ordered']
                              Function will look for result.csv files for mobility if they exist.
                              """
        ),
    )
    parser.add_argument(
        "-s",
        "--sequence",
        required=True,
        type=lambda s: [float(item) for item in s.split(",")],
        default=None,
        help=(
            """Create a figure in the current directory that describes the evolution
                              of the anisotropy/mobility using the specified comma-delimited string
                              as the sequence of x values. For instance -s '1.5,1.75,2.0,2.25,2.5'
                              will assign each of the 5 keys in the args.combine list to these x-values when
                              plotting the errorbar graph."""
        ),
    )
    parser.add_argument(
        "-x", "--x_label", required=True, help=("""Set an x label for the final plot""")
    )
    parser.add_argument(
        "-p",
        "--prop",
        required=True,
        help=(
            """The property, as specified in the results.csv, to use during the aggregation
                              of the data and the error bar calculation."""
        ),
    )
    parser.add_argument(
        "-cp",
        "--cutoff_prop",
        required=False,
        default=None,
        help=(
            """The cutoff property to check. If a system's cutoff_prop is not >=
            the cutoff_val, then it is skipped from the mean and std calculation."""
        ),
    )
    parser.add_argument(
        "-cv",
        "--cutoff_val",
        required=False,
        default=0,
        help=(
            """The value to check against for the cutoff property. If a system's
            cutoff_prop is not >= the cutoff_val, then it is skipped from the mean and
            std calculation."""
        ),
    )
    parser.add_argument(
        "-b",
        "--backend",
        default=None,
        required=False,
        help=(
            "Specify a backend for matplotlib to use when plotting. Default = user defined in the"
            " .matplotlibrc."
        ),
    )
    parser.add_argument(
        "-o",
        "--output_file",
        required=False,
        default="./output/hole_mobility.pdf",
        help=("""Set an x label for the final plot"""),
    )
    args, directory_list = parser.parse_known_args()
    if args.backend is not None:
        import matplotlib

        matplotlib.use(args.backend.strip())
    import matplotlib.pyplot as plt

    print("Input list to combine:")
    for key, val in args.combine.items():
        print(key, val)

    calc_mean_and_dev(
        args.combine, args.sequence, args.x_label, args.output_file, prop=args.prop,
        cutoff_prop=args.cutoff_prop, cutoff_val=args.cutoff_val,
    )


if __name__ == "__main__":
    main()
