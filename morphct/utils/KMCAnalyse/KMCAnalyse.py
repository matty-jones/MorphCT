import os
import sys
import pickle
import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import scipy.stats
from morphct.code import helperFunctions as hf
from collections import OrderedDict
try:
    import mpl_toolkits.mplot3d as p3
except ImportError:
    print("Could not import 3D plotting engine, calling the plotMolecule3D function will result in an error!")
import shutil
import glob
import re
import argparse
import copy

elementary_charge = 1.60217657E-19  # C
kB = 1.3806488E-23  # m^{2} kg s^{-2} K^{-1}
hbar = 1.05457173E-34  # m^{2} kg s^{-1}
temperature = 290  # K


def load_KMC_results_pickle(directory):
    try:
        with open(directory + '/KMC/KMC_results.pickle', 'rb') as pickle_file:
            carrier_data = pickle.load(pickle_file)
    except FileNotFoundError:
        print("No final KMC_results.pickle found. Creating it from incomplete parts...")
        create_results_pickle(directory)
        with open(directory + '/KMC/KMC_results.pickle', 'rb') as pickle_file:
            carrier_data = pickle.load(pickle_file)
    except UnicodeDecodeError:
        with open(directory + '/KMC/KMC_results.pickle', 'rb') as pickle_file:
            carrier_data = pickle.load(pickle_file, encoding='latin1')
    return carrier_data


def split_carriers_by_type(carrier_data):
    # If only one carrier type has been given, call the carriers holes and skip the electron calculations
    list_variables = ['current_time', 'ID', 'no_hops', 'displacement', 'lifetime', 'final_position', 'image',
                      'initial_position']
    try:
        carrier_data_holes = {'carrier_history_matrix': carrier_data['hole_history_matrix'],
                              'seed': carrier_data['seed']}
        carrier_data_electrons = {'carrier_history_matrix': carrier_data['electron_history_matrix'],
                                  'seed': carrier_data['seed']}
        for list_var in list_variables:
            carrier_data_holes[list_var] = []
            carrier_data_electrons[list_var] = []
            for carrier_index, charge_type in enumerate(carrier_data['carrier_type']):
                if charge_type == 'hole':
                    carrier_data_holes[list_var].append(carrier_data[list_var][carrier_index])
                elif charge_type == 'electron':
                    carrier_data_electrons[list_var].append(carrier_data[list_var][carrier_index])
    except:
        print("Multiple charge carriers not found, assuming donor material and holes only")
        try:
            carrier_data_holes = {'carrier_history_matrix': carrier_data['carrier_history_matrix'],
                                  'seed': carrier_data['seed']}
        except KeyError:
            carrier_data_holes = {'carrier_history_matrix': carrier_data['carrier_history_matrix'], 'seed': 0}
        carrier_data_electrons = None
        for list_var in list_variables:
            carrier_data_holes[list_var] = []
            for carrier_index, carrier_ID in enumerate(carrier_data['ID']):
                carrier_data_holes[list_var].append(carrier_data[list_var][carrier_index])
    return carrier_data_holes, carrier_data_electrons


def get_carrier_data(carrier_data):
    try:
        carrier_history = carrier_data['carrier_history_matrix']
    except:
        carrier_history = None
    total_data_points = 0
    total_data_points_averaged_over = 0
    squared_disps = {}
    actual_times = {}
    for carrier_index, displacement in enumerate(carrier_data['displacement']):
        if (carrier_data['current_time'][carrier_index] > carrier_data['lifetime'][carrier_index] * 2) or\
           (carrier_data['current_time'][carrier_index] < carrier_data['lifetime'][carrier_index] / 2.0) or\
           (carrier_data['no_hops'][carrier_index] == 1):
            total_data_points += 1
            continue
        carrier_key = str(carrier_data['lifetime'][carrier_index])
        if carrier_key not in squared_disps:
            squared_disps[carrier_key] = [(carrier_data['displacement'][carrier_index] * 1E-10) ** 2]  # A -> m
            actual_times[carrier_key] = [carrier_data['current_time'][carrier_index]]
        else:
            squared_disps[carrier_key].append((carrier_data['displacement'][carrier_index] * 1E-10) ** 2)  # A -> m
            actual_times[carrier_key].append(carrier_data['current_time'][carrier_index])
        # Also keep track of whether each carrier is a hole or an electron
        total_data_points_averaged_over += 1
        total_data_points += 1
    times = []
    MSDs = []
    time_standard_errors = []
    MSD_standard_errors = []
    for time, disps in squared_disps.items():
        times.append(float(time))
        time_standard_errors.append(np.std(actual_times[time]) / len(actual_times[time]))
        MSDs.append(np.average(disps))
        MSD_standard_errors.append(np.std(disps) / len(disps))
    return carrier_history, times, MSDs, time_standard_errors, MSD_standard_errors


def create_array_for_plot_connections(chromophore_list, carrier_history, sim_dims):
    """
    Function to create an array of with a starting point, a vector
    and the number of hops that occured.
    Requires:
        chromophore_list,
        carrier_history
        sim_dims
    Returns:
        7xN array
    """
    # Create an "empty" array to store data.
    connections_array = np.zeros(7)
    # Iterate through the chromophoreList
    for i, chromo in enumerate(chromophore_list):
        # Iterate through the neighbors of the chromophore
        for neighbor in zip(chromo.neighbours):
            index = neighbor[0][0]  # index of the neighbor
            image = neighbor[0][1]  # check to see if they are in the same relative image
            # Only consider one direction.
            if i < index:
                # Get the vector between the two chromophores.
                if not np.count_nonzero(image):
                    vector = chromophore_list[index].posn - chromo.posn
                # Account for periodic boundary conditions if not in same relative image.
                else:
                    vector = chromophore_list[index].posn - chromo.posn
                    vector += image * np.array([2 * sim_dims[0][1], 2 * sim_dims[1][1], 2 * sim_dims[2][1]])

                # Get the net number of times the path was travelled.
                forward = carrier_history[index, i]
                reverse = carrier_history[i, index]
                times_travelled = abs(forward - reverse)

                # Append the array if the net times travelled is greater than 0
                if times_travelled > 0:
                    datum = np.hstack((chromo.posn, vector, np.array([np.log10(times_travelled)])))
                    connections_array = np.vstack((connections_array, datum))
    return connections_array[1:]  # Return the array excluding the zeros first line.


def plot_connections(chromophore_list, sim_dims, carrier_history, directory, carrier_type):
    # A complicated function that shows connections between carriers in 3D that carriers prefer to hop between.
    # Connections that are frequently used are highlighted in black, whereas rarely used connections are more white.
    # Import matplotlib color modules to set up color bar.
    import matplotlib.colors
    import matplotlib.cm as cmx

    # Create a figure class
    fig = plt.figure(figsize=(7, 6))
    # Make a 3D subplot
    ax = fig.add_subplot(111, projection='3d')

    # Create the array for all the chromophore connections
    connections_array = create_array_for_plot_connections(chromophore_list, carrier_history, sim_dims)

    # Determine the smalles, non-zero number of times two chromophores are connected.
    vmin = np.min(np.array(connections_array)[:, 6])
    # Determine the max number of times two chormophores are connected.
    vmax = np.max(np.array(connections_array)[:, 6])

    # Set up the color bar.
    plasma = plt.get_cmap('plasma')
    c_norm = matplotlib.colors.normalize(vmin=vmin, vmax=vmax)
    scalar_map = cmx.scalar_mappable(norm=c_norm, cmap=plasma)
    hopcolors = scalar_map.to_rgba(connections_array[:, 6])

    # Set up the intensity for the hops so more travelled paths are more intense
    alphas = connections_array[:, 6] / vmax
    hopcolors[:, 3] = alphas

    # Plot the vectors between two chromophores
    ax.quiver(connections_array[:, 0],
              connections_array[:, 1],
              connections_array[:, 2],
              connections_array[:, 3],
              connections_array[:, 4],
              connections_array[:, 5],
              color=hopcolors,
              arrow_length_ratio=0, linewidth=0.7)

    # Plot the color bar
    scalar_map.set_array(connections_array[:, 6])
    tick_location = np.arange(0, int(vmax) + 1, 1)
    fig.colorbar(scalar_map, ticks=tick_location, shrink=0.8, aspect=20)

    # Draw boxlines
    # Varying X
    ax.plot([sim_dims[0][0], sim_dims[0][1]], [sim_dims[1][0], sim_dims[1][0]], [sim_dims[2][0], sim_dims[2][0]],
            c='k', linewidth=1.0)
    ax.plot([sim_dims[0][0], sim_dims[0][1]], [sim_dims[1][1], sim_dims[1][1]], [sim_dims[2][0], sim_dims[2][0]],
            c='k', linewidth=1.0)
    ax.plot([sim_dims[0][0], sim_dims[0][1]], [sim_dims[1][0], sim_dims[1][0]], [sim_dims[2][1], sim_dims[2][1]],
            c='k', linewidth=1.0)
    ax.plot([sim_dims[0][0], sim_dims[0][1]], [sim_dims[1][1], sim_dims[1][1]], [sim_dims[2][1], sim_dims[2][1]],
            c='k', linewidth=1.0)
    # Varying Y
    ax.plot([sim_dims[0][0], sim_dims[0][0]], [sim_dims[1][0], sim_dims[1][1]], [sim_dims[2][0], sim_dims[2][0]],
            c='k', linewidth=1.0)
    ax.plot([sim_dims[0][1], sim_dims[0][1]], [sim_dims[1][0], sim_dims[1][1]], [sim_dims[2][0], sim_dims[2][0]],
            c='k', linewidth=1.0)
    ax.plot([sim_dims[0][0], sim_dims[0][0]], [sim_dims[1][0], sim_dims[1][1]], [sim_dims[2][1], sim_dims[2][1]],
            c='k', linewidth=1.0)
    ax.plot([sim_dims[0][1], sim_dims[0][1]], [sim_dims[1][0], sim_dims[1][1]], [sim_dims[2][1], sim_dims[2][1]],
            c='k', linewidth=1.0)
    # Varying Z
    ax.plot([sim_dims[0][0], sim_dims[0][0]], [sim_dims[1][0], sim_dims[1][0]], [sim_dims[2][0], sim_dims[2][1]],
            c='k', linewidth=1.0)
    ax.plot([sim_dims[0][0], sim_dims[0][0]], [sim_dims[1][1], sim_dims[1][1]], [sim_dims[2][0], sim_dims[2][1]],
            c='k', linewidth=1.0)
    ax.plot([sim_dims[0][1], sim_dims[0][1]], [sim_dims[1][0], sim_dims[1][0]], [sim_dims[2][0], sim_dims[2][1]],
            c='k', linewidth=1.0)
    ax.plot([sim_dims[0][1], sim_dims[0][1]], [sim_dims[1][1], sim_dims[1][1]], [sim_dims[2][0], sim_dims[2][1]],
            c='k', linewidth=1.0)

    # Name and save the figure.
    plt.title(carrier_type.capitalize() + ' Network', y=1.1)
    file_name = '01_3d' + carrier_type + '.pdf'
    plt.savefig(directory + '/figures/' + file_name, bbox_inches='tight')
    print("Figure saved as", directory + "/figures/" + file_name)
    plt.clf()


def calc_mobility(lin_fit_X, lin_fit_Y, av_time_error, av_MSD_error):
    # YVals have a std error avMSDError associated with them
    # XVals have a std error avTimeError assosciated with them
    numerator = lin_fit_Y[-1] - lin_fit_Y[0]
    denominator = lin_fit_X[-1] - lin_fit_X[0]
    diffusion_coeff = numerator / denominator
    # The error in the mobility is the proportionally the same as the error in the diffusion coefficient as the
    # other variables are constants with zero error
    diff_error = diffusion_coeff * np.sqrt((av_MSD_error / numerator)**2 + (av_time_error / denominator)**2)
    # Use Einstein relation (include the factor of 1/6!! It is in the Carbone/Troisi 2014 paper)
    mobility = elementary_charge * diffusion_coeff / (6 * kB * temperature)  # This is in m^{2} / vs
    # Convert to cm^{2}/ Vs
    mobility *= (100**2)
    mob_error = (diff_error / diffusion_coeff) * mobility
    return mobility, mob_error


def plot_MSD(times, MSDs, time_standard_errors, MSD_standard_errors, directory, carrier_type):
    fit_X = np.linspace(np.min(times), np.max(times), 100)
    gradient, intercept, r_val, p_val, std_err = scipy.stats.linregress(times, MSDs)
    print("Standard Error", std_err)
    print("Fitting r_val =", r_val)
    fit_Y = (fit_X * gradient) + intercept
    mobility, mob_error = calc_mobility(fit_X, fit_Y, np.average(time_standard_errors),
                                        np.average(MSD_standard_errors))
    plt.plot(times, MSDs)
    plt.errorbar(times, MSDs, xerr=time_standard_errors, yerr=MSD_standard_errors)
    plt.plot(fit_X, fit_Y, 'r')
    plt.xlabel('time (s)')
    plt.ylabel('MSD (m' + r'$^{2}$)')
    mobility_string = '%.3e' % mobility
    plt.title(r'$\mu_{0,' + carrier_type[0] + r'}$' + ' = ' + mobility_string + ' cm' + r'$^{2}$/vs' % (mobility),
              y=1.1)
    file_name = '18_lin_MSD' + carrier_type + '.pdf'
    plt.savefig(directory + '/figures/' + file_name)
    plt.clf()
    print("Figure saved as", directory + "/figures/" + file_name)
    plt.semilogx(times, MSDs)
    plt.errorbar(times, MSDs, xerr=time_standard_errors, yerr=MSD_standard_errors)
    plt.semilogx(fit_X, fit_Y, 'r')
    plt.xlabel('time (s)')
    plt.ylabel('MSD (m' + r'$^{2}$)')
    mobility_string = '%.3e' % mobility
    plt.title(r'$\mu_{0,' + carrier_type[0] + r'}$' + ' = ' + mobility_string + ' cm' + r'$^{2}$/vs' % (mobility),
              y=1.1)
    file_name = '19_semi_log_MSD' + carrier_type + '.pdf'
    plt.savefig(directory + '/figures/' + file_name)
    plt.clf()
    print("Figure saved as", directory + "/figures/" + file_name)
    plt.plot(times, MSDs)
    plt.errorbar(times, MSDs, xerr=time_standard_errors, yerr=MSD_standard_errors)
    plt.plot(fit_X, fit_Y, 'r')
    plt.xlabel('time (s)')
    plt.ylabel('MSD (m' + r'$^{2}$)')
    plt.xscale('log')
    plt.yscale('log')
    mobility_string = '%.3e' % mobility
    plt.title(r'$\mu_{0,' + carrier_type[0] + r'}$' + ' = ' + mobility_string + ' cm' + r'$^{2}$/vs' % (mobility),
              y=1.1)
    file_name = '20_log_MSD' + carrier_type + '.pdf'
    plt.savefig(directory + '/figures/' + file_name)
    plt.clf()
    print("Figure saved as", directory + "/figures/" + file_name)
    return mobility, mob_error, r_val**2


def calculate_anisotropy(xvals, yvals, zvals):
    # First calculate the `centre of position' for the particles
    centre = [np.mean(xvals), np.mean(yvals), np.mean(zvals)]
    # First calculate the gyration tensor:
    sxx = 0
    sxy = 0
    sxz = 0
    syy = 0
    syz = 0
    szz = 0
    for carrier_ID, raw_xval in enumerate(xvals):
        xval = raw_xval - centre[0]
        yval = yvals[carrier_ID] - centre[1]
        zval = zvals[carrier_ID] - centre[2]
        sxx += xval * xval
        sxy += xval * yval
        sxz += xval * zval
        syy += yval * yval
        syz += yval * zval
        szz += zval * zval
    S = np.array([[sxx, sxy, sxz], [sxy, syy, syz], [sxz, syz, szz]])
    eigenvalues, eigenvectors = np.linalg.eig(S)
    # Diagonalisation of S is the diagonal matrix of the eigenvalues in ascending order
    # diagonalMatrix = np.diag(sorted(eigenValues))
    # We only need the eigenvalues though, no more matrix multiplication
    diagonal = sorted(eigenvalues)
    # Then calculate the relative shape anisotropy (kappa**2)
    anisotropy = (3 / 2) * (((diagonal[0] ** 2) + (diagonal[1] ** 2) + (diagonal[2] ** 2))
                            / ((diagonal[0] + diagonal[1] + diagonal[2]) ** 2)) - (1 / 2)
    return anisotropy


def plot_anisotropy(carrier_data, directory, sim_dims, carrier_type, plot3D_graphs):
    sim_extent = [value[1] - value[0] for value in sim_dims]
    xvals = []
    yvals = []
    zvals = []
    colours = []
    sim_dims_nm = list(map(list, np.array(sim_dims) / 10.))
    # Get the indices of the carriers that travelled the furthest
    if len(carrier_data['final_position']) <= 1000:
        carrier_indices_to_use = range(len(carrier_data['final_position']))
    else:
        displacements = copy.deepcopy(np.array(carrier_data['displacement']))
        carrier_indices_to_use = displacements.argsort()[-1000:][::-1]
    for carrier_no in carrier_indices_to_use:
        posn = carrier_data['final_position'][carrier_no]
        # if bool(sum([x < -3 or x > 3 for x in image])):
        #     continue
        position = [0.0, 0.0, 0.0]
        for axis in range(len(posn)):
            position[axis] = (carrier_data['image'][carrier_no][axis] * sim_extent[axis]) + posn[axis]
        xvals.append(position[0] / 10.)
        yvals.append(position[1] / 10.)
        zvals.append(position[2] / 10.)
        colours.append('b')
    anisotropy = calculate_anisotropy(xvals, yvals, zvals)
    if not plot3D_graphs:
        return anisotropy
    print("----------====================----------")
    print(carrier_type + " charge transport anisotropy calculated as", anisotropy)
    print("----------====================----------")
    # Reduce number of plot markers
    fig = plt.gcf()
    ax = p3.Axes3D(fig)
    if len(xvals) > 1000:
        xvals = xvals[0:len(xvals):len(xvals) // 1000]
        yvals = yvals[0:len(yvals):len(yvals) // 1000]
        zvals = zvals[0:len(zvals):len(zvals) // 1000]
    plt.scatter(xvals, yvals, zs=zvals, c=colours, s=20)
    plt.scatter(0, 0, zs=0, c='r', s=50)
    # Draw boxlines
    # Varying X
    ax.plot([sim_dims_nm[0][0], sim_dims_nm[0][1]], [sim_dims_nm[1][0], sim_dims_nm[1][0]],
            [sim_dims_nm[2][0], sim_dims_nm[2][0]], c='k', linewidth=1.0)
    ax.plot([sim_dims_nm[0][0], sim_dims_nm[0][1]], [sim_dims_nm[1][1], sim_dims_nm[1][1]],
            [sim_dims_nm[2][0], sim_dims_nm[2][0]], c='k', linewidth=1.0)
    ax.plot([sim_dims_nm[0][0], sim_dims_nm[0][1]], [sim_dims_nm[1][0], sim_dims_nm[1][0]],
            [sim_dims_nm[2][1], sim_dims_nm[2][1]], c='k', linewidth=1.0)
    ax.plot([sim_dims_nm[0][0], sim_dims_nm[0][1]], [sim_dims_nm[1][1], sim_dims_nm[1][1]],
            [sim_dims_nm[2][1], sim_dims_nm[2][1]], c='k', linewidth=1.0)
    # Varying Y
    ax.plot([sim_dims_nm[0][0], sim_dims_nm[0][0]], [sim_dims_nm[1][0], sim_dims_nm[1][1]],
            [sim_dims_nm[2][0], sim_dims_nm[2][0]], c='k', linewidth=1.0)
    ax.plot([sim_dims_nm[0][1], sim_dims_nm[0][1]], [sim_dims_nm[1][0], sim_dims_nm[1][1]],
            [sim_dims_nm[2][0], sim_dims_nm[2][0]], c='k', linewidth=1.0)
    ax.plot([sim_dims_nm[0][0], sim_dims_nm[0][0]], [sim_dims_nm[1][0], sim_dims_nm[1][1]],
            [sim_dims_nm[2][1], sim_dims_nm[2][1]], c='k', linewidth=1.0)
    ax.plot([sim_dims_nm[0][1], sim_dims_nm[0][1]], [sim_dims_nm[1][0], sim_dims_nm[1][1]],
            [sim_dims_nm[2][1], sim_dims_nm[2][1]], c='k', linewidth=1.0)
    # Varying Z
    ax.plot([sim_dims_nm[0][0], sim_dims_nm[0][0]], [sim_dims_nm[1][0], sim_dims_nm[1][0]],
            [sim_dims_nm[2][0], sim_dims_nm[2][1]], c='k', linewidth=1.0)
    ax.plot([sim_dims_nm[0][0], sim_dims_nm[0][0]], [sim_dims_nm[1][1], sim_dims_nm[1][1]],
            [sim_dims_nm[2][0], sim_dims_nm[2][1]], c='k', linewidth=1.0)
    ax.plot([sim_dims_nm[0][1], sim_dims_nm[0][1]], [sim_dims_nm[1][0], sim_dims_nm[1][0]],
            [sim_dims_nm[2][0], sim_dims_nm[2][1]], c='k', linewidth=1.0)
    ax.plot([sim_dims_nm[0][1], sim_dims_nm[0][1]], [sim_dims_nm[1][1], sim_dims_nm[1][1]],
            [sim_dims_nm[2][0], sim_dims_nm[2][1]], c='k', linewidth=1.0)
    ax.set_xlabel('X (nm)', fontsize=20, labelpad=40)
    ax.set_ylabel('Y (nm)', fontsize=20, labelpad=40)
    ax.set_zlabel('Z (nm)', fontsize=20, labelpad=40)
    maximum = max([max(xvals), max(yvals), max(zvals)])
    ax.set_xlim([-maximum, maximum])
    ax.set_ylim([-maximum, maximum])
    ax.set_zlim([-maximum, maximum])
    for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks() + ax.zaxis.get_major_ticks():
        tick.label.set_fontsize(16)
    ax.dist = 11
    if carrier_type == 'hole':
        figure_index = '08'
    elif carrier_type == 'electron':
        figure_index = '09'
    plt.title('anisotropy (' + carrier_type + ')', y=1.1)
    file_name = directory + '/figures/' + figure_index + '_anisotropy' + carrier_type + '.pdf'
    plt.savefig(file_name, bbox_inches='tight')
    plt.clf()
    print("Figure saved as", file_name)
    return anisotropy


def get_temp_val(string):
    hyphen_list = hf.find_index(string, '-')
    temp_val = float(string[hyphen_list[-2] + 2: hyphen_list[-1]])
    return temp_val


def get_frame_val(string):
    hyphen_list = hf.find_index(string, '-')
    temp_val = int(string[hyphen_list[0] + 1:hyphen_list[1]])
    return temp_val


def plot_temperature_progression(temp_data, mobility_data, anisotropy_data, carrier_type, x_label):
    plt.gcf()
    xvals = temp_data
    yvals = list(np.array(mobility_data)[:, 0])
    yerrs = list(np.array(mobility_data)[:, 1])
    plt.xlabel(x_label)
    plt.ylabel('mobility (cm' + r'$^{2}$ ' + 'V' + r'$^{-1}$' + r's$^{-1}$)')
    plt.semilogy(xvals, yvals, c='k')
    plt.errorbar(xvals, yvals, xerr=0, yerr=yerrs)
    file_name = './mobility' + carrier_type + '.pdf'
    plt.savefig(file_name)
    plt.clf()
    print("Figure saved as " + file_name)

    plt.plot(temp_data, anisotropy_data, c='r')
    file_name = './anisotropy' + carrier_type + '.pdf'
    plt.xlabel(x_label)
    plt.ylabel(r'$\kappa$' + ' (Arb. U)')
    plt.savefig(file_name)
    plt.clf()
    print("Figure saved as " + file_name)


def calculate_lambdaij(chromo_length):
    # The equation for the internal reorganisation energy was obtained from the data given in
    # Johansson, E and Larsson, S; 2004, Synthetic Metals 144: 183-191.
    # External reorganisation energy obtained from
    # Liu, T and Cheung, D. L. and Troisi, A; 2011, Phys. Chem. Chem. Phys. 13: 21461-21470
    lambda_external = 0.11  # eV
    if chromo_length < 12:
        lambda_internal = 0.20826 - (chromo_length * 0.01196)
    else:
        lambda_internal = 0.06474
    lambdae_V = lambda_external + lambda_internal
    return lambdae_V


def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))


def gauss_fit(data):
    mean = np.mean(data)
    std = np.std(data)
    hist, bin_edges = np.histogram(data, bins=100)
    try:
        fit_args, fit_conv = curve_fit(gaussian, bin_edges[:-1], hist, p0=[1, mean, std])
    except RuntimeError:
        fit_args = None
    return bin_edges, fit_args, mean, std


def split_molecules(input_dictionary):
    # Split the full morphology into individual molecules
    molecule_AAIDs = []
    molecule_lengths = []
    # Create a lookup table `neighbour list' for all connected atoms called {bondedAtoms}
    bonded_atoms = hf.obtain_bonded_list(input_dictionary['bond'])
    molecule_list = [i for i in range(len(input_dictionary['type']))]
    # Recursively add all atoms in the neighbour list to this molecule
    for mol_ID in range(len(molecule_list)):
        molecule_list = update_molecule(mol_ID, molecule_list, bonded_atoms)
    # Create a dictionary of the molecule data
    molecule_data = {}
    for atom_ID in range(len(input_dictionary['type'])):
        if molecule_list[atom_ID] not in molecule_data:
            molecule_data[molecule_list[atom_ID]] = [atom_ID]
        else:
            molecule_data[molecule_list[atom_ID]].append(atom_ID)
    # Return the list of AAIDs and the lengths of the molecules
    for molecule_ID in list(molecule_data.keys()):
        molecule_AAIDs.append(sorted(molecule_data[molecule_ID]))
        molecule_lengths.append(len(molecule_data[molecule_ID]))
    return molecule_AAIDs, molecule_lengths


def update_molecule(atom_ID, molecule_list, bonded_atoms):
    # Recursively add all neighbours of atom number atomID to this molecule
    try:
        for bonded_atom in bonded_atoms[atom_ID]:
            # If the moleculeID of the bonded atom is larger than that of the current one,
            # update the bonded atom's ID to the current one's to put it in this molecule,
            # then iterate through all of the bonded atom's neighbours
            if molecule_list[bonded_atom] > molecule_list[atom_ID]:
                molecule_list[bonded_atom] = molecule_list[atom_ID]
                molecule_list = update_molecule(bonded_atom, molecule_list, bonded_atoms)
            # If the moleculeID of the current atom is larger than that of the bonded one,
            # update the current atom's ID to the bonded one's to put it in this molecule,
            # then iterate through all of the current atom's neighbours
            elif molecule_list[bonded_atom] < molecule_list[atom_ID]:
                molecule_list[atom_ID] = molecule_list[bonded_atom]
                molecule_list = update_molecule(atom_ID, molecule_list, bonded_atoms)
            # Else: both the current and the bonded atom are already known to be in this
            # molecule, so we don't have to do anything else.
    except KeyError:
        # This means that there are no bonded CG sites (i.e. it's a single molecule)
        pass
    return molecule_list


def get_neighbour_cut_off(chromophore_list, morphology_shape, output_dir, periodic=True, specified_cut_off_donor=None,
                          specified_cut_off_acceptor=None):
    specified_cut_offs = {'Donor': specified_cut_off_donor, 'Acceptor': specified_cut_off_acceptor}
    separation_dist_donor = []
    separation_dist_acceptor = []
    for chromo1 in chromophore_list:
        for chromo2_details in chromo1.neighbours:
            if (chromo2_details is None)\
               or ((periodic is False) and (not np.array_equal(chromo2_details[1], [0, 0, 0])))\
               or (chromo1.ID == chromophore_list[chromo2_details[0]].ID):
                continue
            chromo2 = chromophore_list[chromo2_details[0]]
            separation = np.linalg.norm((np.array(chromo2.posn) + (np.array(chromo2_details[1])
                                                                   * np.array(morphology_shape))) - chromo1.posn)
            if chromo1.species == 'Donor':
                separation_dist_donor.append(separation)
            elif chromo1.species == 'Acceptor':
                separation_dist_acceptor.append(separation)
    cut_offs = []
    material = ['Donor', 'Acceptor']
    for material_type, separation_dist in enumerate([separation_dist_donor, separation_dist_acceptor]):
        if len(separation_dist) == 0:
            cut_offs.append(None)
            continue
        plt.figure()
        (n, bin_edges, patches) = plt.hist(separation_dist, bins=20, color='b')
        bins = 0.5 * (bin_edges[1:] + bin_edges[:-1])
        bins = np.insert(bins, 0, 0)
        n = np.insert(n, 0, 0)
        dn = np.diff(n)
        minima_indices = []
        maxima_indices = []
        previous_value = 1E99
        for index, val in enumerate(dn):
            if (previous_value <= 0) and (val > 0):
                minima_indices.append(index)
            if (previous_value >= 0) and (val < 0):
                maxima_indices.append(index)
            previous_value = val
        # Minimum is half way between the first maximum and the first minimum of the distribution
        if specified_cut_offs[material[material_type]] is None:
            cut_off = (bins[maxima_indices[0]] + bins[minima_indices[0]]) / 2.0
        else:
            cut_off = specified_cut_offs[material[material_type]]
        cut_offs.append(cut_off)
        plt.axvline(x=cut_off, ls='dashed', c='k')
        plt.xlabel(material[material_type] + r' r$_{ij}$' + ' (A)')
        plt.ylabel("Frequency (Arb. U.)")
        plt.savefig(output_dir + "/03_neighbour_hist_" + material[material_type].lower() + ".pdf")
        plt.close()
        print("Neighbour histogram figure saved as", output_dir + "/03_neighbourHist_"
              + material[material_type].lower() + ".pdf")
    return cut_offs


def get_stacks(chromophore_list, morphology_shape, cut_off_donor, cut_off_acceptor, periodic=True):
    cut_offs = [cut_off_donor, cut_off_acceptor]
    materials_to_check = ['Donor', 'Acceptor']
    stack_dicts = []
    for type_index, material_type in enumerate(materials_to_check):
        cut_off = cut_offs[type_index]
        if cut_off is None:
            stack_dicts.append(None)
            continue
        chromo_IDs = [chromo.ID for chromo in chromophore_list if chromo.species == material_type]
        # Create a neighbourlist based on the cutoff
        neighbour_dict = create_neighbour_list(chromophore_list, morphology_shape, cut_off, periodic, material_type)
        # Do the usual stackList neighbourList stuff
        stack_list = [_ for _ in range(len(chromophore_list))]
        for stack_ID in range(len(stack_list)):
            stack_list = update_stack(stack_ID, stack_list, neighbour_dict)
        actual_stacks = [stack_list[x] for x in chromo_IDs]
        print("There are", len(set(actual_stacks)), material_type, "stacks in the system")
        stack_dict = {}
        for index, chromophore in enumerate(chromophore_list):
            if chromophore.species != material_type:
                continue
            stack_dict[chromophore.ID] = stack_list[index]
        stack_dicts.append(stack_dict)
    return stack_dicts


def create_neighbour_list(chromophore_list, morphology_shape, cut_off, periodic, material_type):
    neighbour_dict = {}
    for chromo1 in chromophore_list:
        for [chromo2ID, rel_image] in chromo1.neighbours:
            if periodic is False:
                if not np.array_equal(rel_image, [0, 0, 0]):
                    continue
            if chromo1.species != material_type:
                continue
            chromo1posn = chromo1.posn
            chromo2posn = np.array(chromophore_list[chromo2ID].posn)\
                + (np.array(rel_image) * np.array(morphology_shape))
            separation = np.linalg.norm(chromo2posn - chromo1posn)
            if separation < cut_off:
                if chromo1.ID in neighbour_dict.keys():
                    neighbour_dict[chromo1.ID].append(chromo2ID)
                else:
                    neighbour_dict[chromo1.ID] = [chromo2ID]
    return neighbour_dict


def update_stack(atom_ID, cluster_list, neighbour_dict):
    try:
        for neighbour in neighbour_dict[atom_ID]:
            if cluster_list[neighbour] > cluster_list[atom_ID]:
                cluster_list[neighbour] = cluster_list[atom_ID]
                cluster_list = update_stack(neighbour, cluster_list, neighbour_dict)
            elif cluster_list[neighbour] < cluster_list[atom_ID]:
                cluster_list[atom_ID] = cluster_list[neighbour]
                cluster_list = update_stack(neighbour, cluster_list, neighbour_dict)
    except KeyError:
        pass
    return cluster_list


def plot_Stacks3D(output_dir, chromophore_list, stack_dicts, sim_dims):
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    colours = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'w']
    stack_dict = {}
    for dictionary in stack_dicts:
        if dictionary is not None:
            stack_dict.update(dictionary)
    stack_list = {}
    for chromophore in chromophore_list:
        stack_ID = stack_dict[chromophore.ID]
        if stack_ID not in stack_list.keys():
            stack_list[stack_ID] = [chromophore]
        else:
            stack_list[stack_ID].append(chromophore)
    for stack_ID, chromos in enumerate(stack_list.values()):
        for chromo in chromos:
            if chromo.species == 'Donor':
                ax.scatter(chromo.posn[0], chromo.posn[1], chromo.posn[2], facecolors='w',
                           edgecolors=colours[stack_ID % 8], s=40)
            elif chromo.species == 'Acceptor':
                ax.scatter(chromo.posn[0], chromo.posn[1], chromo.posn[2], c=colours[stack_ID % 8],
                           edgecolors=None, s=40)
    # Draw boxlines
    # Varying X
    ax.plot([sim_dims[0][0], sim_dims[0][1]], [sim_dims[1][0], sim_dims[1][0]], [sim_dims[2][0], sim_dims[2][0]],
            c='k', linewidth=1.0)
    ax.plot([sim_dims[0][0], sim_dims[0][1]], [sim_dims[1][1], sim_dims[1][1]], [sim_dims[2][0], sim_dims[2][0]],
            c='k', linewidth=1.0)
    ax.plot([sim_dims[0][0], sim_dims[0][1]], [sim_dims[1][0], sim_dims[1][0]], [sim_dims[2][1], sim_dims[2][1]],
            c='k', linewidth=1.0)
    ax.plot([sim_dims[0][0], sim_dims[0][1]], [sim_dims[1][1], sim_dims[1][1]], [sim_dims[2][1], sim_dims[2][1]],
            c='k', linewidth=1.0)
    # Varying Y
    ax.plot([sim_dims[0][0], sim_dims[0][0]], [sim_dims[1][0], sim_dims[1][1]], [sim_dims[2][0], sim_dims[2][0]],
            c='k', linewidth=1.0)
    ax.plot([sim_dims[0][1], sim_dims[0][1]], [sim_dims[1][0], sim_dims[1][1]], [sim_dims[2][0], sim_dims[2][0]],
            c='k', linewidth=1.0)
    ax.plot([sim_dims[0][0], sim_dims[0][0]], [sim_dims[1][0], sim_dims[1][1]], [sim_dims[2][1], sim_dims[2][1]],
            c='k', linewidth=1.0)
    ax.plot([sim_dims[0][1], sim_dims[0][1]], [sim_dims[1][0], sim_dims[1][1]], [sim_dims[2][1], sim_dims[2][1]],
            c='k', linewidth=1.0)
    # Varying Z
    ax.plot([sim_dims[0][0], sim_dims[0][0]], [sim_dims[1][0], sim_dims[1][0]], [sim_dims[2][0], sim_dims[2][1]],
            c='k', linewidth=1.0)
    ax.plot([sim_dims[0][0], sim_dims[0][0]], [sim_dims[1][1], sim_dims[1][1]], [sim_dims[2][0], sim_dims[2][1]],
            c='k', linewidth=1.0)
    ax.plot([sim_dims[0][1], sim_dims[0][1]], [sim_dims[1][0], sim_dims[1][0]], [sim_dims[2][0], sim_dims[2][1]],
            c='k', linewidth=1.0)
    ax.plot([sim_dims[0][1], sim_dims[0][1]], [sim_dims[1][1], sim_dims[1][1]], [sim_dims[2][0], sim_dims[2][1]],
            c='k', linewidth=1.0)
    ax.set_xlim([sim_dims[0][0], sim_dims[0][1]])
    ax.set_ylim([sim_dims[1][0], sim_dims[1][1]])
    ax.set_zlim([sim_dims[2][0], sim_dims[2][1]])
    plt.savefig(output_dir + "/02_stacks.pdf", bbox_inches='tight')
    plt.close()
    print("3D Stack figure saved as", output_dir + "/02_stacks.pdf")


def determine_molecule_IDs(CG_to_AAID_master, AA_morphology_dict, parameter_dict, chromophore_list):
    print("Determining molecule IDs...")
    CGID_to_mol_ID = {}
    if CG_to_AAID_master is not None:
        # Normal operation with a CGMorphology defined (fine-graining was performed)
        for mol_ID, mol_dict in enumerate(CG_to_AAID_master):
            for CGID in list(mol_dict.keys()):
                CGID_to_mol_ID[CGID] = mol_ID
    elif (len(parameter_dict['CG_site_species']) == 1)\
            and (('AA_rigid_body_species' not in parameter_dict)
                 or (len(parameter_dict['AA_rigid_body_species']) == 0)):
        print("Small-molecule system detected, assuming each chromophore is its own molecule...")
        # When CGMorphology doesn't exist, and no rigid body species have been specified, then
        # every chromophore is its own molecule)
        for index, chromo in enumerate(chromophore_list):
            for CGID in chromo.CGIDs:
                CGID_to_mol_ID[CGID] = chromo.ID
    else:
        # No CGMorphology, but not small molecules either, so determine molecules based on bonds
        print("Polymeric system detected, determining molecules based on AA bonds (slow calculation)...")
        molecule_AAIDs, molecule_lengths = split_molecules(AA_morphology_dict)
        for index, molecule_AAID_list in enumerate(molecule_AAIDs):
            for AAID in molecule_AAID_list:
                CGID_to_mol_ID[AAID] = index
    return CGID_to_mol_ID


def plot_energy_levels(output_dir, chromophore_list, data_dict):
    HOMO_levels = []
    LUMO_levels = []
    donor_delta_Eij = []
    acceptor_delta_Eij = []
    for chromo in chromophore_list:
        if chromo.species == 'Donor':
            HOMO_levels.append(chromo.HOMO)
            for neighbour_index, delta_Eij in enumerate(chromo.neighbours_delta_E):
                if (delta_Eij is not None) and (chromo.neighbours_TI[neighbour_index] is not None):
                    donor_delta_Eij.append(delta_Eij)
        else:
            LUMO_levels.append(chromo.LUMO)
            for neighbour_index, delta_Eij in enumerate(chromo.neighbours_delta_E):
                if (delta_Eij is not None) and (chromo.neighbours_TI[neighbour_index] is not None):
                    acceptor_delta_Eij.append(delta_Eij)
    if len(donor_delta_Eij) > 0:
        donor_bin_edges, donor_fit_args, donor_mean, donor_STD = gauss_fit(donor_delta_Eij)
        data_dict['donor_delta_Eij_mean'] = donor_mean
        data_dict['donor_delta_Eij_std'] = donor_STD
        data_dict['donor_delta_Eij_err'] = donor_STD / np.sqrt(len(donor_delta_Eij))
        HOMO_av = np.average(HOMO_levels)
        HOMO_std = np.std(HOMO_levels)
        HOMO_err = HOMO_std / np.sqrt(len(HOMO_levels))
        data_dict['donor_frontier_MO_mean'] = HOMO_av
        data_dict['donor_frontier_MO_std'] = HOMO_std
        data_dict['donor_frontier_MO_err'] = HOMO_err
        print("Donor HOMO Level =", HOMO_av, "+/-", HOMO_err)
        print("Donor Delta Eij stats: mean =", donor_mean, "+/-", donor_STD / np.sqrt(len(donor_delta_Eij)))
        plot_delta_Eij(donor_delta_Eij, donor_bin_edges, donor_fit_args, 'Donor',
                       output_dir + '/05_donor_delta_Eij.pdf')
    if len(acceptor_delta_Eij) > 0:
        acceptor_bin_edges, acceptor_fit_args, acceptor_mean, acceptor_STD = gauss_fit(acceptor_delta_Eij)
        data_dict['acceptor_delta_Eij_mean'] = acceptor_mean
        data_dict['acceptor_delta_Eij_std'] = acceptor_STD
        data_dict['acceptor_delta_Eij_err'] = acceptor_STD / np.sqrt(len(acceptor_delta_Eij))
        LUMO_av = np.average(LUMO_levels)
        LUMO_std = np.std(LUMO_levels)
        LUMO_err = LUMO_std / np.sqrt(len(LUMO_levels))
        data_dict['acceptor_frontier_MO_mean'] = LUMO_av
        data_dict['acceptor_frontier_MO_std'] = LUMO_std
        data_dict['acceptor_frontier_MO_err'] = LUMO_err
        print("Acceptor LUMO Level =", LUMO_av, "+/-", LUMO_err)
        print("Acceptor Delta Eij stats: mean =", acceptor_mean, "+/-",
              acceptor_STD / np.sqrt(len(acceptor_delta_Eij)))
        plot_delta_Eij(acceptor_delta_Eij, acceptor_bin_edges, acceptor_fit_args, 'Acceptor',
                       output_dir + '/06_acceptor_delta_Eij.pdf')
    return data_dict


def generate_data_dict():
        materials = ['donor', 'acceptor']
        material_inspecific_properties = ['name', 'density']
        hop_types = ['intra', 'inter']
        hop_targets = ['mol', 'stack']
        hop_dependent_properties = ['hops', 'proportion']
        error_properties = ['frontier_MO', 'delta_Eij']
        dictionary_elements = [(prop, '---') for prop in material_inspecific_properties]
        dictionary_elements += [(material + '_' + hop_type + '_' + hop_target + '_' + hop_property, '---')
                                for material in materials for hop_type in hop_types for hop_target in hop_targets
                                for hop_property in hop_dependent_properties]
        dictionary_elements += [(material + '_' + error_property + '_' + stat, '---') for material in materials
                                for error_property in error_properties for stat in ['mean', 'std', 'err']]
        data_dict = OrderedDict(dictionary_elements)
        return data_dict


def plot_delta_Eij(delta_Eij, gauss_bins, fit_args, data_type, file_name):
    plt.figure()
    n, bins, patches = plt.hist(delta_Eij, np.linspace(-0.5, 0.5, 20), color=['b'])
    if fit_args is not None:
        gauss_Y = gaussian(gauss_bins[:-1], *fit_args)
        scale_factor = max(n) / max(gauss_Y)
        plt.plot(gauss_bins[:-1], gauss_Y * scale_factor, 'ro:')
    else:
        print("No Gaussian found (probably zero-width delta function)")
    plt.ylabel('Frequency (Arb. U.)')
    plt.xlabel(data_type + r' $\Delta E_{ij}$ (eV)')
    plt.xlim([-0.5, 0.5])
    plt.savefig(file_name)
    plt.close()
    print("Figure saved as", file_name)


def plot_mixed_hopping_rates(output_dir, chromophore_list, parameter_dict, stack_dicts, CG_to_mol_ID, data_dict,
                             AA_morphology_dict):
    # Create all the empty lists we need
    hop_types = ['intra', 'inter']
    hop_targets = ['stack', 'mol']
    hop_properties = ['rates', 'TIs']
    chromo_species = ['donor', 'acceptor']
    property_lists = {}
    for property_name in [hop_type + '_' + hop_target + '_' + hop_property + '_' + species for hop_type in hop_types
                          for hop_target in hop_targets for hop_property in hop_properties
                          for species in chromo_species]:
        property_lists[property_name] = []
    try:
        if parameter_dict['reorganisation_energy_donor'] is not None:
            donor_lambdaij = parameter_dict['reorganisation_energy_donor']
        if parameter_dict['reorganisation_energy_acceptor'] is not None:
            acceptor_lambdaij = parameter_dict['reorganisation_energy_acceptor']
    except KeyError:  # Old MorphCT fix
        print("Only one reorganisation energy found, assuming donor and continuing")
        donor_lambdaij = parameter_dict['reorganisation_energy']
    T = 290
    for chromo in chromophore_list:
        mol1ID = CG_to_mol_ID[chromo.CGIDs[0]]
        for index, Tij in enumerate(chromo.neighbours_TI):
            if (Tij is None) or (Tij == 0):
                continue
            chromo2 = chromophore_list[chromo.neighbours[index][0]]
            mol2ID = CG_to_mol_ID[chromo2.CGIDs[0]]
            delta_E = chromo.neighbours_delta_E[index]
            if chromo.species == 'acceptor':
                lambdaij = acceptor_lambdaij
            else:
                lambdaij = donor_lambdaij
            # Now take into account the various behaviours we can have from the parameter file
            prefactor = 1.0
            # Apply the koopmans prefactor
            try:
                use_koop = parameter_dict['use_koopmans_approximation']
                if use_koop:
                    prefactor *= parameter_dict['koopmans_hopping_prefactor']
            except KeyError:
                pass
            # Apply the simple energetic penalty model
            try:
                boltz_pen = parameter_dict['use_simple_energetic_penalty']
            except KeyError:
                boltz_pen = False
            # Apply the distance penalty due to VRH
            try:
                VRH = parameter_dict['use_VRH']
                if VRH is True:
                    VRH_prefactor = 1.0 / parameter_dict['VRH_delocalisation']
            except KeyError:
                VRH = False
            if VRH is True:
                relative_image = chromo.neighbours[index][1]
                neighbour_chromo_posn = chromo2.posn + (np.array(relative_image)
                                                        * np.array([AA_morphology_dict[axis]
                                                                    for axis in ['lx', 'ly', 'lz']]))
                chromophore_separation = hf.calculate_separation(chromo.posn, neighbour_chromo_posn) * 1E-10
                rate = hf.calculate_carrier_hop_rate(lambdaij * elementary_charge, Tij * elementary_charge,
                                                     delta_E * elementary_charge, prefactor, T, use_VRH=VRH,
                                                     rij=chromophore_separation, VRH_prefactor=VRH_prefactor,
                                                     boltz_pen=boltz_pen)
            else:
                rate = hf.calculate_carrier_hop_rate(lambdaij * elementary_charge, Tij * elementary_charge,
                                                     delta_E * elementary_charge, prefactor, T, boltz_pen=boltz_pen)
            if chromo2.ID < chromo.ID:
                continue
            # Do intra- / inter- stacks
            if chromo.species == 'acceptor':
                if stack_dicts[1][chromo.ID] == stack_dicts[1][chromo.neighbours[index][0]]:
                    property_lists['intra_stack_rates_acceptor'].append(rate)
                    property_lists['intra_stack_TIs_acceptor'].append(Tij)
                else:
                    property_lists['inter_stack_rates_acceptor'].append(rate)
                    property_lists['inter_stack_TIs_acceptor'].append(Tij)
            else:
                if stack_dicts[0][chromo.ID] == stack_dicts[0][chromo.neighbours[index][0]]:
                    property_lists['intra_stack_rates_donor'].append(rate)
                    property_lists['intra_stack_TIs_donor'].append(Tij)
                else:
                    property_lists['inter_stack_rates_donor'].append(rate)
                    property_lists['inter_stack_TIs_donor'].append(Tij)
            # Now do intra- / inter- molecules
            if mol1ID == mol2ID:
                if chromo.species == 'acceptor':
                    property_lists['intra_mol_rates_acceptor'].append(rate)
                    property_lists['intra_mol_TIs_acceptor'].append(Tij)
                else:
                    property_lists['intra_mol_rates_donor'].append(rate)
                    property_lists['intra_mol_TIs_donor'].append(Tij)
            else:
                if chromo.species == 'acceptor':
                    property_lists['inter_mol_rates_acceptor'].append(rate)
                    property_lists['inter_mol_TIs_acceptor'].append(Tij)
                else:
                    property_lists['inter_mol_rates_donor'].append(rate)
                    property_lists['inter_mol_TIs_donor'].append(Tij)
    # Donor Stack Plots:
    if (len(property_lists['intra_stack_rates_donor']) > 0) or (len(property_lists['inter_stack_rates_donor']) > 0):
        print("Mean intra-stack Donor rate =", np.mean(property_lists['intra_stack_rates_donor']), "+/-",
              np.std(property_lists['intra_stack_rates_donor'])
              / float(len(property_lists['intra_stack_rates_donor'])))
        print("Mean inter-stack Donor rate =", np.mean(property_lists['inter_stack_rates_donor']), "+/-",
              np.std(property_lists['inter_stack_rates_donor'])
              / float(len(property_lists['inter_stack_rates_donor'])))
        plot_stacked_hist_rates(property_lists['intra_stack_rates_donor'], property_lists['inter_stack_rates_donor'],
                                ['Intra-stack', 'Inter-stack'], 'Donor',
                                output_dir + '/16_donor_hopping_rate_stacks.pdf')
        plot_stacked_hist_TIs(property_lists['intra_stack_TIs_donor'], property_lists['inter_stack_TIs_donor'],
                              ['Intra-stack', 'Inter-stack'], 'Donor',
                              output_dir + '/12_donor_transfer_integral_stacks.pdf')
    # Acceptor Stack Plots:
    if (len(property_lists['intra_stack_rates_acceptor']) > 0)\
       or (len(property_lists['inter_stack_rates_acceptor']) > 0):
        print("Mean intra-stack Acceptor rate =", np.mean(property_lists['intra_stack_rates_acceptor']), "+/-",
              np.std(property_lists['intra_stack_rates_acceptor'])
              / float(len(property_lists['intra_stack_rates_acceptor'])))
        print("Mean inter-stack Acceptor rate =", np.mean(property_lists['inter_stack_rates_acceptor']), "+/-",
              np.std(property_lists['inter_stack_rates_acceptor'])
              / float(len(property_lists['inter_stack_rates_acceptor'])))
        plot_stacked_hist_rates(property_lists['intra_stack_rates_acceptor'],
                                property_lists['inter_stack_rates_acceptor'], ['Intra-stack', 'Inter-stack'],
                                'Acceptor', output_dir + '/18_acceptor_hopping_rate_stacks.pdf')
        plot_stacked_hist_TIs(property_lists['intra_stack_TIs_acceptor'], property_lists['inter_stack_TIs_acceptor'],
                              ['Intra-stack', 'Inter-stack'], 'Acceptor',
                              output_dir + '/14_acceptor_transfer_integral_stacks.pdf')
    # Donor Mol Plots:
    if (len(property_lists['intra_mol_rates_donor']) > 0) or (len(property_lists['inter_mol_rates_donor']) > 0):
        print("Mean intra-molecular donor rate =", np.mean(property_lists['intra_mol_rates_donor']), "+/-",
              np.std(property_lists['intra_mol_rates_donor']) / float(len(property_lists['intra_mol_rates_donor'])))
        print("Mean inter-molecular donor rate =", np.mean(property_lists['inter_mol_rates_donor']), "+/-",
              np.std(property_lists['inter_mol_rates_donor']) / float(len(property_lists['inter_mol_rates_donor'])))
        plot_stacked_hist_rates(property_lists['intra_mol_rates_donor'], property_lists['inter_mol_rates_donor'],
                                ['Intra-mol', 'Inter-mol'], 'Donor', output_dir + '/15_donor_hopping_rate_mols.pdf')
        plot_stacked_hist_TIs(property_lists['intra_mol_TIs_donor'], property_lists['inter_mol_TIs_donor'],
                              ['Intra-mol', 'Inter-mol'], 'Donor', output_dir + '/11_donor_transfer_integral_mols.pdf')
    # Acceptor Mol Plots:
    if (len(property_lists['intra_mol_rates_acceptor']) > 0) or (len(property_lists['inter_mol_rates_acceptor']) > 0):
        print("Mean intra-molecular acceptor rate =", np.mean(property_lists['intra_mol_rates_acceptor']), "+/-",
              np.std(property_lists['intra_mol_rates_acceptor'])
              / float(len(property_lists['intra_mol_rates_acceptor'])))
        print("Mean inter-molecular acceptor rate =", np.mean(property_lists['inter_mol_rates_acceptor']), "+/-",
              np.std(property_lists['inter_mol_rates_acceptor'])
              / float(len(property_lists['inter_mol_rates_acceptor'])))
        plot_stacked_hist_rates(property_lists['intra_mol_rates_acceptor'], property_lists['inter_mol_rates_acceptor'],
                                ['Intra-mol', 'Inter-mol'], 'Acceptor',
                                output_dir + '/17_acceptor_hopping_rate_mols.pdf')
        plot_stacked_hist_TIs(property_lists['intra_mol_TIs_acceptor'], property_lists['inter_mol_TIs_acceptor'],
                              ['Intra-mol', 'Inter-mol'], 'Acceptor',
                              output_dir + '/13_acceptor_transfer_integral_mols.pdf')
    # Update the dataDict
    for material in chromo_species:
        for hop_type in hop_types:
            for hop_target in hop_targets:
                number_of_hops = len(property_lists[hop_type + '_' + hop_target + '_rates_' + material])
                if number_of_hops == 0:
                    continue
                other_hop_type = hop_types[int((hop_types.index(hop_type) * -1) + 1)]
                proportion = number_of_hops / (number_of_hops
                                               + len(property_lists[other_hop_type + '_' + hop_target + '_rates_'
                                                                    + material]))
                mean_rate = np.mean(property_lists[hop_type + '_' + hop_target + '_rates_' + material])
                data_dict[material.lower() + '_' + hop_type + '_' + hop_target.lower() + "_hops"] = number_of_hops
                data_dict[material.lower() + '_' + hop_type + '_' + hop_target.lower() + "_proportion"] = proportion
                data_dict[material.lower() + '_' + hop_type + '_' + hop_target.lower() + "_rate_mean"] = mean_rate
    return data_dict


def plot_stacked_hist_rates(data1, data2, labels, data_type, file_name):
    plt.figure()
    (n, bins, patches) = plt.hist([data1, data2], bins=np.logspace(1, 18, 40), stacked=True, color=['r', 'b'],
                                  label=labels)
    plt.ylabel('Frequency (Arb. U.)')
    plt.xlabel(data_type + r' k$_{ij}$ (s' + r'$^{-1}$' + ')')
    plt.xlim([1, 1E18])
    plt.xticks([1E0, 1E3, 1E6, 1E9, 1E12, 1E15, 1E18])
    plt.ylim([0, np.max(n) * 1.02])
    plt.legend(loc=0, prop={'size': 18})
    plt.gca().set_xscale('log')
    plt.savefig(file_name)
    plt.close()
    print("Figure saved as", file_name)


def plot_stacked_hist_TIs(data1, data2, labels, data_type, file_name):
    plt.figure()
    (n, bins, patches) = plt.hist([data1, data2], bins=np.linspace(0, 1.2, 20), stacked=True, color=['r', 'b'],
                                  label=labels)
    plt.ylabel('Frequency (Arb. U.)')
    plt.xlabel(data_type + r' T$_{ij}$ (eV)')
    plt.xlim([0, 1.2])
    plt.ylim([0, np.max(n) * 1.02])
    plt.legend(loc=0, prop={'size': 18})
    plt.savefig(file_name)
    plt.close()
    print("Figure saved as", file_name)


def write_CSV(data_dict, directory):
    CSV_file_name = directory + '/results.csv'
    with open(CSV_file_name, 'w+') as CSV_file:
        CSV_writer = csv.writer(CSV_file)
        for key, val in data_dict.items():
            CSV_writer.writerow([key, val])
    print("CSV file written to " + CSV_file_name)


def create_results_pickle(directory):
    cores_list = []
    for file_name in glob.glob(directory + '/KMC/*'):
        try:
            cores_list.append(re.findall("([_])(..)([\.])", file_name)[0][1])
        except IndexError:
            pass
    cores_list = sorted(list(set(cores_list)))
    results_pickles_list = []
    keep_list = []
    for core in cores_list:
        # Check if there is already a finished KMC_results pickle
        main = directory + '/KMC/KMC_results_%02d.pickle' % (int(core))
        if os.path.exists(main):
            results_pickles_list.append(main)
            keep_list.append(None)
            continue
        # If not, find the slot1 and slot2 pickle that is most recent
        slot1 = directory + '/KMC/KMC_slot1_results_%02d.pickle' % (int(core))
        slot2 = directory + '/KMC/KMC_slot2_results_%02d.pickle' % (int(core))
        if (os.path.exists(slot1) and not os.path.exists(slot2)):
            keep_list.append(slot1)
        elif (os.path.exists(slot2) and not os.path.exists(slot1)):
            keep_list.append(slot2)
        elif os.path.getsize(slot1) >= os.path.getsize(slot2):
            keep_list.append(slot1)
        else:
            keep_list.append(slot2)
    print("%d pickle files found to combine!" % (len(keep_list)))
    print("Combining", keep_list)
    for keeper in zip(cores_list, keep_list):
        # Skip this core if we already have a finished KMC_results for it
        if keeper[1] is None:
            continue
        new_name = directory + '/KMC/KMC_results_' + str(keeper[0]) + '.pickle'
        shutil.copyfile(str(keeper[1]), new_name)
        results_pickles_list.append(new_name)
    combine_results_pickles(directory, results_pickles_list)


def combine_results_pickles(directory, pickle_files):
    combined_data = {}
    pickle_files = sorted(pickle_files)
    for file_name in pickle_files:
        # The pickle was repeatedly dumped to, in order to save time.
        # Each dump stream is self-contained, so iteratively unpickle to add the new data.
        with open(file_name, 'rb') as pickle_file:
            pickled_data = pickle.load(pickle_file)
            for key, val in pickled_data.items():
                if val is None:
                    continue
                if key not in combined_data:
                    combined_data[key] = val
                else:
                    combined_data[key] += val
    # Write out the combined data
    print("Writing out the combined pickle file...")
    with open(directory + '/KMC/KMC_results.pickle', 'wb+') as pickle_file:
        pickle.dump(combined_data, pickle_file)
    print("Complete data written to", directory + "/KMC_results.pickle.")


def calculate_mobility(directory, current_carrier_type, carrier_data, sim_dims, plot3D_graphs, chromophore_list):
    print("Considering the transport of", current_carrier_type + "...")
    print("Obtaining mean squared displacements...")
    carrier_history, times, MSDs, time_standard_errors, MSD_standard_errors = get_carrier_data(carrier_data)
    print("MSDs obtained")
    # Create the first figure that will be replotted each time
    plt.figure()
    anisotropy = plot_anisotropy(carrier_data, directory, sim_dims, current_carrier_type, plot3D_graphs)
    if (carrier_history is not None) and plot3D_graphs:
        print("Determining carrier hopping connections...")
        plot_connections(chromophore_list, sim_dims, carrier_history, directory, current_carrier_type)
    times, MSDs = hf.parallel_sort(times, MSDs)
    print("Calculating MSD...")
    mobility, mob_error, r_squared = plot_MSD(times, MSDs, time_standard_errors, MSD_standard_errors, directory,
                                              current_carrier_type)
    print("----------====================----------")
    print(current_carrier_type, "mobility for", directory,
          "= %.2E +- %.2E cm^{2} V^{-1} s^{-1}" % (mobility, mob_error))
    print("----------====================----------")
    return mobility, mob_error, r_squared, anisotropy


def KMCAnalyse():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--three_D", action="store_true", required=False,
                        help=('If present, use matplotlib to plot the 3D graphs (3D network, anisotropy and stack'
                              ' positions. This takes a while (usually a couple of minutes) to plot.'
                              ' Defaults to False.'))
    parser.add_argument("-p", "--periodic_stacks", action="store_true", required=False,
                        help=('If present, allow periodic connections to add chromophores to stacks, as well as'
                              ' non-periodic connections (this usually reduces the number of stacks in the system).'
                              ' Defaults to False.'))
    parser.add_argument("-d", "--cut_off_donor", type=float, default=None, required=False,
                        help=('Specify a manual cut-off for the determination of stacks within the donor material.'
                              'Connections with separation > cut-off will be classed as inter-stack. If omitted, stack'
                              'cut-off will be determined automatically as the first minimum of the RDF.'))
    parser.add_argument("-a", "--cut_off_acceptor", type=float, default=None, required=False,
                        help=('Specify a manual cut-off for the determination of stacks within the acceptor material.'
                              'Connections with separation > cut-off will be classed as inter-stack. if omitted, stack'
                              'cut-off will be determined automatically as the first minimum of the RDF.'))
    parser.add_argument("-s", "--sequence", type=lambda s: [float(item) for item in s.split(',')], default=None,
                        required=False, help=('Create a figure in the current directory that describes the evolution'
                                              ' of the anisotropy/mobility using the specified comma-delimited string'
                                              ' as the sequence of x values. For instance -s "1.5,1.75,2.0,2.25,2.5"'
                                              ' will assign each of the 5 following directories these x-values when'
                                              ' plotting the mobility evolution.'))
    parser.add_argument("-x", "--xlabel", default="temperature (Arb. U.)", required=False,
                        help=('Specify an x-label for the combined plot (only used if -s is specified). Default ='
                              '"temperature (Arb. U.)"'))
    args, directory_list = parser.parse_known_args()

    sys.setrecursionlimit(10000)
    hole_mobility_data = []
    hole_anisotropy_data = []
    electron_mobility_data = []
    electron_anisotropy_data = []
    for directory in directory_list:
        # Create the figures directory if it doesn't already exist
        os.makedirs(directory + '/figures', exist_ok=True)
        # Now create the data dictionary
        data_dict = generate_data_dict()
        print("\n")
        print("Getting carrier data...")
        carrier_data = load_KMC_results_pickle(directory)
        print("Carrier Data obtained")
        # Now need to split up the carrierData into both electrons and holes
        carrier_data_holes, carrier_data_electrons = split_carriers_by_type(carrier_data)
        print("Loading chromophoreList...")
        data = hf.load_pickle(directory + '/code/' + directory + '.pickle')
        AA_morphology_dict = data[0]
        CG_morphology_dict = data[1]
        CG_to_AAID_master = data[2]
        parameter_dict = data[3]
        chromophore_list = data[4]
        print("ChromophoreList obtained")
        morphology_shape = np.array([AA_morphology_dict[axis] for axis in ['lx', 'ly', 'lz']])
        sim_dims = [[-AA_morphology_dict[axis] / 2.0, AA_morphology_dict[axis] / 2.0] for axis in ['lx', 'ly', 'lz']]
        # Calculate the mobilities
        complete_carrier_types = []
        complete_carrier_data = []
        if (carrier_data_holes is not None) and (len(carrier_data_holes['ID']) > 0):
            complete_carrier_types.append('hole')
            complete_carrier_data.append(carrier_data_holes)
        if (carrier_data_electrons is not None) and (len(carrier_data_electrons['ID']) > 0):
            complete_carrier_types.append('electron')
            complete_carrier_data.append(carrier_data_electrons)
        for carrier_type_index, carrier_data in enumerate(complete_carrier_data):
            current_carrier_type = complete_carrier_types[carrier_type_index]
            mobility, mob_error, r_squared, anisotropy = calculate_mobility(directory, current_carrier_type,
                                                                            carrier_data, sim_dims, args.three_D,
                                                                            chromophore_list)
            if current_carrier_type == 'hole':
                hole_anisotropy_data.append(anisotropy)
                hole_mobility_data.append([mobility, mob_error])
            elif current_carrier_type == 'electron':
                electron_anisotropy_data.append(anisotropy)
                electron_mobility_data.append([mobility, mob_error])
            data_dict['name'] = directory
            data_dict[current_carrier_type.lower() + '_anisotropy'] = anisotropy
            data_dict[current_carrier_type.lower() + '_mobility'] = mobility
            data_dict[current_carrier_type.lower() + '_mobility_r_squared'] = r_squared
        # Now plot the distributions!
        temp_dir = directory + '/figures'
        CG_to_mol_ID = determine_molecule_IDs(CG_to_AAID_master, AA_morphology_dict, parameter_dict, chromophore_list)
        data_dict = plot_energy_levels(temp_dir, chromophore_list, data_dict)
        print("Finding the relevant stack cut off as the midpoint between the first maxmimum and the first minimum",
              "of the neighbour distance distribution...")
        print("Considering periodic neighbours is", args.periodic_stacks)
        cut_offs = get_neighbour_cut_off(chromophore_list, morphology_shape, temp_dir, periodic=args.periodic_stacks,
                                         specified_cut_off_donor=args.cut_off_donor,
                                         specified_cut_off_acceptor=args.cut_off_acceptor)
        calculated_cut_off_donor = cut_offs[0]
        calculated_cut_off_acceptor = cut_offs[1]
        if args.cut_off_donor is None:
            print("No donor cut off has been specified. Using calculated value...")
            cut_off_donor = calculated_cut_off_donor
        else:
            cut_off_donor = args.cut_off_donor
        if args.cut_off_acceptor is None:
            print("No acceptor cut off has been specified. Using calculated value...")
            cut_off_acceptor = calculated_cut_off_acceptor
        else:
            cut_off_acceptor = args.cut_off_acceptor
        print("Cut off in Angstroems (Donor) =", cut_off_donor)
        print("Cut off in Angstroems (Acceptor) =", cut_off_acceptor)
        stack_dicts = get_stacks(chromophore_list, morphology_shape, cut_off_donor, cut_off_acceptor,
                                 periodic=args.periodic_stacks)
        if args.three_D:
            plot_Stacks3D(temp_dir, chromophore_list, stack_dicts, sim_dims)
        data_dict = plot_mixed_hopping_rates(temp_dir, chromophore_list, parameter_dict, stack_dicts, CG_to_mol_ID,
                                             data_dict, AA_morphology_dict)
        print("\n")
        print("Writing CSV Output File...")
        write_CSV(data_dict, directory)
    print("Plotting Mobility and Anisotropy progressions...")
    if args.sequence is not None:
        if len(hole_anisotropy_data) > 0:
            plot_temperature_progression(args.sequence, hole_mobility_data, hole_anisotropy_data, 'hole', args.xlabel)
        if len(electron_anisotropy_data) > 0:
            plot_temperature_progression(args.sequence, electron_mobility_data, electron_anisotropy_data, 'electron',
                                         args.xlabel)
    else:
        print("Skipping plotting mobility evolution.")
