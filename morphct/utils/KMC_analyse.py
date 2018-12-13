import copy
import os
import sys
import pickle
import csv
import numpy as np
from scipy.optimize import curve_fit
import scipy.stats
from morphct.code import helper_functions as hf
from collections import OrderedDict
import shutil
import glob
import re
import argparse
import copy
from scipy.signal import argrelextrema
from scipy.ndimage import gaussian_filter
from scipy.sparse import lil_matrix


plt = None
p3 = None
elementary_charge = 1.60217657E-19  # C
kB = 1.3806488E-23  # m^{2} kg s^{-2} K^{-1}
hbar = 1.05457173E-34  # m^{2} kg s^{-1}
temperature = 290  # K


def load_KMC_results_pickle(directory):
    KMC_pickle = os.path.join(directory, "KMC", "KMC_results.pickle")
    try:
        with open(KMC_pickle, "rb") as pickle_file:
            carrier_data = pickle.load(pickle_file)
    except FileNotFoundError:
        print("No final KMC_results.pickle found. Creating it from incomplete parts...")
        create_results_pickle(directory)
        with open(KMC_pickle, "rb") as pickle_file:
            carrier_data = pickle.load(pickle_file)
    except UnicodeDecodeError:
        with open(KMC_pickle, "rb") as pickle_file:
            carrier_data = pickle.load(pickle_file, encoding="latin1")
    return carrier_data


def split_carriers_by_type(carrier_data):
    # If only one carrier type has been given, call the carriers holes and skip the electron calculations
    list_variables = [
        "current_time",
        "ID",
        "no_hops",
        "displacement",
        "lifetime",
        "final_position",
        "image",
        "initial_position",
    ]
    try:
        carrier_data_holes = {
            "carrier_history_matrix": carrier_data["hole_history_matrix"],
            "seed": carrier_data["seed"],
        }
        carrier_data_electrons = {
            "carrier_history_matrix": carrier_data["electron_history_matrix"],
            "seed": carrier_data["seed"],
        }
        for list_var in list_variables:
            carrier_data_holes[list_var] = []
            carrier_data_electrons[list_var] = []
            for carrier_index, charge_type in enumerate(carrier_data["carrier_type"]):
                if charge_type.lower() == "hole":
                    carrier_data_holes[list_var].append(
                        carrier_data[list_var][carrier_index]
                    )
                elif charge_type.lower() == "electron":
                    carrier_data_electrons[list_var].append(
                        carrier_data[list_var][carrier_index]
                    )
    except KeyError:
        print("This is an old-style pickle!")
        print(
            "Multiple charge carriers not found, assuming donor material and holes only."
        )
        try:
            carrier_data_holes = {
                "carrier_history_matrix": carrier_data["carrier_history_matrix"],
                "seed": carrier_data["seed"],
            }
        except KeyError:
            carrier_data_holes = {
                "carrier_history_matrix": carrier_data["carrier_history_matrix"],
                "seed": 0,
            }
        carrier_data_electrons = None
        for list_var in list_variables:
            carrier_data_holes[list_var] = []
            for carrier_index, carrier_ID in enumerate(carrier_data["ID"]):
                carrier_data_holes[list_var].append(
                    carrier_data[list_var][carrier_index]
                )
    return carrier_data_holes, carrier_data_electrons


def get_carrier_data(carrier_data):
    try:
        carrier_history = carrier_data["carrier_history_matrix"]
    except:
        carrier_history = None
    total_data_points = 0
    total_data_points_averaged_over = 0
    squared_disps = {}
    actual_times = {}
    for carrier_index, displacement in enumerate(carrier_data["displacement"]):
        if (
            (
                carrier_data["current_time"][carrier_index]
                > carrier_data["lifetime"][carrier_index] * 2
            )
            or (
                carrier_data["current_time"][carrier_index]
                < carrier_data["lifetime"][carrier_index] / 2.0
            )
            or (carrier_data["no_hops"][carrier_index] == 1)
        ):
            total_data_points += 1
            continue
        carrier_key = str(carrier_data["lifetime"][carrier_index])
        if carrier_key not in squared_disps:
            squared_disps[carrier_key] = [
                (carrier_data["displacement"][carrier_index] * 1E-10) ** 2
            ]  # A -> m
            actual_times[carrier_key] = [carrier_data["current_time"][carrier_index]]
        else:
            squared_disps[carrier_key].append(
                (carrier_data["displacement"][carrier_index] * 1E-10) ** 2
            )  # A -> m
            actual_times[carrier_key].append(
                carrier_data["current_time"][carrier_index]
            )
        # Also keep track of whether each carrier is a hole or an electron
        total_data_points_averaged_over += 1
        total_data_points += 1
    if total_data_points > total_data_points_averaged_over:
        print(
            "Notice: The data from",
            total_data_points - total_data_points_averaged_over,
            "carriers was discarded"
            "due to the carrier lifetime being more than double (or less than half of) the specified carrier lifetime.",
        )
    times = []
    MSDs = []
    time_standard_errors = []
    MSD_standard_errors = []
    for time, disps in squared_disps.items():
        times.append(float(time))
        time_standard_errors.append(
            np.std(actual_times[time]) / len(actual_times[time])
        )
        MSDs.append(np.average(disps))
        MSD_standard_errors.append(np.std(disps) / len(disps))
    return carrier_history, times, MSDs, time_standard_errors, MSD_standard_errors


def plot_displacement_dist(carrier_data, directory, carrier_type):
    carrier_types = ["hole", "electron"]
    plt.figure()
    plt.hist(np.array(carrier_data["displacement"]) * 0.1, bins=60, color="b")
    plt.xlabel("".join([carrier_type, "Displacement (nm)"]))
    plt.ylabel("Frequency (Arb. U.)")
    # 30 for hole displacement dist, 31 for electron displacement dist
    file_name = "".join(
        [
            "{:02}_".format(30 + carrier_types.index(carrier_type)),
            carrier_type,
            "_displacement_dist.png",
        ]
    )
    plt.savefig(os.path.join(directory, "figures", file_name), dpi=300)
    print("Figure saved as", os.path.join(directory, "figures", file_name))


def plot_cluster_size_dist(cluster_freqs, directory):
    carrier_types = ["hole", "electron"]
    for carrier_type_index, carrier_type in enumerate(carrier_types):
        try:
            sizes = list(cluster_freqs[carrier_type_index].values())
            sizes = [np.log10(size) for size in sizes if size > 5]
            if len(sizes) == 0:
                raise IndexError
        except IndexError:
            return
        plt.figure()
        try:
            plt.hist(sizes, bins=np.logspace(0, np.ceil(np.max(sizes)), 20), color="b")
            plt.xscale("log")
            # tick_labels = np.arange(0, np.ceil(np.max(sizes)) + 1, 1)
            plt.xlim([0, 10 ** (np.ceil(np.max(sizes)))])
            # plt.xticks(tick_labels, [r'10$^{{{}}}$'.format(int(x)) for x in tick_labels])
        except ValueError:
            print(
                "EXCEPTION: No clusters found. Are you sure the cluster criteria are correct?"
            )
        plt.xlabel("Cluster Size (Arb. U.)")
        plt.ylabel("Frequency (Arb. U.)")
        # 32 for hole cluster size dist, 33 for electron cluster size dist
        file_name = "".join(
            [
                "{:02}_".format(32 + carrier_types.index(carrier_type)),
                carrier_type,
                "_cluster_dist.png",
            ]
        )
        plt.savefig(os.path.join(directory, "figures", file_name), dpi=300)
        print("Figure saved as", os.path.join(directory, "figures", file_name))


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
    # Iterate through the chromophore_list
    for i, chromo in enumerate(chromophore_list):
        # Iterate through the neighbors of the chromophore
        for neighbor in zip(chromo.neighbours):
            index = neighbor[0][0]  # index of the neighbor
            image = neighbor[0][1]  # check to see if they are in the same rel image
            # Only consider one direction.
            if i < index:
                # Get the vector between the two chromophores.
                if not np.count_nonzero(image):
                    vector = chromophore_list[index].posn - chromo.posn
                # Account for periodic boundary conditions if not in same relative image.
                else:
                    vector = chromophore_list[index].posn - chromo.posn
                    vector += image * np.array(
                        [2 * sim_dims[0][1], 2 * sim_dims[1][1], 2 * sim_dims[2][1]]
                    )

                # Get the net number of times the path was travelled.
                forward = carrier_history[index, i]
                reverse = carrier_history[i, index]
                times_travelled = abs(forward - reverse)

                # Append the array if the net times travelled is greater than 0
                if times_travelled > 0:
                    datum = np.hstack(
                        (chromo.posn, vector, np.array([np.log10(times_travelled)]))
                    )
                    connections_array = np.vstack((connections_array, datum))
    return connections_array[1:]  # Return the array excluding the zeros first line.


def plot_connections(
    chromophore_list, sim_dims, carrier_history, directory, carrier_type
):
    # A complicated function that shows connections between carriers in 3D that carriers prefer to hop between.
    # Connections that are frequently used are highlighted in black, whereas rarely used connections are more white.
    # Import matplotlib color modules to set up color bar.
    import matplotlib.colors
    import matplotlib.cm as cmx

    # Create a figure class
    fig = plt.figure(figsize=(7, 6))
    # Make a 3D subplot
    ax = fig.add_subplot(111, projection="3d")

    # Create the array for all the chromophore connections
    connections_array = create_array_for_plot_connections(
        chromophore_list, carrier_history, sim_dims
    )

    # Determine the smalles, non-zero number of times two chromophores are connected.
    vmin = np.min(np.array(connections_array)[:, 6])
    # Determine the max number of times two chormophores are connected.
    vmax = np.max(np.array(connections_array)[:, 6])

    # Set up the color bar.
    colour_map = plt.get_cmap("inferno")
    c_norm = matplotlib.colors.Normalize(vmin=np.floor(vmin), vmax=np.ceil(vmax))
    scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=colour_map)
    hop_colors = scalar_map.to_rgba(connections_array[:, 6])

    # Set up the intensity for the hops so more travelled paths are more intense
    alphas = connections_array[:, 6] / vmax
    hop_colors[:, 3] = alphas

    # Plot the vectors between two chromophores
    ax.quiver(
        connections_array[:, 0],
        connections_array[:, 1],
        connections_array[:, 2],
        connections_array[:, 3],
        connections_array[:, 4],
        connections_array[:, 5],
        color=hop_colors,
        arrow_length_ratio=0,
        linewidth=0.7,
    )

    # Draw boxlines
    # Varying X
    ax.plot(
        [sim_dims[0][0], sim_dims[0][1]],
        [sim_dims[1][0], sim_dims[1][0]],
        [sim_dims[2][0], sim_dims[2][0]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims[0][0], sim_dims[0][1]],
        [sim_dims[1][1], sim_dims[1][1]],
        [sim_dims[2][0], sim_dims[2][0]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims[0][0], sim_dims[0][1]],
        [sim_dims[1][0], sim_dims[1][0]],
        [sim_dims[2][1], sim_dims[2][1]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims[0][0], sim_dims[0][1]],
        [sim_dims[1][1], sim_dims[1][1]],
        [sim_dims[2][1], sim_dims[2][1]],
        c="k",
        linewidth=1.0,
    )
    # Varying Y
    ax.plot(
        [sim_dims[0][0], sim_dims[0][0]],
        [sim_dims[1][0], sim_dims[1][1]],
        [sim_dims[2][0], sim_dims[2][0]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims[0][1], sim_dims[0][1]],
        [sim_dims[1][0], sim_dims[1][1]],
        [sim_dims[2][0], sim_dims[2][0]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims[0][0], sim_dims[0][0]],
        [sim_dims[1][0], sim_dims[1][1]],
        [sim_dims[2][1], sim_dims[2][1]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims[0][1], sim_dims[0][1]],
        [sim_dims[1][0], sim_dims[1][1]],
        [sim_dims[2][1], sim_dims[2][1]],
        c="k",
        linewidth=1.0,
    )
    # Varying Z
    ax.plot(
        [sim_dims[0][0], sim_dims[0][0]],
        [sim_dims[1][0], sim_dims[1][0]],
        [sim_dims[2][0], sim_dims[2][1]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims[0][0], sim_dims[0][0]],
        [sim_dims[1][1], sim_dims[1][1]],
        [sim_dims[2][0], sim_dims[2][1]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims[0][1], sim_dims[0][1]],
        [sim_dims[1][0], sim_dims[1][0]],
        [sim_dims[2][0], sim_dims[2][1]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims[0][1], sim_dims[0][1]],
        [sim_dims[1][1], sim_dims[1][1]],
        [sim_dims[2][0], sim_dims[2][1]],
        c="k",
        linewidth=1.0,
    )

    # Make the colour bar
    scalar_map.set_array(connections_array[:, 6])
    # tick_location = np.arange(0, int(np.ceil(vmax)) + 1, 1)
    tick_location = np.arange(0, vmax, 1)
    # Plot the colour bar
    cbar = plt.colorbar(
        scalar_map, ticks=tick_location, shrink=0.8, aspect=20
    )  # , fraction=0.046, pad=0.04)
    cbar.ax.set_yticklabels([r"10$^{{{}}}$".format(int(x)) for x in tick_location])

    # Name and save the figure.
    carrier_types = ["hole", "electron"]
    materials_types = ["donor", "acceptor"]
    carrier_index = carrier_types.index(carrier_type)
    figure_title = " ".join(
        [
            materials_types[carrier_index].capitalize(),
            "".join(["(", carrier_type.capitalize(), ")"]),
            "Network",
        ]
    )
    # plt.title(figure_title, y=1.1)
    # 01 for donor 3d network, 02 for acceptor 3d network
    file_name = "".join(["{:02}_3d_".format(1 + carrier_index), carrier_type, ".png"])
    plt.savefig(
        os.path.join(directory, "figures", file_name), bbox_inches="tight", dpi=300
    )
    print("Figure saved as", os.path.join(directory, "figures", file_name))

    plt.clf()


def calc_mobility(lin_fit_X, lin_fit_Y, av_time_error, av_MSD_error):
    # YVals have a std error avMSDError associated with them
    # XVals have a std error avTimeError assosciated with them
    numerator = lin_fit_Y[-1] - lin_fit_Y[0]
    denominator = lin_fit_X[-1] - lin_fit_X[0]
    # Diffusion coeff D = d(MSD)/dt * 1/2n (n = 3 = number of dimensions)
    # Ref: Carbone2014a (Carbone and Troisi)
    diffusion_coeff = numerator / (6 * denominator)
    # The error in the mobility is proportionally the same as the error in the diffusion coefficient as the
    # other variables are constants with zero error
    diff_error = diffusion_coeff * np.sqrt(
        (av_MSD_error / numerator) ** 2 + (av_time_error / denominator) ** 2
    )
    # Use Einstein-Smoluchowski relation
    mobility = (
        elementary_charge * diffusion_coeff / (kB * temperature)
    )  # This is in m^{2} / Vs
    # Convert to cm^{2}/ Vs
    mobility *= 100 ** 2
    mob_error = (diff_error / diffusion_coeff) * mobility
    return mobility, mob_error


def plot_MSD(
    times, MSDs, time_standard_errors, MSD_standard_errors, directory, carrier_type
):
    carrier_types = ["hole", "electron"]
    fit_X = np.linspace(np.min(times), np.max(times), 100)
    gradient, intercept, r_val, p_val, std_err = scipy.stats.linregress(times, MSDs)
    print("Standard Error", std_err)
    print("Fitting r_val =", r_val)
    fit_Y = (fit_X * gradient) + intercept
    mobility, mob_error = calc_mobility(
        fit_X, fit_Y, np.average(time_standard_errors), np.average(MSD_standard_errors)
    )
    plt.plot(times, MSDs)
    plt.errorbar(times, MSDs, xerr=time_standard_errors, yerr=MSD_standard_errors)
    plt.plot(fit_X, fit_Y, "r")
    plt.xlabel("Time (s)")
    plt.ylabel(r"MSD (m$^{2}$)")
    plt.title(
        "".join(
            [
                r"$\mu_{{{0}, {1}}}$ = {2:.3e} cm$^{{{3}}}$/Vs".format(
                    0, carrier_type[0], mobility, 2
                )
            ]
        ),
        y=1.1,
    )
    # 18 for hole linear MSD, 19 for electron linear MSD
    file_name = "".join(
        [
            "{:02}_lin_MSD_".format(18 + carrier_types.index(carrier_type)),
            carrier_type,
            ".png",
        ]
    )
    plt.savefig(os.path.join(directory, "figures", file_name), dpi=300)
    plt.clf()
    print("Figure saved as", os.path.join(directory, "figures", file_name))
    plt.semilogx(times, MSDs)
    plt.errorbar(times, MSDs, xerr=time_standard_errors, yerr=MSD_standard_errors)
    plt.semilogx(fit_X, fit_Y, "r")
    plt.xlabel("Time (s)")
    plt.ylabel(r"MSD (m$^{2}$)")
    plt.title(
        "".join(
            [
                r"$\mu_{{{0}, {1}}}$ = {2:.3e} cm$^{{{3}}}$/Vs".format(
                    0, carrier_type[0], mobility, 2
                )
            ]
        ),
        y=1.1,
    )
    # 20 for hole semilog MSD, 21 for electron semilog MSD
    file_name = "".join(
        [
            "{:02}_semi_log_MSD_".format(20 + carrier_types.index(carrier_type)),
            carrier_type,
            ".png",
        ]
    )
    plt.savefig(os.path.join(directory, "figures", file_name), dpi=300)
    plt.clf()
    print("Figure saved as", os.path.join(directory, "figures", file_name))
    plt.plot(times, MSDs)
    plt.errorbar(times, MSDs, xerr=time_standard_errors, yerr=MSD_standard_errors)
    plt.plot(fit_X, fit_Y, "r")
    plt.xlabel("Time (s)")
    plt.ylabel(r"MSD (m$^{2}$)")
    plt.xscale("log")
    plt.yscale("log")
    plt.title(
        "".join(
            [
                r"$\mu_{{{0}, {1}}}$ = {2:.3e} cm$^{{{3}}}$/Vs".format(
                    0, carrier_type[0], mobility, 2
                )
            ]
        ),
        y=1.1,
    )
    # 22 for hole semilog MSD, 23 for electron semilog MSD
    file_name = "".join(
        [
            "{:02}_log_MSD_".format(22 + carrier_types.index(carrier_type)),
            carrier_type,
            ".png",
        ]
    )
    plt.savefig(os.path.join(directory, "figures", file_name), dpi=300)
    plt.clf()
    print("Figure saved as", os.path.join(directory, "figures", file_name))
    return mobility, mob_error, r_val ** 2


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
    anisotropy = (3 / 2) * (
        ((diagonal[0] ** 2) + (diagonal[1] ** 2) + (diagonal[2] ** 2))
        / ((diagonal[0] + diagonal[1] + diagonal[2]) ** 2)
    ) - (1 / 2)
    return anisotropy


def plot_hop_vectors(
    carrier_data,
    AA_morphology_dict,
    chromophore_list,
    directory,
    sim_dims,
    carrier_type,
    plot3D_graphs,
):
    import matplotlib.colors
    import matplotlib.cm as cmx

    if not plot3D_graphs:
        return

    box_lengths = np.array([axis[1] - axis[0] for axis in sim_dims])
    carrier_history = carrier_data["carrier_history_matrix"]
    non_zero_coords = list(zip(*carrier_history.nonzero()))
    hop_vectors = []
    intensities = []
    for (chromo1, chromo2) in non_zero_coords:
        # Find the normal of chromo1 plane (first three atoms)
        triplet_posn = [
            np.array(AA_morphology_dict["position"][AAID])
            for AAID in chromophore_list[chromo1].AAIDs[:3]
        ]
        AB = triplet_posn[1] - triplet_posn[0]
        AC = triplet_posn[2] - triplet_posn[0]
        normal = np.cross(AB, AC)
        normal /= np.linalg.norm(normal)
        # Calculate rotation matrix required to map this onto baseline (0, 0, 1)
        rotation_matrix = hf.get_rotation_matrix(normal, np.array([0, 0, 1]))
        # Find the vector from chromo1 to chromo2 (using the correct rleative image)
        neighbour_IDs = [
            neighbour[0] for neighbour in chromophore_list[chromo1].neighbours
        ]
        nlist_index = neighbour_IDs.index(chromo2)
        relative_image = np.array(chromophore_list[chromo1].neighbours[nlist_index][1])
        chromo1_posn = chromophore_list[chromo1].posn
        chromo2_posn = chromophore_list[chromo2].posn + (relative_image * box_lengths)
        unrotated_hop_vector = chromo2_posn - chromo1_posn
        # Apply rotation matrix to chromo1-chromo2 vector
        hop_vector = np.matmul(rotation_matrix, unrotated_hop_vector)
        # Store the hop vector and the intensity
        intensity = np.log10(carrier_history[chromo1, chromo2])
        hop_vectors.append(hop_vector)
        intensities.append(intensity)

    # Convert hop_vectors to an np.array
    hop_vectors = np.array(hop_vectors)
    # Create a figure class
    fig = plt.figure(figsize=(7, 6))
    # Make a 3D subplot
    ax = fig.add_subplot(111, projection="3d")
    # Determine the smallest, non-zero number of times two chromophores are connected.
    I_min = np.min(intensities)
    # Determine the max number of times two chormophores are connected.
    I_max = np.max(intensities)

    # Set up the color bar.
    colour_map = plt.get_cmap("Greys")
    c_norm = matplotlib.colors.Normalize(vmin=np.floor(I_min), vmax=np.ceil(I_max))
    scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=colour_map)
    hop_colors = scalar_map.to_rgba(intensities)

    # # Set up the intensity for the hops so more travelled paths are more intense
    # alphas = intensities / I_max
    # hop_colors[:, 3] = alphas

    # Plot the vectors between two chromophores
    ax.quiver(
        # X-vals for start of each arrow
        np.array([0] * len(hop_colors)),
        # Y-vals for start of each arrow
        np.array([0] * len(hop_colors)),
        # Z-vals for start of each arrow
        np.array([0] * len(hop_colors)),
        # X-vecs for direction of each arrow
        hop_vectors[:, 0],
        # Y-vecs for direction of each arrow
        hop_vectors[:, 1],
        # Z-vecs for direction of each arrow
        hop_vectors[:, 2],
        color=hop_colors,
        arrow_length_ratio=0,
        linewidth=0.7,
    )
    ax.set_xlim([-np.max(hop_vectors[:, 0]), np.max(hop_vectors[:, 0])])
    ax.set_ylim([-np.max(hop_vectors[:, 1]), np.max(hop_vectors[:, 1])])
    ax.set_zlim([-np.max(hop_vectors[:, 2]), np.max(hop_vectors[:, 2])])

    # Make the colour bar
    scalar_map.set_array(intensities)
    # tick_location = np.arange(0, int(np.ceil(vmax)) + 1, 1)
    tick_location = np.arange(0, I_max, 1)
    # Plot the colour bar
    cbar = plt.colorbar(
        scalar_map, ticks=tick_location, shrink=0.8, aspect=20
    )  # , fraction=0.046, pad=0.04)
    cbar.ax.set_yticklabels([r"10$^{{{}}}$".format(int(x)) for x in tick_location])

    # Name and save the figure.
    carrier_types = ["hole", "electron"]
    materials_types = ["donor", "acceptor"]
    carrier_index = carrier_types.index(carrier_type)
    figure_title = " ".join(
        [
            materials_types[carrier_index].capitalize(),
            "".join(["(", carrier_type.capitalize(), ")"]),
            "Network",
        ]
    )
    # plt.title(figure_title, y=1.1)
    # 36 for donor 3d network, 37 for acceptor 3d network
    file_name = "".join(
        ["{:02}_hop_vec_".format(36 + carrier_index), carrier_type, ".png"]
    )
    plt.savefig(
        os.path.join(directory, "figures", file_name), bbox_inches="tight", dpi=300
    )
    print("Figure saved as", os.path.join(directory, "figures", file_name))
    plt.clf()


def plot_anisotropy(carrier_data, directory, sim_dims, carrier_type, plot3D_graphs):
    sim_extent = [value[1] - value[0] for value in sim_dims]
    xvals = []
    yvals = []
    zvals = []
    colours = []
    sim_dims_nm = list(map(list, np.array(sim_dims) / 10.))
    # Get the indices of the carriers that travelled the furthest
    if len(carrier_data["final_position"]) <= 1000:
        carrier_indices_to_use = range(len(carrier_data["final_position"]))
    else:
        displacements = copy.deepcopy(np.array(carrier_data["displacement"]))
        carrier_indices_to_use = displacements.argsort()[-1000:][::-1]
    for carrier_no in carrier_indices_to_use:
        posn = carrier_data["final_position"][carrier_no]
        # if bool(sum([x < -3 or x > 3 for x in image])):
        #     continue
        position = [0.0, 0.0, 0.0]
        for axis in range(len(posn)):
            position[axis] = (
                carrier_data["image"][carrier_no][axis] * sim_extent[axis]
            ) + posn[axis]
        xvals.append(position[0] / 10.)
        yvals.append(position[1] / 10.)
        zvals.append(position[2] / 10.)
        colours.append("b")
    anisotropy = calculate_anisotropy(xvals, yvals, zvals)
    if not plot3D_graphs:
        return anisotropy
    print("----------====================----------")
    print(
        "{0:s} charge transport anisotropy calculated as {1:f}".format(
            carrier_type.capitalize(), anisotropy
        )
    )
    print("----------====================----------")
    # Reduce number of plot markers
    fig = plt.gcf()
    ax = p3.Axes3D(fig)
    if len(xvals) > 1000:
        xvals = xvals[0 : len(xvals) : len(xvals) // 1000]
        yvals = yvals[0 : len(yvals) : len(yvals) // 1000]
        zvals = zvals[0 : len(zvals) : len(zvals) // 1000]
    plt.scatter(xvals, yvals, zs=zvals, c=colours, s=20)
    plt.scatter(0, 0, zs=0, c="r", s=50)
    # Draw boxlines
    # Varying X
    ax.plot(
        [sim_dims_nm[0][0], sim_dims_nm[0][1]],
        [sim_dims_nm[1][0], sim_dims_nm[1][0]],
        [sim_dims_nm[2][0], sim_dims_nm[2][0]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims_nm[0][0], sim_dims_nm[0][1]],
        [sim_dims_nm[1][1], sim_dims_nm[1][1]],
        [sim_dims_nm[2][0], sim_dims_nm[2][0]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims_nm[0][0], sim_dims_nm[0][1]],
        [sim_dims_nm[1][0], sim_dims_nm[1][0]],
        [sim_dims_nm[2][1], sim_dims_nm[2][1]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims_nm[0][0], sim_dims_nm[0][1]],
        [sim_dims_nm[1][1], sim_dims_nm[1][1]],
        [sim_dims_nm[2][1], sim_dims_nm[2][1]],
        c="k",
        linewidth=1.0,
    )
    # Varying Y
    ax.plot(
        [sim_dims_nm[0][0], sim_dims_nm[0][0]],
        [sim_dims_nm[1][0], sim_dims_nm[1][1]],
        [sim_dims_nm[2][0], sim_dims_nm[2][0]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims_nm[0][1], sim_dims_nm[0][1]],
        [sim_dims_nm[1][0], sim_dims_nm[1][1]],
        [sim_dims_nm[2][0], sim_dims_nm[2][0]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims_nm[0][0], sim_dims_nm[0][0]],
        [sim_dims_nm[1][0], sim_dims_nm[1][1]],
        [sim_dims_nm[2][1], sim_dims_nm[2][1]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims_nm[0][1], sim_dims_nm[0][1]],
        [sim_dims_nm[1][0], sim_dims_nm[1][1]],
        [sim_dims_nm[2][1], sim_dims_nm[2][1]],
        c="k",
        linewidth=1.0,
    )
    # Varying Z
    ax.plot(
        [sim_dims_nm[0][0], sim_dims_nm[0][0]],
        [sim_dims_nm[1][0], sim_dims_nm[1][0]],
        [sim_dims_nm[2][0], sim_dims_nm[2][1]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims_nm[0][0], sim_dims_nm[0][0]],
        [sim_dims_nm[1][1], sim_dims_nm[1][1]],
        [sim_dims_nm[2][0], sim_dims_nm[2][1]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims_nm[0][1], sim_dims_nm[0][1]],
        [sim_dims_nm[1][0], sim_dims_nm[1][0]],
        [sim_dims_nm[2][0], sim_dims_nm[2][1]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims_nm[0][1], sim_dims_nm[0][1]],
        [sim_dims_nm[1][1], sim_dims_nm[1][1]],
        [sim_dims_nm[2][0], sim_dims_nm[2][1]],
        c="k",
        linewidth=1.0,
    )
    ax.set_xlabel("X (nm)", fontsize=20, labelpad=40)
    ax.set_ylabel("Y (nm)", fontsize=20, labelpad=40)
    ax.set_zlabel("Z (nm)", fontsize=20, labelpad=40)
    maximum = max([max(xvals), max(yvals), max(zvals)])
    ax.set_xlim([-maximum, maximum])
    ax.set_ylim([-maximum, maximum])
    ax.set_zlim([-maximum, maximum])
    for tick in (
        ax.xaxis.get_major_ticks()
        + ax.yaxis.get_major_ticks()
        + ax.zaxis.get_major_ticks()
    ):
        tick.label.set_fontsize(16)
    ax.dist = 11
    # 08 for hole anisotropy, 09 for electron anisotropy
    if carrier_type == "hole":
        figure_index = "08"
    elif carrier_type == "electron":
        figure_index = "09"
    plt.title("Anisotropy ({:s})".format(carrier_type.capitalize()), y=1.1)
    file_name = os.path.join(
        directory,
        "figures",
        "".join([figure_index, "_anisotropy_", carrier_type, ".png"]),
    )
    plt.savefig(file_name, bbox_inches="tight", dpi=300)
    plt.clf()
    print("Figure saved as", file_name)
    return anisotropy


def get_temp_val(string):
    return string.split("-")[-2][2:]


def get_frame_val(string):
    return string.split("-")[0][1:]


def plot_temperature_progression(
    temp_data, mobility_data, anisotropy_data, carrier_type, x_label
):
    plt.gcf()
    plt.clf()
    xvals = temp_data
    yvals = list(np.array(mobility_data)[:, 0])
    yerrs = list(np.array(mobility_data)[:, 1])
    plt.xlabel(x_label)
    plt.ylabel(r"Mobility (cm$^{2}$ / Vs)")
    plt.errorbar(xvals, yvals, xerr=0, yerr=yerrs)
    plt.yscale("log")
    file_name = "".join(["mobility_", carrier_type, ".png"])
    plt.savefig(file_name, dpi=300)
    plt.clf()
    print("Figure saved as {:s}".format(file_name))

    plt.plot(temp_data, anisotropy_data, c="r")
    plt.xlabel(x_label)
    plt.ylabel(r"$\kappa$ (Arb. U)")
    file_name = "".join(["anisotropy_", carrier_type, ".png"])
    plt.savefig(file_name, dpi=300)
    plt.clf()
    print("Figure saved as {:s}".format(file_name))


def calculate_lambda_ij(chromo_length):
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
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def gauss_fit(data):
    mean = np.mean(data)
    std = np.std(data)
    hist, bin_edges = np.histogram(data, bins=100)
    try:
        fit_args, fit_conv = curve_fit(
            gaussian, bin_edges[:-1], hist, p0=[1, mean, std]
        )
    except RuntimeError:
        fit_args = None
    return bin_edges, fit_args, mean, std


def split_molecules(input_dictionary):
    # Split the full morphology into individual molecules
    molecule_AAIDs = []
    molecule_lengths = []
    # Create a lookup table `neighbour list' for all connected atoms called {bondedAtoms}
    bonded_atoms = hf.obtain_bonded_list(input_dictionary["bond"])
    molecule_list = [i for i in range(len(input_dictionary["type"]))]
    # Recursively add all atoms in the neighbour list to this molecule
    for mol_ID in range(len(molecule_list)):
        molecule_list = update_molecule(mol_ID, molecule_list, bonded_atoms)
    # Create a dictionary of the molecule data
    molecule_data = {}
    for atom_ID in range(len(input_dictionary["type"])):
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
                molecule_list = update_molecule(
                    bonded_atom, molecule_list, bonded_atoms
                )
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


def plot_neighbour_hist(
    chromophore_list,
    chromo_to_mol_ID,
    morphology_shape,
    output_dir,
    sep_cut_donor,
    sep_cut_acceptor,
):
    separation_dist_donor = []
    separation_dist_acceptor = []
    for chromo1 in chromophore_list:
        for chromo2_details in chromo1.neighbours:
            if (chromo2_details is None) or (
                chromo1.ID == chromophore_list[chromo2_details[0]].ID
            ):
                continue
            chromo2 = chromophore_list[chromo2_details[0]]
            # Skip any chromophores that are part of the same molecule
            if chromo_to_mol_ID[chromo1.ID] == chromo_to_mol_ID[chromo2.ID]:
                continue
            separation = np.linalg.norm(
                (
                    np.array(chromo2.posn)
                    + (np.array(chromo2_details[1]) * np.array(morphology_shape))
                )
                - chromo1.posn
            )
            if chromo1.species == "donor":
                separation_dist_donor.append(separation)
            elif chromo1.species == "acceptor":
                separation_dist_acceptor.append(separation)
    material = ["donor", "acceptor"]
    sep_cuts = [sep_cut_donor, sep_cut_acceptor]
    for material_type, separation_dist in enumerate(
        [separation_dist_donor, separation_dist_acceptor]
    ):
        if len(separation_dist) == 0:
            continue
        plt.figure()
        n, bin_edges, _ = plt.hist(separation_dist, bins=40, color="b")
        bin_centres = (bin_edges[1:] + bin_edges[:-1]) / 2.0
        bin_centres = np.insert(bin_centres, 0, 0)
        n = np.insert(n, 0, 0)
        smoothed_n = gaussian_filter(n, 1.0)
        plt.plot(bin_centres, smoothed_n, color="r")
        if (sep_cuts[material_type] is not None) and (
            sep_cuts[material_type].lower() == "auto"
        ):
            sep_cuts[material_type] = calculate_cut_off_from_dist(
                bin_centres,
                smoothed_n,
                minimum_index=0,
                value_at_least=100,
                logarithmic=False,
            )
        if sep_cuts[material_type] is not None:
            print(
                "Cluster cut-off based on",
                material[material_type],
                "chromophore separation set to",
                sep_cuts[material_type],
            )
            plt.axvline(float(sep_cuts[material_type]), c="k")
        plt.xlabel("".join([material[material_type].capitalize(), r" r$_{i,j}$ (\AA)"]))
        plt.ylabel("Frequency (Arb. U.)")
        # 04 for donor neighbour hist, 05 for acceptor neighbour hist
        file_name = "".join(
            [
                "{:02}_neighbour_hist_".format(4 + material_type),
                material[material_type].lower(),
                ".png",
            ]
        )
        plt.savefig(os.path.join(output_dir, file_name), dpi=300)
        plt.close()
        print(
            "Neighbour histogram figure saved as", os.path.join(output_dir, file_name)
        )
    return sep_cuts[0], sep_cuts[1]


def plot_orientation_hist(
    chromophore_list,
    chromo_to_mol_ID,
    orientations_data,
    output_dir,
    o_cut_donor,
    o_cut_acceptor,
):
    orientation_dist_donor = []
    orientation_dist_acceptor = []
    for chromo1 in chromophore_list:
        orientation_1 = orientations_data[chromo1.ID]
        for chromo2_details in chromo1.neighbours:
            if (chromo2_details is None) or (
                chromo1.ID == chromophore_list[chromo2_details[0]].ID
            ):
                continue
            chromo2 = chromophore_list[chromo2_details[0]]
            # Skip any chromophores that are part of the same molecule
            if chromo_to_mol_ID[chromo1.ID] == chromo_to_mol_ID[chromo2.ID]:
                continue
            orientation_2 = orientations_data[chromo2_details[0]]
            dot_product = np.dot(orientation_1, orientation_2)
            separation_angle = np.arccos(np.abs(dot_product)) * 180 / np.pi
            if chromo1.species == "donor":
                orientation_dist_donor.append(separation_angle)
            elif chromo1.species == "acceptor":
                orientation_dist_acceptor.append(separation_angle)
    material = ["donor", "acceptor"]
    o_cuts = [o_cut_donor, o_cut_acceptor]
    for material_type, orientation_dist in enumerate(
        [orientation_dist_donor, orientation_dist_acceptor]
    ):
        if len(orientation_dist) == 0:
            continue
        plt.figure()
        n, bin_edges, _ = plt.hist(orientation_dist, bins=40, color="b")
        bin_centres = (bin_edges[1:] + bin_edges[:-1]) / 2.0
        bin_centres = np.insert(bin_centres, 0, 0)
        n = np.insert(n, 0, 0)
        smoothed_n = gaussian_filter(n, 1.0)
        plt.plot(bin_centres, smoothed_n, color="r")
        if (o_cuts[material_type] is not None) and (
            o_cuts[material_type].lower() == "auto"
        ):
            o_cuts[material_type] = calculate_cut_off_from_dist(
                bin_centres,
                smoothed_n,
                maximum_index=0,
                value_at_least=100,
                logarithmic=False,
            )
        if o_cuts[material_type] is not None:
            print(
                "Cluster cut-off based on",
                material[material_type],
                "relative chromophore orientations set to",
                o_cuts[material_type],
            )
            plt.axvline(float(o_cuts[material_type]), c="k")
            plt.xlabel(
                r"{:s} Orientations (Deg)".format(material[material_type].capitalize())
            )
        plt.xlim([0, 90])
        plt.xticks(np.arange(0, 91, 15))
        plt.ylabel("Frequency (Arb. U.)")
        # 34 for donor neighbour hist, 35 for acceptor neighbour hist
        file_name = "".join(
            [
                "{:02}_orientation_hist_".format(34 + material_type),
                material[material_type].lower(),
                ".png",
            ]
        )
        plt.savefig(os.path.join(output_dir, file_name), dpi=300)
        plt.close()
        print(
            "Orientation histogram figure saved as", os.path.join(output_dir, file_name)
        )
    return o_cuts[0], o_cuts[1]


def create_cut_off_dict(
    sep_cut_donor,
    sep_cut_acceptor,
    o_cut_donor,
    o_cut_acceptor,
    ti_cut_donor,
    ti_cut_acceptor,
    freq_cut_donor,
    freq_cut_acceptor,
):
    cut_off_dict = {
        "separation": [sep_cut_donor, sep_cut_acceptor],
        "orientation": [o_cut_donor, o_cut_acceptor],
        "TI": [ti_cut_donor, ti_cut_acceptor],
        "freq": [freq_cut_donor, freq_cut_acceptor],
    }
    # Convert all cuts to floats (leaving NoneType intact)
    for cut_type, cuts in cut_off_dict.items():
        for material_index, cut in enumerate(cuts):
            try:
                cut_off_dict[cut_type][material_index] = float(cut)
            except TypeError:
                pass
            except ValueError:
                # Most likely both donor and acceptor cuts have been specified
                # for this KMCA job, but only one type of carrier is available in
                # this instance, so we can safely set the opposing carrier type cut
                # to None.
                cut_off_dict[cut_type][material_index] = None
    return cut_off_dict


def get_clusters(
    chromophore_list,
    carrier_history_dict,
    orientations_all,
    cut_off_dict,
    AA_morphology_dict,
):
    freud_box = [[AA_morphology_dict[coord] for coord in ["lx", "ly", "lz"]]]
    materials_to_check = ["donor", "acceptor"]
    carriers_to_check = ["hole", "electron"]
    cluster_freqs = [{}, {}]
    cluster_dicts = [{}, {}]
    clusters_total = [0, 0]
    clusters_large = [0, 0]
    clusters_biggest = [0, 0]
    species_psi = [0, 0]
    clusters_cutoffs = [[], []]
    if any(cut is not None for cut in cut_off_dict["separation"]):
        separations = get_separations(chromophore_list, AA_morphology_dict)
    else:
        separations = None
    for type_index, material_type in enumerate(materials_to_check):
        material_chromos = [
            chromo for chromo in chromophore_list if chromo.species == material_type
        ]
        print("Examining the", material_type, "material...")
        positions = np.array([chromo.posn for chromo in material_chromos])
        if len(positions) == 0:
            print("No material found. Continuing...")
            continue
        print("Obtaining orientations of each chromophore...")
        if cut_off_dict["orientation"][type_index] is not None:
            orientations = [orientations_all[chromo.ID] for chromo in material_chromos]
        else:
            orientations = None
        print("Calculating clusters...")
        n_list = get_n_list(
            chromophore_list,
            carrier_history_dict[carriers_to_check[type_index]],
            separations=separations,
            r_cut=cut_off_dict["separation"][type_index],
            orientations=orientations,
            o_cut=cut_off_dict["orientation"][type_index],
            ti_cut=cut_off_dict["TI"][type_index],
            freq_cut=cut_off_dict["freq"][type_index],
            sim_dims=freud_box,
        )
        clusters_list = make_clusters(n_list)
        cluster_dict = {}
        for chromo_ID, cluster_ID in enumerate(clusters_list):
            cluster_dict[chromo_ID] = cluster_ID
        cluster_freq = {}
        for cluster_ID in set(clusters_list):
            cluster_freq[cluster_ID] = clusters_list.count(cluster_ID)
        species_psi[type_index] = sum(
            [pair for key, pair in cluster_freq.items() if pair > 6]
        ) / len(clusters_list)
        clusters_total[type_index] = len([key for key, val in cluster_freq.items()])
        clusters_large[type_index] = len(
            [key for key, val in cluster_freq.items() if val > 6]
        )
        clusters_biggest[type_index] = np.max(list(cluster_freq.values()))
        clusters_cutoffs[type_index] = [
            cut_off_dict[cut_type][type_index]
            for cut_type in ["separation", "orientation", "TI", "freq"]
        ]
        print("----------====================----------")
        print(
            "Detected", clusters_total[type_index], material_type, "clusters in total."
        )
        print(
            "Detected",
            clusters_large[type_index],
            material_type,
            "clusters with size > 6.",
        )
        print("Largest cluster size =", clusters_biggest[type_index], "chromophores.")
        print(
            'Ratio of chromophores in "large" clusters (psi): {:.3f}'.format(
                species_psi[type_index]
            )
        )
        print("----------====================----------")
        cluster_dicts[type_index] = cluster_dict
        cluster_freqs[type_index] = cluster_freq
    return (
        cluster_dicts,
        cluster_freqs,
        clusters_total,
        clusters_large,
        clusters_biggest,
        clusters_cutoffs,
    )


def make_clusters(n_list):
    """
    Function to call for the creation
    of turning the neighbor list into a
    cluster list.
    Requires:
        n_list - neighbor list
    Returns:
        c_list - cluster list
    """
    sys.setrecursionlimit(int(5e4))
    c_list = [i for i in range(len(n_list))]
    for i in range(len(c_list)):
        n_list, c_list = update_neighbors(i, c_list, n_list)
    return c_list


def update_neighbors(particle, cluster_list, neighbor_list):
    """Recursive function to convert neighborlist into cluster list"""
    for n in neighbor_list[particle]:
        if cluster_list[n] > cluster_list[particle]:
            cluster_list[n] = cluster_list[particle]
            neighbor_list, cluster_list = update_neighbors(
                n, cluster_list, neighbor_list
            )
        elif cluster_list[n] < cluster_list[particle]:
            cluster_list[particle] = cluster_list[n]
            neighbor_list, cluster_list = update_neighbors(
                particle, cluster_list, neighbor_list
            )
    return neighbor_list, cluster_list


def get_n_list(
    chromophore_list,
    carrier_history,
    separations,
    r_cut,
    orientations,
    o_cut,
    ti_cut,
    freq_cut,
    sim_dims,
):
    if o_cut is not None:
        # o_cut is currently an angle in degrees.
        # Need to calculate the corresponding dot-product cut off for this
        # angle
        dotcut = np.cos(o_cut * np.pi / 180)
    n_list = []
    for chromophore in chromophore_list:
        n_list.append([neighbour[0] for neighbour in chromophore.neighbours])
    printing = False
    for chromo_ID, neighbours in enumerate(n_list):
        # if chromo_ID == 47:
        #    printing = True
        remove_list = []
        if printing is True:
            print(chromo_ID)
            print(chromophore_list[chromo_ID].neighbours, "==", n_list[chromo_ID])
            print(chromophore_list[chromo_ID].neighbours_TI)
        for neighbour_index, neighbour_ID in enumerate(neighbours):
            ti = chromophore_list[chromo_ID].neighbours_TI[neighbour_index]
            if printing is True:
                print(
                    "Examining neighbour_index",
                    neighbour_index,
                    "which corresponds to chromo ID",
                    neighbour_ID,
                )
                print("TI =", ti)
            # Check the cut_offs in order of ease of calculation
            # Simple TI lookup first
            if (ti_cut is not None) and (ti < ti_cut):
                remove_list.append(neighbour_ID)
                if printing is True:
                    print("Adding", neighbour_ID, "to remove list as TI =", ti)
                    print("Remove list now =", remove_list)
                continue
            # Simple hop_frequency cutoff lookup in the carrier history
            if (freq_cut is not None) and (carrier_history is not None):
                total_hops = (
                    carrier_history[chromo_ID, neighbour_ID]
                    + carrier_history[neighbour_ID, chromo_ID]
                )
                if total_hops < freq_cut:
                    remove_list.append(neighbour_ID)
                    continue
            # Separation cutoff lookup in the separations matrix depends on the chromo_IDs
            if (r_cut is not None) and (separations is not None):
                if chromo_ID < neighbour_ID:
                    separation = separations[chromo_ID, neighbour_ID]
                else:
                    separation = separations[neighbour_ID, chromo_ID]
                if separation > r_cut:
                    remove_list.append(neighbour_ID)
                    continue
            # Some dot product manipulation is required to get the orientations right
            if (o_cut is not None) and (orientations is not None):
                chromo1_normal = orientations[chromo_ID]
                chromo2_normal = orientations[neighbour_ID]
                rotation_dot_product = abs(np.dot(chromo1_normal, chromo2_normal))
                if rotation_dot_product < dotcut:
                    remove_list.append(neighbour_ID)
                    if printing is True:
                        print(
                            "Adding",
                            neighbour_ID,
                            "to remove list as dot_product =",
                            rotation_dot_product,
                        )
                        print("Remove list now =", remove_list)
                    continue
        if printing is True:
            print("n_list_current =", n_list[chromo_ID])
            print("Remove list final =", remove_list)
        for neighbour_ID in remove_list:
            n_list[chromo_ID].remove(neighbour_ID)
            ##Will need to check the reverse remove.
            # n_list[neighbour_ID].remove(chromo_ID)
        if printing is True:
            print("n_list_after remove =", n_list[chromo_ID])
            exit()
        # input("".join(map(str, [chromo_ID, n_list[chromo_ID]])))
    return n_list


def get_orientations(
    chromophore_list,
    CG_morphology_dict,
    AA_morphology_dict,
    CG_to_AAID_master,
    parameter_dict,
):
    orientations = []
    for index, chromophore in enumerate(chromophore_list):
        positions = get_electronic_atom_positions(
            chromophore,
            CG_morphology_dict,
            AA_morphology_dict,
            CG_to_AAID_master,
            parameter_dict,
        )
        # There is a really cool way to do this with single value composition
        # but on the time crunch I didn't have time to learn how to implement
        # it properly. Check https://goo.gl/jxuhvJ for more details.
        plane = calculate_plane(positions)
        orientations.append(np.array(plane) / np.linalg.norm(plane))

        # colours = {"S": "y", "CA": "r", "H1": "b"}
        # for i in range(len(positions)):
        #    print(i, types[i], positions[i])
        # plt.figure()
        # fig = plt.gcf()
        # ax = p3.Axes3D(fig)
        # for AAID, position in enumerate(positions):
        #    plt.scatter([coord[0] for coord in positions], [coord[1] for coord in positions], zs=[coord[2] for coord in positions], facecolor=colours[types[AAID]], s=50)
        # plane = np.array(plane) / np.linalg.norm(plane)
        # com = hf.calc_COM(positions, types)
        # for i in range(20):
        #    new_coord = com + (i * plane / 10)
        #    plt.scatter(new_coord[0], new_coord[1], zs=new_coord[2], c='k', s=20)
        # plot_size = max([axis[1] - axis[0] for axis in [ax.get_xlim(), ax.get_ylim(), ax.get_zlim()]])
        # plot_av = [np.mean([axis[1], axis[0]]) for axis in [ax.get_xlim(), ax.get_ylim(), ax.get_zlim()]]
        # lim = [[mean - (plot_size / 2.0), mean + (plot_size / 2.0)] for mean in plot_av]
        # ax.set_xlim(lim[0])
        # ax.set_ylim(lim[1])
        # ax.set_zlim(lim[2])
        # plt.show()
    return orientations


def calculate_plane(positions):

    # tmp_A = []
    # tmp_b = []
    # for i in range(len(positions)):
    #    tmp_A.append([positions[i][0], positions[i][1], 1])
    #    tmp_b.append(positions[i][2])
    # b = np.matrix(tmp_b).T
    # A = np.matrix(tmp_A)
    # normal_vec = (A.T * A).I * A.T * b

    ## See https://goo.gl/jxuhvJ for details on this methodology.
    # LHS = np.matrix([[x, y, 1] for [x, y, _] in positions])
    # RHS = np.matrix([[z] for [_, _, z] in positions])
    # normal_vec = np.linalg.inv(LHS.T * LHS) * LHS.T * RHS
    # return np.array([value for val in normal_vec.tolist() for value in val])

    vec1 = hf.find_axis(positions[0], positions[1])
    vec2 = hf.find_axis(positions[0], positions[2])
    return np.cross(vec1, vec2)


def get_separations(chromophore_list, AA_morphology_dict):
    box_dims = [AA_morphology_dict[axis] for axis in ["lx", "ly", "lz"]]
    separations = lil_matrix(
        (len(chromophore_list), len(chromophore_list)), dtype=float
    )
    for chromo_ID, chromophore in enumerate(chromophore_list):
        for neighbour_details in chromophore.neighbours:
            neighbour_ID = neighbour_details[0]
            relative_image = neighbour_details[1]
            neighbour_posn = chromophore_list[neighbour_ID].posn
            neighbour_chromo_posn = neighbour_posn + (
                np.array(relative_image) * np.array(box_dims)
            )
            separation = hf.calculate_separation(
                chromophore.posn, neighbour_chromo_posn
            )
            if chromo_ID < neighbour_ID:
                separations[chromo_ID, neighbour_ID] = separation
            else:
                separations[neighbour_ID, chromo_ID] = separation
    return separations


def get_electronic_atom_positions(
    chromophore,
    CG_morphology_dict,
    AA_morphology_dict,
    CG_to_AAID_master,
    parameter_dict,
):
    # We don't save this in the chromophore info so we'll have to calculate it
    # again.
    # Determine whether this chromophore is a donor or an acceptor, as well
    # as the site types that have been defined as the electronically active
    # in the chromophore
    if CG_morphology_dict is not None:
        # Normal operation
        CG_types = sorted(
            list(set([CG_morphology_dict["type"][CGID] for CGID in chromophore.CGIDs]))
        )
        active_CG_sites, _ = obtain_electronic_species(
            chromophore.CGIDs,
            CG_morphology_dict["type"],
            parameter_dict["CG_site_species"],
        )
        # CG_to_AAID_master is a list of dictionaries where each list
        # element corresponds to a new molecule. Firstly, flatten this out
        # so that it becomes a single CG:AAID dictionary
        flattened_CG_to_AAID_master = {
            dict_key: dict_val[1]
            for dictionary in CG_to_AAID_master
            for dict_key, dict_val in dictionary.items()
        }
        # By using active_CG_sites, determine the AAIDs for
        # the electrically active proportion of the chromophore, so that we
        # can calculate its proper position. Again each element corresponds
        # to each CG site so the list needs to be flattened afterwards.
        electronically_active_AAIDs = [
            AAID
            for AAIDs in [flattened_CG_to_AAID_master[CGID] for CGID in active_CG_sites]
            for AAID in AAIDs
        ]
    else:
        # No fine-graining has been performed by MorphCT, so we know that
        # the input morphology is already atomistic.
        if len(parameter_dict["CG_site_species"]) == 1:
            # If the morphology contains only a single type of electronic
            # species, then the parameter_dict['CG_site_species'] should
            # only have one entry, and we can set all chromophores to be
            # this species.
            active_CG_sites = chromophore.CGIDs
            electronically_active_AAIDs = chromophore.CGIDs
        elif (len(parameter_dict["CG_site_species"]) == 0) and (
            len(parameter_dict["AA_rigid_body_species"]) > 0
        ):
            # If the CG_site_species have not been specified, then look to
            # the AA_rigid_body_species dictionary to determine which rigid
            # bodies are donors and which are acceptors
            electronically_active_AAIDs = []
            for AAID in chromophore.CGIDs:
                if AA_morphology_dict["body"][AAID] != -1:
                    electronically_active_AAIDs.append(AAID)
        else:
            raise SystemError(
                "Multiple electronic species defined, but no way to map them"
                " without a coarse-grained morphology (no CG morph has been given)"
            )
    # The position of the chromophore can be calculated easily. Note that
    # here, the `self.image' is the periodic image that the
    # unwrapped_position of the chromophore is located in, relative to the
    # original simulation volume.
    electronically_active_unwrapped_posns = [
        AA_morphology_dict["unwrapped_position"][AAID]
        for AAID in electronically_active_AAIDs
    ]
    return electronically_active_unwrapped_posns


def obtain_electronic_species(chromophore_CG_sites, CG_site_types, CG_to_species):
    electronically_active_sites = []
    current_chromophore_species = None
    for CG_site_ID in chromophore_CG_sites:
        site_type = CG_site_types[CG_site_ID]
        site_species = CG_to_species[site_type]
        if site_species.lower() != "none":
            if (current_chromophore_species is not None) and (
                current_chromophore_species != site_species
            ):
                raise SystemError(
                    "Problem - multiple electronic species defined in the same "
                    " chromophore. Please modify the chromophore generation code "
                    " to fix this issue for your molecule!"
                )
            else:
                current_chromophore_species = site_species
                electronically_active_sites.append(CG_site_ID)
    return electronically_active_sites, current_chromophore_species


def update_cluster(atom_ID, cluster_list, neighbour_dict):
    try:
        for neighbour in neighbour_dict[atom_ID]:
            if cluster_list[neighbour] > cluster_list[atom_ID]:
                cluster_list[neighbour] = cluster_list[atom_ID]
                cluster_list = update_cluster(neighbour, cluster_list, neighbour_dict)
            elif cluster_list[neighbour] < cluster_list[atom_ID]:
                cluster_list[atom_ID] = cluster_list[neighbour]
                cluster_list = update_cluster(neighbour, cluster_list, neighbour_dict)
    except KeyError:
        pass
    return cluster_list


def write_cluster_tcl_script(output_dir, cluster_lookup, large_cluster):
    """
    Create a tcl script for each identified cluster.
    """
    # Obtain the IDs of the cluster sizes, sorted by largest first
    cluster_order = list(
        zip(
            *sorted(
                zip(
                    [len(val) for val in cluster_lookup.values()], cluster_lookup.keys()
                ),
                reverse=True,
            )
        )
    )[1]
    colors = list(range(int(1e6)))
    count = 0

    tcl_text = ["mol delrep 0 0;"]
    tcl_text += ["pbc wrap -center origin;"]
    tcl_text += ["pbc box -color black -center origin -width 4;"]
    tcl_text += ["display resetview;"]
    tcl_text += ["color change rgb 9 1.0 0.29 0.5;"]
    tcl_text += ["color change rgb 16 0.25 0.25 0.25;"]

    for cluster_ID in cluster_order:
        chromos = cluster_lookup[cluster_ID]
        chromo_IDs = [chromo.ID for chromo in chromos if chromo.species == "donor"]
        if len(chromo_IDs) > large_cluster:  # Only make clusters that are ``large''
            inclust = ""
            for chromo in chromo_IDs:
                inclust = "".join([inclust, "{:d} ".format(chromo)])
            tcl_text += ["mol material AOEdgy;"]  # Use AOEdgy if donor
            # The +1 makes the largest cluster red rather than blue (looks better
            # with AO, DoF, shadows)
            tcl_text += ["mol color ColorID {:d};".format(colors[count + 1 % 32])]
            # VMD has 32 unique colors
            tcl_text += ["mol representation VDW 4.0 8.0;"]
            tcl_text += ["mol selection resid {:s};".format(inclust)]
            tcl_text += ["mol addrep 0;"]
            count += 1
        chromo_IDs = [chromo.ID for chromo in chromos if chromo.species == "acceptor"]
        if len(chromo_IDs) > large_cluster:
            inclust = ""
            for chromo in chromo_IDs:
                inclust += "".join([inclust, "{:d} ".format(chromo)])
            tcl_text += ["mol material Glass2;"]  # Use Glass2 if acceptor
            # The +1 makes the largest cluster red rather than blue (looks better
            # with AO, DoF, shadows)
            tcl_text += ["mol color ColorID {:d};".format(colors[count + 1 % 32])]
            tcl_text += ["mol representation VDW 4.0 8.0;"]
            tcl_text += ["mol selection resid {:s};".format(inclust)]
            tcl_text += ["mol addrep 0;"]
            count += 1
    tcl_file_path = os.path.join(
        output_dir.replace("figures", "morphology"), "cluster_colors.tcl"
    )
    with open(tcl_file_path, "w+") as tcl_file:
        tcl_file.writelines("".join(tcl_text))
    print("Clusters coloring written to {:s}".format(tcl_file_path))


def generate_lists_for_3d_clusters(cluster_lookup, colours, large_cluster):
    data = []
    for cluster_ID, chromos in cluster_lookup.items():
        if len(chromos) > large_cluster:
            for chromo in chromos:
                if chromo.species == "donor":
                    data.append(
                        [
                            chromo.posn[0],
                            chromo.posn[1],
                            chromo.posn[2],
                            "w",
                            colours[cluster_ID % 7],
                        ]
                    )
                if chromo.species == "acceptor":
                    data.append(
                        [
                            chromo.posn[0],
                            chromo.posn[1],
                            chromo.posn[2],
                            colours[cluster_ID % 7],
                            "none",
                        ]
                    )
    # Needs a sort so that final image is layered properly
    data = list(sorted(data, key=lambda x: x[0]))
    # Split up list into sublists
    xs = [row[0] for row in data]
    ys = [row[1] for row in data]
    zs = [row[2] for row in data]
    face_colors = [row[3] for row in data]
    edge_colors = [row[4] for row in data]
    return xs, ys, zs, face_colors, edge_colors


def plot_clusters_3D(
    output_dir, chromophore_list, cluster_dicts, sim_dims, generate_tcl
):
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    colours = ["r", "g", "b", "c", "m", "y", "k"]
    large_cluster = 6
    cluster_dict = {}
    for dictionary in cluster_dicts:
        if dictionary is not None:
            cluster_dict.update(dictionary)
    cluster_lookup = {}
    for chromo_ID, cluster_ID in cluster_dict.items():
        if cluster_ID not in cluster_lookup.keys():
            cluster_lookup[cluster_ID] = []
        else:
            cluster_lookup[cluster_ID].append(chromophore_list[chromo_ID])
    if generate_tcl:
        write_cluster_tcl_script(output_dir, cluster_lookup, large_cluster)

    xs, ys, zs, face_colors, edge_colors = generate_lists_for_3d_clusters(
        cluster_lookup, colours, large_cluster
    )
    ax.scatter(
        xs, ys, zs, facecolors=face_colors, edgecolors=edge_colors, alpha=0.6, s=40
    )
    # Draw boxlines
    # Varying X
    ax.plot(
        [sim_dims[0][0], sim_dims[0][1]],
        [sim_dims[1][0], sim_dims[1][0]],
        [sim_dims[2][0], sim_dims[2][0]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims[0][0], sim_dims[0][1]],
        [sim_dims[1][1], sim_dims[1][1]],
        [sim_dims[2][0], sim_dims[2][0]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims[0][0], sim_dims[0][1]],
        [sim_dims[1][0], sim_dims[1][0]],
        [sim_dims[2][1], sim_dims[2][1]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims[0][0], sim_dims[0][1]],
        [sim_dims[1][1], sim_dims[1][1]],
        [sim_dims[2][1], sim_dims[2][1]],
        c="k",
        linewidth=1.0,
    )
    # Varying Y
    ax.plot(
        [sim_dims[0][0], sim_dims[0][0]],
        [sim_dims[1][0], sim_dims[1][1]],
        [sim_dims[2][0], sim_dims[2][0]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims[0][1], sim_dims[0][1]],
        [sim_dims[1][0], sim_dims[1][1]],
        [sim_dims[2][0], sim_dims[2][0]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims[0][0], sim_dims[0][0]],
        [sim_dims[1][0], sim_dims[1][1]],
        [sim_dims[2][1], sim_dims[2][1]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims[0][1], sim_dims[0][1]],
        [sim_dims[1][0], sim_dims[1][1]],
        [sim_dims[2][1], sim_dims[2][1]],
        c="k",
        linewidth=1.0,
    )
    # Varying Z
    ax.plot(
        [sim_dims[0][0], sim_dims[0][0]],
        [sim_dims[1][0], sim_dims[1][0]],
        [sim_dims[2][0], sim_dims[2][1]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims[0][0], sim_dims[0][0]],
        [sim_dims[1][1], sim_dims[1][1]],
        [sim_dims[2][0], sim_dims[2][1]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims[0][1], sim_dims[0][1]],
        [sim_dims[1][0], sim_dims[1][0]],
        [sim_dims[2][0], sim_dims[2][1]],
        c="k",
        linewidth=1.0,
    )
    ax.plot(
        [sim_dims[0][1], sim_dims[0][1]],
        [sim_dims[1][1], sim_dims[1][1]],
        [sim_dims[2][0], sim_dims[2][1]],
        c="k",
        linewidth=1.0,
    )
    ax.set_xlim([sim_dims[0][0], sim_dims[0][1]])
    ax.set_ylim([sim_dims[1][0], sim_dims[1][1]])
    ax.set_zlim([sim_dims[2][0], sim_dims[2][1]])
    # 03 for clusters (material agnostic)
    plt.savefig(
        os.path.join(output_dir, "03_clusters.png"), bbox_inches="tight", dpi=300
    )
    print("3D cluster figure saved as", os.path.join(output_dir, "03_clusters.png"))
    plt.clf()


def determine_molecule_IDs(
    CG_to_AAID_master, AA_morphology_dict, parameter_dict, chromophore_list
):
    print("Determining molecule IDs...")
    CGID_to_mol_ID = {}
    if CG_to_AAID_master is not None:
        # Normal operation with a CGMorphology defined (fine-graining was performed)
        for mol_ID, mol_dict in enumerate(CG_to_AAID_master):
            for CGID in list(mol_dict.keys()):
                CGID_to_mol_ID[CGID] = mol_ID
    elif (len(parameter_dict["CG_site_species"]) == 1) and (
        ("AA_rigid_body_species" not in parameter_dict)
        or (len(parameter_dict["AA_rigid_body_species"]) == 0)
    ):
        print(
            "Small-molecule system detected, assuming each chromophore is its own molecule..."
        )
        # When CGMorphology doesn't exist, and no rigid body species have been specified, then
        # every chromophore is its own molecule)
        for index, chromo in enumerate(chromophore_list):
            for CGID in chromo.CGIDs:
                CGID_to_mol_ID[CGID] = chromo.ID
    else:
        # No CGMorphology, but not small molecules either, so determine molecules based on bonds
        print(
            "Polymeric system detected, determining molecules based on AA bonds (slow calculation)..."
        )
        molecule_AAIDs, molecule_lengths = split_molecules(AA_morphology_dict)
        for index, molecule_AAID_list in enumerate(molecule_AAIDs):
            for AAID in molecule_AAID_list:
                CGID_to_mol_ID[AAID] = index
    # Convert to chromo_to_mol_ID
    chromo_to_mol_ID = {}
    for chromophore in chromophore_list:
        first_CGID = chromophore.CGIDs[0]
        chromo_to_mol_ID[chromophore.ID] = CGID_to_mol_ID[first_CGID]
    return chromo_to_mol_ID


def plot_energy_levels(output_dir, chromophore_list, data_dict):
    HOMO_levels = []
    LUMO_levels = []
    donor_delta_E_ij = []
    acceptor_delta_E_ij = []
    donor_lambda_ij = None
    acceptor_lambda_ij = None
    for chromo in chromophore_list:
        if chromo.species == "donor":
            HOMO_levels.append(chromo.HOMO)
            for neighbour_index, delta_E_ij in enumerate(chromo.neighbours_delta_E):
                if (delta_E_ij is not None) and (
                    chromo.neighbours_TI[neighbour_index] is not None
                ):
                    donor_delta_E_ij.append(delta_E_ij)
                if "reorganisation_energy" in chromo.__dict__:
                    donor_lambda_ij = chromo.reorganisation_energy
        else:
            LUMO_levels.append(chromo.LUMO)
            for neighbour_index, delta_E_ij in enumerate(chromo.neighbours_delta_E):
                if (delta_E_ij is not None) and (
                    chromo.neighbours_TI[neighbour_index] is not None
                ):
                    acceptor_delta_E_ij.append(delta_E_ij)
                if "reorganisation_energy" in chromo.__dict__:
                    acceptor_lambda_ij = chromo.reorganisation_energy
    if len(donor_delta_E_ij) > 0:
        donor_bin_edges, donor_fit_args, donor_mean, donor_std = gauss_fit(
            donor_delta_E_ij
        )
        data_dict["donor_delta_E_ij_mean"] = donor_mean
        data_dict["donor_delta_E_ij_std"] = donor_std
        data_dict["donor_delta_E_ij_err"] = donor_std / np.sqrt(len(donor_delta_E_ij))
        HOMO_av = np.average(HOMO_levels)
        HOMO_std = np.std(HOMO_levels)
        HOMO_err = HOMO_std / np.sqrt(len(HOMO_levels))
        data_dict["donor_frontier_MO_mean"] = HOMO_av
        data_dict["donor_frontier_MO_std"] = HOMO_std
        data_dict["donor_frontier_MO_err"] = HOMO_err
        print("Donor HOMO Level =", HOMO_av, "+/-", HOMO_err)
        print(
            "Donor Delta E_ij stats: mean =",
            donor_mean,
            "+/-",
            donor_std / np.sqrt(len(donor_delta_E_ij)),
        )
        # 06 for donor delta Eij
        plot_delta_E_ij(
            donor_delta_E_ij,
            donor_bin_edges,
            donor_fit_args,
            "donor",
            os.path.join(output_dir, "06_donor_delta_E_ij.png"),
            donor_lambda_ij,
        )
    if len(acceptor_delta_E_ij) > 0:
        acceptor_bin_edges, acceptor_fit_args, acceptor_mean, acceptor_std = gauss_fit(
            acceptor_delta_E_ij
        )
        data_dict["acceptor_delta_E_ij_mean"] = acceptor_mean
        data_dict["acceptor_delta_E_ij_std"] = acceptor_std
        data_dict["acceptor_delta_E_ij_err"] = acceptor_std / np.sqrt(
            len(acceptor_delta_E_ij)
        )
        LUMO_av = np.average(LUMO_levels)
        LUMO_std = np.std(LUMO_levels)
        LUMO_err = LUMO_std / np.sqrt(len(LUMO_levels))
        data_dict["acceptor_frontier_MO_mean"] = LUMO_av
        data_dict["acceptor_frontier_MO_std"] = LUMO_std
        data_dict["acceptor_frontier_MO_err"] = LUMO_err
        print("Acceptor LUMO Level =", LUMO_av, "+/-", LUMO_err)
        print(
            "Acceptor Delta E_ij stats: mean =",
            acceptor_mean,
            "+/-",
            acceptor_std / np.sqrt(len(acceptor_delta_E_ij)),
        )
        # 07 for acceptor delta Eij
        plot_delta_E_ij(
            acceptor_delta_E_ij,
            acceptor_bin_edges,
            acceptor_fit_args,
            "acceptor",
            os.path.join(output_dir, "07_acceptor_delta_E_ij.png"),
            acceptor_lambda_ij,
        )
    return data_dict


def plot_delta_E_ij(delta_E_ij, gauss_bins, fit_args, data_type, file_name, lambda_ij):
    plt.figure()
    n, bins, patches = plt.hist(
        delta_E_ij, np.arange(np.min(delta_E_ij), np.max(delta_E_ij), 0.05), color=["b"]
    )
    if fit_args is not None:
        gauss_Y = gaussian(gauss_bins[:-1], *fit_args)
        scale_factor = max(n) / max(gauss_Y)
        plt.plot(gauss_bins[:-1], gauss_Y * scale_factor, "ro:")
    else:
        print("No Gaussian found (probably zero-width delta function)")
    if lambda_ij is not None:
        plt.axvline(-float(lambda_ij), c="k")
    plt.ylabel("Frequency (Arb. U.)")
    plt.xlabel("".join([data_type.capitalize(), r" $\Delta E_{i,j}$ (eV)"]))
    # plt.xlim([-0.5, 0.5])
    plt.savefig(file_name, dpi=300)
    plt.close()
    print("Figure saved as", file_name)


def plot_mixed_hopping_rates(
    output_dir,
    chromophore_list,
    parameter_dict,
    cluster_dicts,
    chromo_to_mol_ID,
    data_dict,
    AA_morphology_dict,
    cut_off_dict,
):
    # Create all the empty lists we need
    hop_types = ["intra", "inter"]
    hop_targets = ["cluster", "mol"]
    hop_properties = ["rates", "TIs"]
    chromo_species = ["donor", "acceptor"]
    property_lists = {}
    for property_name in [
        "_".join([hop_type, hop_target, hop_property, species])
        for hop_type in hop_types
        for hop_target in hop_targets
        for hop_property in hop_properties
        for species in chromo_species
    ]:
        property_lists[property_name] = []
    T = 290
    for chromo in chromophore_list:
        mol1ID = chromo_to_mol_ID[chromo.ID]
        for index, T_ij in enumerate(chromo.neighbours_TI):
            if (T_ij is None) or (T_ij == 0):
                continue
            chromo2 = chromophore_list[chromo.neighbours[index][0]]
            mol2ID = chromo_to_mol_ID[chromo2.ID]
            delta_E = chromo.neighbours_delta_E[index]
            if chromo.sub_species == chromo2.sub_species:
                lambda_ij = chromo.reorganisation_energy
            else:
                lambda_ij = (
                    chromo.reorganisation_energy + chromo2.reorganisation_energy
                ) / 2
            # Now take into account the various behaviours we can have from the parameter file
            prefactor = 1.0
            # Apply the koopmans prefactor
            try:
                use_koop = parameter_dict["use_koopmans_approximation"]
                if use_koop:
                    prefactor *= parameter_dict["koopmans_hopping_prefactor"]
            except KeyError:
                pass
            # Apply the simple energetic penalty model
            try:
                boltz_pen = parameter_dict["use_simple_energetic_penalty"]
            except KeyError:
                boltz_pen = False
            # Apply the distance penalty due to VRH
            try:
                VRH = parameter_dict["use_VRH"]
                if VRH is True:
                    VRH_delocalisation = 1.0 / chromo.VRH_delocalisation
            except KeyError:
                VRH = False
            if VRH is True:
                relative_image = chromo.neighbours[index][1]
                neighbour_chromo_posn = chromo2.posn + (
                    np.array(relative_image)
                    * np.array(
                        [AA_morphology_dict[axis] for axis in ["lx", "ly", "lz"]]
                    )
                )
                chromophore_separation = (
                    hf.calculate_separation(chromo.posn, neighbour_chromo_posn) * 1E-10
                )
                rate = hf.calculate_carrier_hop_rate(
                    lambda_ij * elementary_charge,
                    T_ij * elementary_charge,
                    delta_E * elementary_charge,
                    prefactor,
                    T,
                    use_VRH=VRH,
                    rij=chromophore_separation,
                    VRH_delocalisation=VRH_delocalisation,
                    boltz_pen=boltz_pen,
                )
            else:
                rate = hf.calculate_carrier_hop_rate(
                    lambda_ij * elementary_charge,
                    T_ij * elementary_charge,
                    delta_E * elementary_charge,
                    prefactor,
                    T,
                    boltz_pen=boltz_pen,
                )
            if chromo2.ID < chromo.ID:
                continue
            # Do intra- / inter- clusters
            if chromo.species == "acceptor":
                if (
                    cluster_dicts[1][chromo.ID]
                    == cluster_dicts[1][chromo.neighbours[index][0]]
                ):
                    property_lists["intra_cluster_rates_acceptor"].append(rate)
                    property_lists["intra_cluster_TIs_acceptor"].append(T_ij)
                else:
                    property_lists["inter_cluster_rates_acceptor"].append(rate)
                    property_lists["inter_cluster_TIs_acceptor"].append(T_ij)
            else:
                if (
                    cluster_dicts[0][chromo.ID]
                    == cluster_dicts[0][chromo.neighbours[index][0]]
                ):
                    property_lists["intra_cluster_rates_donor"].append(rate)
                    property_lists["intra_cluster_TIs_donor"].append(T_ij)
                else:
                    property_lists["inter_cluster_rates_donor"].append(rate)
                    property_lists["inter_cluster_TIs_donor"].append(T_ij)
            # Now do intra- / inter- molecules
            if mol1ID == mol2ID:
                if chromo.species == "acceptor":
                    property_lists["intra_mol_rates_acceptor"].append(rate)
                    property_lists["intra_mol_TIs_acceptor"].append(T_ij)
                else:
                    property_lists["intra_mol_rates_donor"].append(rate)
                    property_lists["intra_mol_TIs_donor"].append(T_ij)
            else:
                if chromo.species == "acceptor":
                    property_lists["inter_mol_rates_acceptor"].append(rate)
                    property_lists["inter_mol_TIs_acceptor"].append(T_ij)
                else:
                    property_lists["inter_mol_rates_donor"].append(rate)
                    property_lists["inter_mol_TIs_donor"].append(T_ij)
    # 12 for the donor cluster TI, 13 for the acceptor cluster TI, 14 for
    # the donor mol kij, 15 for the acceptor mol kij, 16 for the donor
    # cluster kij, 17 for the acceptor cluster kij
    # Donor cluster Plots:
    if (len(property_lists["intra_cluster_rates_donor"]) > 0) or (
        len(property_lists["inter_cluster_rates_donor"]) > 0
    ):
        print(
            "Mean intra-cluster donor rate =",
            np.mean(property_lists["intra_cluster_rates_donor"]),
            "+/-",
            np.std(property_lists["intra_cluster_rates_donor"])
            / float(len(property_lists["intra_cluster_rates_donor"])),
        )
        print(
            "Mean inter-cluster donor rate =",
            np.mean(property_lists["inter_cluster_rates_donor"]),
            "+/-",
            np.std(property_lists["inter_cluster_rates_donor"])
            / float(len(property_lists["inter_cluster_rates_donor"])),
        )
        plot_stacked_hist_rates(
            property_lists["intra_cluster_rates_donor"],
            property_lists["inter_cluster_rates_donor"],
            ["Intra-cluster", "Inter-cluster"],
            "donor",
            os.path.join(output_dir, "16_donor_hopping_rate_clusters.png"),
        )
        plot_stacked_hist_TIs(
            property_lists["intra_cluster_TIs_donor"],
            property_lists["inter_cluster_TIs_donor"],
            ["Intra-cluster", "Inter-cluster"],
            "donor",
            os.path.join(output_dir, "12_donor_transfer_integral_clusters.png"),
            cut_off_dict["TI"][0],
        )
    # Acceptor cluster Plots:
    if (len(property_lists["intra_cluster_rates_acceptor"]) > 0) or (
        len(property_lists["inter_cluster_rates_acceptor"]) > 0
    ):
        print(
            "Mean intra-cluster acceptor rate =",
            np.mean(property_lists["intra_cluster_rates_acceptor"]),
            "+/-",
            np.std(property_lists["intra_cluster_rates_acceptor"])
            / float(len(property_lists["intra_cluster_rates_acceptor"])),
        )
        print(
            "Mean inter-cluster acceptor rate =",
            np.mean(property_lists["inter_cluster_rates_acceptor"]),
            "+/-",
            np.std(property_lists["inter_cluster_rates_acceptor"])
            / float(len(property_lists["inter_cluster_rates_acceptor"])),
        )
        plot_stacked_hist_rates(
            property_lists["intra_cluster_rates_acceptor"],
            property_lists["inter_cluster_rates_acceptor"],
            ["Intra-cluster", "Inter-cluster"],
            "acceptor",
            os.path.join(output_dir, "17_acceptor_hopping_rate_clusters.png"),
        )
        plot_stacked_hist_TIs(
            property_lists["intra_cluster_TIs_acceptor"],
            property_lists["inter_cluster_TIs_acceptor"],
            ["Intra-cluster", "Inter-cluster"],
            "acceptor",
            os.path.join(output_dir, "13_acceptor_transfer_integral_clusters.png"),
            cut_off_dict["TI"][1],
        )
    # Donor Mol Plots:
    if (len(property_lists["intra_mol_rates_donor"]) > 0) or (
        len(property_lists["inter_mol_rates_donor"]) > 0
    ):
        print(
            "Mean intra-molecular donor rate =",
            np.mean(property_lists["intra_mol_rates_donor"]),
            "+/-",
            np.std(property_lists["intra_mol_rates_donor"])
            / float(len(property_lists["intra_mol_rates_donor"])),
        )
        print(
            "Mean inter-molecular donor rate =",
            np.mean(property_lists["inter_mol_rates_donor"]),
            "+/-",
            np.std(property_lists["inter_mol_rates_donor"])
            / float(len(property_lists["inter_mol_rates_donor"])),
        )
        plot_stacked_hist_rates(
            property_lists["intra_mol_rates_donor"],
            property_lists["inter_mol_rates_donor"],
            ["Intra-mol", "Inter-mol"],
            "donor",
            os.path.join(output_dir, "14_donor_hopping_rate_mols.png"),
        )
    # Acceptor Mol Plots:
    if (len(property_lists["intra_mol_rates_acceptor"]) > 0) or (
        len(property_lists["inter_mol_rates_acceptor"]) > 0
    ):
        print(
            "Mean intra-molecular acceptor rate =",
            np.mean(property_lists["intra_mol_rates_acceptor"]),
            "+/-",
            np.std(property_lists["intra_mol_rates_acceptor"])
            / float(len(property_lists["intra_mol_rates_acceptor"])),
        )
        print(
            "Mean inter-molecular acceptor rate =",
            np.mean(property_lists["inter_mol_rates_acceptor"]),
            "+/-",
            np.std(property_lists["inter_mol_rates_acceptor"])
            / float(len(property_lists["inter_mol_rates_acceptor"])),
        )
        plot_stacked_hist_rates(
            property_lists["intra_mol_rates_acceptor"],
            property_lists["inter_mol_rates_acceptor"],
            ["Intra-mol", "Inter-mol"],
            "acceptor",
            os.path.join(output_dir, "15_acceptor_hopping_rate_mols.png"),
        )
    # Update the dataDict
    for material in chromo_species:
        for hop_type in hop_types:
            for hop_target in hop_targets:
                number_of_hops = len(
                    property_lists[
                        "".join([hop_type, "_", hop_target, "_rates_", material])
                    ]
                )
                if number_of_hops == 0:
                    continue
                other_hop_type = hop_types[int((hop_types.index(hop_type) * -1) + 1)]
                proportion = number_of_hops / (
                    number_of_hops
                    + len(
                        property_lists[
                            "".join(
                                [other_hop_type, "_", hop_target, "_rates_", material]
                            )
                        ]
                    )
                )
                mean_rate = np.mean(
                    property_lists[
                        "".join([hop_type, "_", hop_target, "_rates_", material])
                    ]
                )
                deviation_rate = np.std(
                    property_lists[
                        "".join([hop_type, "_", hop_target, "_rates_", material])
                    ]
                )
                data_dict[
                    "".join(
                        [
                            material.lower(),
                            "_",
                            hop_type,
                            "_",
                            hop_target.lower(),
                            "_hops",
                        ]
                    )
                ] = number_of_hops
                data_dict[
                    "".join(
                        [
                            material.lower(),
                            "_",
                            hop_type,
                            "_",
                            hop_target.lower(),
                            "_proportion",
                        ]
                    )
                ] = proportion
                data_dict[
                    "".join(
                        [
                            material.lower(),
                            "_",
                            hop_type,
                            "_",
                            hop_target.lower(),
                            "_rate_mean",
                        ]
                    )
                ] = mean_rate
                data_dict[
                    "".join(
                        [
                            material.lower(),
                            "_",
                            hop_type,
                            "_",
                            hop_target.lower(),
                            "_rate_std",
                        ]
                    )
                ] = deviation_rate
    return data_dict


def plot_stacked_hist_rates(data1, data2, labels, data_type, file_name):
    plt.figure()
    (n, bins, patches) = plt.hist(
        [data1, data2],
        bins=np.logspace(1, 18, 40),
        stacked=True,
        color=["r", "b"],
        label=labels,
    )
    plt.ylabel("Frequency (Arb. U.)")
    plt.xlabel("".join([data_type.capitalize(), r" k$_{i,j}$ (s$^{-1}$)"]))
    plt.xlim([1, 1E18])
    plt.xticks([1E0, 1E3, 1E6, 1E9, 1E12, 1E15, 1E18])
    # plt.ylim([0, 14000])
    plt.ylim([0, np.max(n) * 1.02])
    plt.legend(loc=2, prop={"size": 18})
    plt.gca().set_xscale("log")
    plt.savefig(file_name, dpi=300)
    plt.close()
    print("Figure saved as", file_name)


def plot_stacked_hist_TIs(data1, data2, labels, data_type, file_name, cut_off):
    plt.figure()
    (n, bins, patches) = plt.hist(
        [data1, data2],
        bins=np.linspace(0, 1.2, 20),
        stacked=True,
        color=["r", "b"],
        label=labels,
    )
    plt.ylabel("Frequency (Arb. U.)")
    plt.xlabel("".join([data_type.capitalize(), r" J$_{i,j}$ (eV)"]))
    # plt.xlim([0, 1.2])
    plt.ylim([0, np.max(n) * 1.02])
    if cut_off is not None:
        plt.axvline(float(cut_off), c="k")
    plt.legend(loc=0, prop={"size": 18})
    plt.savefig(file_name, dpi=300)
    plt.close()
    print("Figure saved as", file_name)


def write_CSV(data_dict, directory):
    CSV_file_name = os.path.join(directory, "results.csv")
    with open(CSV_file_name, "w+") as CSV_file:
        CSV_writer = csv.writer(CSV_file)
        for key in sorted(data_dict.keys()):
            CSV_writer.writerow([key, data_dict[key]])
    print("CSV file written to {:s}".format(CSV_file_name))


def create_results_pickle(directory):
    cores_list = []
    for file_name in glob.glob(os.path.join(directory, "KMC", "*")):
        if "log" not in file_name:
            continue
        try:
            cores_list.append(os.path.split(file_name)[1].split(".")[0].split("_")[-1])
        except IndexError:
            pass
    cores_list = sorted(list(set(cores_list)))
    results_pickles_list = []
    keep_list = []
    for core in cores_list:
        # Check if there is already a finished KMC_results pickle
        main = os.path.join(
            directory, "KMC", "KMC_results_{:02d}.pickle".format(int(core))
        )
        if os.path.exists(main):
            results_pickles_list.append(main)
            keep_list.append(None)
            continue
        # If not, find the slot1 and slot2 pickle that is most recent
        slot1 = os.path.join(
            directory, "KMC", "KMC_slot1_results_{:02d}.pickle".format(int(core))
        )
        slot2 = os.path.join(
            directory, "KMC", "KMC_slot2_results_{:02d}.pickle".format(int(core))
        )
        if os.path.exists(slot1) and not os.path.exists(slot2):
            keep_list.append(slot1)
        elif os.path.exists(slot2) and not os.path.exists(slot1):
            keep_list.append(slot2)
        elif os.path.getsize(slot1) >= os.path.getsize(slot2):
            keep_list.append(slot1)
        else:
            keep_list.append(slot2)
    print("{:d} pickle files found to combine!".format(len(keep_list)))
    print("Combining", keep_list)
    for keeper in zip(cores_list, keep_list):
        # Skip this core if we already have a finished KMC_results for it
        if keeper[1] is None:
            continue
        new_name = os.path.join(
            directory, "KMC", "KMC_results_{}.pickle".format(keeper[0])
        )
        shutil.copyfile(str(keeper[1]), new_name)
        results_pickles_list.append(new_name)
    combine_results_pickles(directory, results_pickles_list)


def combine_results_pickles(directory, pickle_files):
    combined_data = {}
    pickle_files = sorted(pickle_files)
    for file_name in pickle_files:
        # The pickle was repeatedly dumped to, in order to save time.
        # Each dump stream is self-contained, so iteratively unpickle to add the new data.
        with open(file_name, "rb") as pickle_file:
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
    combined_file_loc = os.path.join(directory, "KMC", "KMC_results.pickle")
    with open(combined_file_loc, "wb+") as pickle_file:
        pickle.dump(combined_data, pickle_file)
    print("Complete data written to", combined_file_loc)


def calculate_cut_off_from_dist(
    bin_centres,
    frequencies,
    minimum_index=None,
    maximum_index=None,
    value_at_least=100,
    logarithmic=False,
):
    try:
        if minimum_index is not None:
            # Looking for minima
            minima = argrelextrema(frequencies, np.less)[0]
            if minimum_index < 0:
                # Sometimes a tiny minimum at super RHS breaks this
                cut_off = 0.0
                while True:
                    selected_minimum = minima[minimum_index]
                    cut_off = bin_centres[selected_minimum]
                    if frequencies[selected_minimum] > value_at_least:
                        break
                    minimum_index -= 1
            else:
                # Sometimes a tiny maximum at super LHS breaks this
                cut_off = 0.0
                while True:
                    selected_minimum = minima[minimum_index]
                    cut_off = bin_centres[selected_minimum]
                    if frequencies[selected_minimum] > value_at_least:
                        break
                    minimum_index += 1
        elif maximum_index is not None:
            # Looking for maxima
            maxima = argrelextrema(frequencies, np.greater)[0]
            if maximum_index < 0:
                # Sometimes a tiny maximum at super RHS breaks this
                cut_off = 0.0
                while True:
                    selected_maximum = maxima[maximum_index]
                    cut_off = bin_centres[selected_maximum]
                    if frequencies[selected_maximum] > value_at_least:
                        break
                    maximum_index -= 1
            else:
                # Sometimes a tiny maximum at super LHS breaks this
                cut_off = 0.0
                while True:
                    selected_maximum = maxima[maximum_index]
                    cut_off = bin_centres[selected_maximum]
                    if frequencies[selected_maximum] > value_at_least:
                        break
                    maximum_index += 1
        if logarithmic is True:
            cut_off = 10 ** cut_off
        # Return as string, as it will be converted to a float later
        return str(cut_off)
    except IndexError:
        print(
            "EXCEPTION: No minima found in frequency distribution. Setting cut_off to None."
        )
        return None


def plot_TI_hist(
    chromophore_list, chromo_to_mol_ID, output_dir, TI_cut_donor, TI_cut_acceptor
):
    # TI_dist [[DONOR], [ACCEPTOR]]
    TI_dist_intra = [[], []]
    TI_dist_inter = [[], []]
    material = ["donor", "acceptor"]
    labels = ["Intra-mol", "Inter-mol"]
    TI_cuts = [TI_cut_donor, TI_cut_acceptor]
    for material_index, material_type in enumerate(material):
        for chromo1 in chromophore_list:
            for neighbour_index, chromo2_details in enumerate(chromo1.neighbours):
                chromo2 = chromophore_list[chromo2_details[0]]
                if (chromo2_details is None) or (chromo1.ID >= chromo2.ID):
                    continue
                if chromo_to_mol_ID[chromo1.ID] == chromo_to_mol_ID[chromo2.ID]:
                    TI_dist_intra[material_index].append(
                        chromo1.neighbours_TI[neighbour_index]
                    )
                else:
                    TI_dist_inter[material_index].append(
                        chromo1.neighbours_TI[neighbour_index]
                    )
        if (len(TI_dist_intra[material_index]) == 0) and (
            len(TI_dist_inter[material_index])
        ):
            continue
        plt.figure()
        n, bin_edges, _ = plt.hist(
            [TI_dist_intra[material_index], TI_dist_inter[material_index]],
            bins=np.linspace(
                0,
                np.max(TI_dist_intra[material_index] + TI_dist_inter[material_index]),
                20,
            ),
            color=["r", "b"],
            stacked=True,
            label=labels,
        )
        bin_centres = (bin_edges[1:] + bin_edges[:-1]) / 2.0
        smoothed_n = gaussian_filter(n[0] + n[1], 1.0)
        plt.plot(bin_centres, smoothed_n, color="r")
        if (TI_cuts[material_index] is not None) and (
            TI_cuts[material_index].lower() == "auto"
        ):
            TI_cuts[material_index] = calculate_cut_off_from_dist(
                bin_centres, smoothed_n, minimum_index=-1, value_at_least=100
            )
        if TI_cuts[material_index] is not None:
            print(
                "Cluster cut-off based on",
                material[material_index],
                "transfer integrals set to",
                TI_cuts[material_index],
            )
            plt.axvline(float(TI_cuts[material_index]), c="k")
        plt.xlim(
            [0, np.max(TI_dist_intra[material_index] + TI_dist_inter[material_index])]
        )
        plt.ylim([0, np.max(n) * 1.02])
        plt.ylabel("Frequency (Arb. U.)")
        plt.xlabel("".join([material[material_index].capitalize(), r" J$_{i,j}$ (eV)"]))
        plt.legend(loc=1, prop={"size": 18})
        # 10 for donor TI mols dist, 11 for acceptor TI mols dist,
        file_name = "".join(
            [
                "{:02}_".format(10 + material_index),
                material_type,
                "_transfer_integral_mols",
                ".png",
            ]
        )
        plt.savefig(os.path.join(output_dir, file_name), dpi=300)
        plt.clf()
        print("Figure saved as", os.path.join(output_dir, file_name))
    return TI_cuts[0], TI_cuts[1]


def plot_frequency_dist(directory, carrier_type, carrier_history, cut_off):
    carrier_types = ["hole", "electron"]
    non_zero_indices = carrier_history.nonzero()
    coordinates = list(zip(non_zero_indices[0], non_zero_indices[1]))
    frequencies = []
    for coords in coordinates:
        if coords[1] < coords[0]:
            # Only consider hops in one direction
            continue
        frequency = carrier_history[coords]
        frequencies.append(np.log10(frequency))
    plt.figure()
    n, bin_edges, _ = plt.hist(frequencies, bins=60, color="b")
    bin_centres = (bin_edges[1:] + bin_edges[:-1]) / 2.0
    smoothed_n = gaussian_filter(n, 1.0)
    plt.plot(bin_centres, smoothed_n, color="r")
    if (cut_off is not None) and (cut_off.lower() == "auto"):
        print("DYNAMIC CUT")
        cut_off = calculate_cut_off_from_dist(
            bin_centres,
            smoothed_n,
            minimum_index=-1,
            value_at_least=100,
            logarithmic=True,
        )
    if cut_off is not None:
        print("Cluster cut-off based on hop frequency set to", cut_off)
        plt.axvline(np.log10(float(cut_off)), c="k")
    plt.xlabel("".join(["Total ", carrier_type, " hops (Arb. U.)"]))
    ax = plt.gca()
    # tick_labels = np.arange(0, 7, 1)
    # plt.xlim([0, 6])
    tick_labels = np.arange(0, np.ceil(np.max(frequencies)) + 1, 1)
    plt.xlim([0, np.ceil(np.max(frequencies))])
    plt.xticks(tick_labels, [r"10$^{{{}}}$".format(int(x)) for x in tick_labels])
    plt.ylabel("Frequency (Arb. U.)")
    # 24 for hole hop frequency dist, 25 for electron hop frequency dist
    file_name = "".join(
        [
            "{:02}_total_hop_freq_".format(24 + carrier_types.index(carrier_type)),
            carrier_type,
            ".png",
        ]
    )
    plt.savefig(os.path.join(directory, "figures", file_name), dpi=300)
    plt.clf()
    print("Figure saved as", os.path.join(directory, "figures", file_name))
    return cut_off


def plot_net_frequency_dist(directory, carrier_type, carrier_history):
    carrier_types = ["hole", "electron"]
    non_zero_indices = carrier_history.nonzero()
    coordinates = list(zip(non_zero_indices[0], non_zero_indices[1]))
    frequencies = []
    for coords in coordinates:
        if coords[1] < coords[0]:
            # Only consider hops in one direction
            continue
        frequency = np.abs(carrier_history[coords] - carrier_history[coords[::-1]])
        if frequency > 0:
            frequencies.append(np.log10(frequency))
    plt.figure()
    plt.hist(frequencies, bins=60, color="b")
    plt.xlabel("".join(["Net ", carrier_type, " hops (Arb. U.)"]))
    ax = plt.gca()
    tick_labels = np.arange(0, np.ceil(np.max(frequencies)) + 1, 1)
    plt.xlim([0, np.ceil(np.max(frequencies))])
    plt.xticks(tick_labels, [r"10$^{{{}}}$".format(int(x)) for x in tick_labels])
    # plt.ylim([0, 2500])
    plt.ylabel("Frequency (Arb. U.)")
    # 26 for hole hop frequency dist, 27 for electron hop frequency dist
    file_name = "".join(
        [
            "{:02}_net_hop_freq_".format(26 + carrier_types.index(carrier_type)),
            carrier_type,
            ".png",
        ]
    )
    plt.savefig(os.path.join(directory, "figures", file_name), dpi=300)
    plt.clf()
    print("Figure saved as", os.path.join(directory, "figures", file_name))


def plot_discrepancy_frequency_dist(directory, carrier_type, carrier_history):
    carrier_types = ["hole", "electron"]
    non_zero_indices = carrier_history.nonzero()
    coordinates = list(zip(non_zero_indices[0], non_zero_indices[1]))
    frequencies = []
    net_equals_total = 0
    net_near_total = 0
    for coords in coordinates:
        if coords[1] < coords[0]:
            # Only consider hops in one direction
            continue
        total_hops = carrier_history[coords] + carrier_history[coords[::-1]]
        net_hops = np.abs(carrier_history[coords] - carrier_history[coords[::-1]])
        frequency = total_hops - net_hops
        if frequency == 0:
            net_equals_total += 1
            net_near_total += 1
        else:
            if frequency < 10:
                net_near_total += 1
            frequencies.append(np.log10(frequency))
    plt.figure()
    plt.hist(frequencies, bins=60, color="b")
    plt.xlabel("".join(["Discrepancy (Arb. U.)"]))
    ax = plt.gca()
    tick_labels = np.arange(0, np.ceil(np.max(frequencies)) + 1, 1)
    plt.xlim([0, np.ceil(np.max(frequencies))])
    plt.xticks(tick_labels, [r"10$^{{{}}}$".format(int(x)) for x in tick_labels])
    # plt.ylim([0, 13000])
    plt.ylabel("Frequency (Arb. U.)")
    # 28 for hole hop frequency dist, 29 for electron hop frequency dist
    file_name = "".join(
        [
            "{:02}_hop_discrepancy_".format(28 + carrier_types.index(carrier_type)),
            carrier_type,
            ".png",
        ]
    )
    plt.savefig(os.path.join(directory, "figures", file_name), dpi=300)
    plt.clf()
    print(
        "There are",
        net_equals_total,
        "paths in this morphology with one-way transport.",
    )
    print(
        "There are", net_near_total, "paths in this morphology with total - net < 10."
    )
    print("Figure saved as", os.path.join(directory, "figures", file_name))


def calculate_mobility(
    directory,
    current_carrier_type,
    times,
    MSDs,
    time_standard_errors,
    MSD_standard_errors,
):
    # Create the first figure that will be replotted each time
    plt.figure()
    times, MSDs = hf.parallel_sort(times, MSDs)
    mobility, mob_error, r_squared = plot_MSD(
        times,
        MSDs,
        time_standard_errors,
        MSD_standard_errors,
        directory,
        current_carrier_type,
    )
    print("----------====================----------")
    print(
        current_carrier_type.capitalize(),
        "mobility for",
        directory,
        "= {0:.2E} +- {1:.2E}".format(mobility, mob_error),
        "cm^{2} V^{-1} s^{-1}",
    )
    print("----------====================----------")
    return mobility, mob_error, r_squared


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-t",
        "--three_D",
        action="store_true",
        required=False,
        help=(
            "If present, use matplotlib to plot the 3D graphs (3D network, anisotropy and cluster"
            " positions. This takes a while (usually a couple of minutes) to plot."
            " Default = False"
        ),
    )
    parser.add_argument(
        "-sd",
        "--sequence_donor",
        type=lambda s: [float(item) for item in s.split(",")],
        default=None,
        required=False,
        help=(
            "Create a figure in the current directory that describes the evolution"
            " of the hole anisotropy/mobility using the specified comma-delimited"
            " string as the sequence of x values. For instance -s"
            " '1.5,1.75,2.0,2.25,2.5' will assign each of the 5 following directories"
            " these x-values when plotting the hole mobility evolution."
        ),
    )
    parser.add_argument(
        "-sa",
        "--sequence_acceptor",
        type=lambda s: [float(item) for item in s.split(",")],
        default=None,
        required=False,
        help=(
            "Create a figure in the current directory that describes the evolution"
            " of the electron anisotropy/mobility using the specified comma-delimited"
            " string as the sequence of x values. For instance -s"
            " '1.5,1.75,2.0,2.25,2.5' will assign each of the 5 following directories"
            " these x-values when plotting the electron mobility evolution."
        ),
    )
    parser.add_argument(
        "-cod",
        "--o_cut_donor",
        default=None,
        required=False,
        type=str,
        help=(
            "Specify the orientation cut-off (in degrees) for the donor material for determining"
            " which chromophores belong to the same cluster. Chromophores with angle between normal"
            " vectors > o_cut_donor will be considered as different crystals. Default = None (Note: if all"
            " cut-offs are specified as None, then the entire morphology will be considered"
            " as a single crystal)."
        ),
    )
    parser.add_argument(
        "-coa",
        "--o_cut_acceptor",
        default=None,
        required=False,
        type=str,
        help=(
            "Specify the orientation cut-off (in degrees) for the acceptor material for determining"
            " which chromophores belong to the same cluster. Chromophores with angle between normal"
            " vectors > o_cut_acceptor will be considered as different crystals. Default = None (Note: if all"
            " cut-offs are specified as None, then the entire morphology will be considered"
            " as a single crystal)."
        ),
    )
    parser.add_argument(
        "-ctd",
        "--ti_cut_donor",
        default=None,
        required=False,
        type=str,
        help=(
            "Specify the transfer integral cut-off (in eV) for the donor material for determining"
            " which chromophores belong to the same cluster. Chromophores with hopping transfer"
            " integral < ti_cut_donor will be considered as different crystals. Default = None (Note: if all"
            " cut-offs are specified as None, then the entire morphology will be considered"
            " as a single crystal)."
        ),
    )
    parser.add_argument(
        "-cta",
        "--ti_cut_acceptor",
        default=None,
        required=False,
        type=str,
        help=(
            "Specify the transfer integral cut-off (in eV) for the acceptor material for determining"
            " which chromophores belong to the same cluster. Chromophores with hopping transfer"
            " integral < ti_cut_acceptor will be considered as different crystals. Default = None (Note: if all"
            " cut-offs are specified as None, then the entire morphology will be considered"
            " as a single crystal."
        ),
    )
    parser.add_argument(
        "-cfd",
        "--freq_cut_donor",
        default=None,
        required=False,
        type=str,
        help=(
            "Specify the hopping frequency cut-off for the donor material for determining"
            " which chromophores belong to the same cluster. Chromophores connected by fewer than"
            " freq_cut_donor total hops in the KMC simulation will be considered as different crystals."
            " Use 'auto' to have KMCAnalyse set the cut-off to the final minimum of the frequency distribution."
            " Default = None (Note: if all cut-offs are specified as None, then the entire morphology will"
            " be considered as a single crystal."
        ),
    )
    parser.add_argument(
        "-cfa",
        "--freq_cut_acceptor",
        default=None,
        required=False,
        type=str,
        help=(
            "Specify the hopping frequency cut-off for the acceptor material for determining"
            " which chromophores belong to the same cluster. Chromophores connected by fewer than"
            " freq_cut_acceptor total hops in the KMC simulation will be considered as different crystals."
            " Use 'auto' to have KMCAnalyse set the cut-off to the final minimum of the frequency distribution."
            " Default = None (Note: if all cut-offs are specified as None, then the entire morphology will"
            " be considered as a single crystal."
        ),
    )
    parser.add_argument(
        "-crd",
        "--sep_cut_donor",
        default=None,
        required=False,
        type=str,
        help=(
            "Specify the separation cut-off (in Angstroms) for the donor material for determining"
            " which chromophores belong to the same cluster. Chromophores with separation"
            " r > sep_cut_donor will be considered as different crystals. Default = None (Note: if all"
            " cut-offs are specified as None, then the entire morphology will be considered"
            " as a single crystal)."
        ),
    )
    parser.add_argument(
        "-cra",
        "--sep_cut_acceptor",
        default=None,
        required=False,
        type=str,
        help=(
            "Specify the separation cut-off (in Angstroms) for the acceptor material for determining"
            " which chromophores belong to the same cluster. Chromophores with separation"
            " r > sep_cut_acceptor will be considered as different crystals. Default = None (Note: if all"
            " cut-offs are specified as None, then the entire morphology will be considered"
            " as a single crystal)."
        ),
    )
    parser.add_argument(
        "-x",
        "--xlabel",
        default="Temperature (Arb. U.)",
        required=False,
        help=(
            "Specify an x-label for the combined plot (only used if -s is specified). Default ="
            ' "Temperature (Arb. U.)"'
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
        "-tcl",
        "--generate_tcl",
        action="store_true",
        required=False,
        help=(
            "Specify to create a tcl file to use with VMD in which each resid is colored"
            'based off of its cluster. Only functions when used in conjunction with the "-t"'
            "flag. Default = False"
        ),
    )
    args, directory_list = parser.parse_known_args()

    # Load the matplotlib backend and the plotting subroutines
    global plt
    global p3
    if args.backend is not None:
        import matplotlib

        matplotlib.use(args.backend.strip())
    import matplotlib.pyplot as plt

    try:
        import mpl_toolkits.mplot3d as p3
    except ImportError:
        print(
            "Could not import 3D plotting engine, calling the plotMolecule3D function will result in an error!"
        )
    hole_mobility_data = []
    hole_anisotropy_data = []
    electron_mobility_data = []
    electron_anisotropy_data = []
    for directory in directory_list:
        # Create the figures directory if it doesn't already exist
        os.makedirs(os.path.join(directory, "figures"), exist_ok=True)
        # Load in all the required data
        data_dict = {}
        print("\n")
        print("-----===== KMC_ANALYSE =====-----")
        print(directory)
        print("-----=======================-----")
        print("Getting carrier data...")
        carrier_data = load_KMC_results_pickle(directory)
        print("Carrier Data obtained")
        # Now need to split up the carrierData into both electrons and holes
        carrier_data_holes, carrier_data_electrons = split_carriers_by_type(
            carrier_data
        )
        print("Loading chromophore_list...")
        pickle_location = os.path.join(
            directory, "code", "".join([os.path.split(directory)[1], ".pickle"])
        )
        data = hf.load_pickle(pickle_location)
        AA_morphology_dict = data[0]
        if "unwrapped_position" not in AA_morphology_dict.keys():
            AA_morphology_dict = hf.add_unwrapped_positions(AA_morphology_dict)
        CG_morphology_dict = data[1]
        CG_to_AAID_master = data[2]
        parameter_dict = data[3]
        chromophore_list = data[4]
        print("Chromophore_list obtained")
        morphology_shape = np.array(
            [AA_morphology_dict[axis] for axis in ["lx", "ly", "lz"]]
        )
        sim_dims = [
            [-AA_morphology_dict[axis] / 2.0, AA_morphology_dict[axis] / 2.0]
            for axis in ["lx", "ly", "lz"]
        ]

        # Calculate the mobilities
        complete_carrier_types = []
        complete_carrier_data = []
        if (carrier_data_holes is not None) and (len(carrier_data_holes["ID"]) > 0):
            complete_carrier_types.append("hole")
            complete_carrier_data.append(carrier_data_holes)
        if (carrier_data_electrons is not None) and (
            len(carrier_data_electrons["ID"]) > 0
        ):
            complete_carrier_types.append("electron")
            complete_carrier_data.append(carrier_data_electrons)
        carrier_history_dict = {}
        for carrier_type_index, carrier_data in enumerate(complete_carrier_data):
            current_carrier_type = complete_carrier_types[carrier_type_index]
            print("Considering the transport of {:s}...".format(current_carrier_type))
            print("Obtaining mean squared displacements...")
            carrier_history, times, MSDs, time_standard_errors, MSD_standard_errors = get_carrier_data(
                carrier_data
            )
            carrier_history_dict[current_carrier_type] = carrier_history
            print("Plotting distribution of carrier displacements")
            plot_displacement_dist(carrier_data, directory, current_carrier_type)
            print("Calculating mobility...")
            mobility, mob_error, r_squared = calculate_mobility(
                directory,
                current_carrier_type,
                times,
                MSDs,
                time_standard_errors,
                MSD_standard_errors,
            )
            print("Plotting hop vector distribution")
            plot_hop_vectors(
                carrier_data,
                AA_morphology_dict,
                chromophore_list,
                directory,
                sim_dims,
                current_carrier_type,
                args.three_D,
            )
            print("Calculating carrier trajectory anisotropy...")
            anisotropy = plot_anisotropy(
                carrier_data, directory, sim_dims, current_carrier_type, args.three_D
            )
            print("Plotting carrier hop frequency distribution...")
            if current_carrier_type == "hole":
                args.freq_cut_donor = plot_frequency_dist(
                    directory,
                    current_carrier_type,
                    carrier_history,
                    args.freq_cut_donor,
                )
            else:
                args.freq_cut_acceptor = plot_frequency_dist(
                    directory,
                    current_carrier_type,
                    carrier_history,
                    args.freq_cut_acceptor,
                )
            print("Plotting carrier net hop frequency distribution...")
            plot_net_frequency_dist(directory, current_carrier_type, carrier_history)
            print("Plotting ((total hops) - (net hops)) discrepancy distribution...")
            plot_discrepancy_frequency_dist(
                directory, current_carrier_type, carrier_history
            )
            if (carrier_history is not None) and args.three_D:
                print("Determining carrier hopping connections (network graph)...")
                plot_connections(
                    chromophore_list,
                    sim_dims,
                    carrier_history,
                    directory,
                    current_carrier_type,
                )
            if current_carrier_type == "hole":
                hole_anisotropy_data.append(anisotropy)
                hole_mobility_data.append([mobility, mob_error])
            elif current_carrier_type == "electron":
                electron_anisotropy_data.append(anisotropy)
                electron_mobility_data.append([mobility, mob_error])
            lower_carrier_type = current_carrier_type.lower()
            data_dict["name"] = os.path.split(directory)[1]
            data_dict["{:s}_anisotropy".format(lower_carrier_type)] = anisotropy
            data_dict["{:s}_mobility".format(lower_carrier_type)] = mobility
            data_dict["{:s}_mobility_r_squared".format(lower_carrier_type)] = r_squared
        # Now plot the distributions!
        temp_dir = os.path.join(directory, "figures")
        chromo_to_mol_ID = determine_molecule_IDs(
            CG_to_AAID_master, AA_morphology_dict, parameter_dict, chromophore_list
        )
        data_dict = plot_energy_levels(temp_dir, chromophore_list, data_dict)
        orientations = get_orientations(
            chromophore_list,
            CG_morphology_dict,
            AA_morphology_dict,
            CG_to_AAID_master,
            parameter_dict,
        )
        args.sep_cut_donor, args.sep_cut_acceptor = plot_neighbour_hist(
            chromophore_list,
            chromo_to_mol_ID,
            morphology_shape,
            temp_dir,
            args.sep_cut_donor,
            args.sep_cut_acceptor,
        )
        args.o_cut_donor, args.o_cut_acceptor = plot_orientation_hist(
            chromophore_list,
            chromo_to_mol_ID,
            orientations,
            temp_dir,
            args.o_cut_donor,
            args.o_cut_acceptor,
        )
        args.ti_cut_donor, args.ti_cut_acceptor = plot_TI_hist(
            chromophore_list,
            chromo_to_mol_ID,
            temp_dir,
            args.ti_cut_donor,
            args.ti_cut_acceptor,
        )
        cut_off_dict = create_cut_off_dict(
            args.sep_cut_donor,
            args.sep_cut_acceptor,
            args.o_cut_donor,
            args.o_cut_acceptor,
            args.ti_cut_donor,
            args.ti_cut_acceptor,
            args.freq_cut_donor,
            args.freq_cut_acceptor,
        )
        print("Cut-offs specified (value format: [donor, acceptor]) =", cut_off_dict)
        cluster_dicts, cluster_freqs, clusters_total, clusters_large, clusters_biggest, clusters_cutoffs = get_clusters(
            chromophore_list,
            carrier_history_dict,
            orientations,
            cut_off_dict,
            AA_morphology_dict,
        )
        if clusters_total[0] > 0:
            data_dict["donor_clusters_total"] = clusters_total[0]
            data_dict["donor_clusters_large"] = clusters_large[0]
            data_dict["donor_clusters_biggest"] = clusters_biggest[0]
            data_dict["donor_clusters_separation_cut"] = repr(clusters_cutoffs[0][0])
            data_dict["donor_clusters_orientation_cut"] = repr(clusters_cutoffs[0][1])
            data_dict["donor_clusters_transfer_integral_cut"] = repr(
                clusters_cutoffs[0][2]
            )
            data_dict["donor_clusters_hop_freq_cut"] = repr(clusters_cutoffs[0][3])
        if clusters_total[1] > 0:
            data_dict["acceptor_clusters_total"] = clusters_total[1]
            data_dict["acceptor_clusters_large"] = clusters_large[1]
            data_dict["acceptor_clusters_biggest"] = clusters_biggest[1]
            data_dict["acceptor_clusters_separation_cut"] = repr(clusters_cutoffs[1][0])
            data_dict["acceptor_clusters_orientation_cut"] = repr(
                clusters_cutoffs[1][1]
            )
            data_dict["acceptor_clusters_transfer_integral_cut"] = repr(
                clusters_cutoffs[1][2]
            )
            data_dict["acceptor_clusters_hop_freq_cut"] = repr(clusters_cutoffs[1][3])
        if args.three_D:
            print("Plotting 3D cluster location plot...")
            plot_clusters_3D(
                temp_dir, chromophore_list, cluster_dicts, sim_dims, args.generate_tcl
            )
        data_dict = plot_mixed_hopping_rates(
            temp_dir,
            chromophore_list,
            parameter_dict,
            cluster_dicts,
            chromo_to_mol_ID,
            data_dict,
            AA_morphology_dict,
            cut_off_dict,
        )
        print("Plotting cluster size distribution...")
        plot_cluster_size_dist(cluster_freqs, directory)
        print("Writing CSV Output File...")
        write_CSV(data_dict, directory)
    print("Plotting Mobility and Anisotropy progressions...")
    if args.sequence_donor is not None:
        if len(hole_anisotropy_data) > 0:
            plot_temperature_progression(
                args.sequence_donor,
                hole_mobility_data,
                hole_anisotropy_data,
                "hole",
                args.xlabel,
            )
    if args.sequence_acceptor is not None:
        if len(electron_anisotropy_data) > 0:
            plot_temperature_progression(
                args.sequence_acceptor,
                electron_mobility_data,
                electron_anisotropy_data,
                "electron",
                args.xlabel,
            )
    else:
        print("Skipping plotting mobility evolution.")


if __name__ == "__main__":
    main()
