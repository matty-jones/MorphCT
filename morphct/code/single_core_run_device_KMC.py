import copy
import heapq
import os
import pickle
import sys
import matplotlib.pyplot as plt
import numpy as np
import time as T

try:
    import mpl_toolkits.mplot3d as p3
except ImportError:
    print(
        "Could not import 3D plotting engine, calling the plot_device_components "
        "function"
        " will result in an error!"
    )
from morphct.code import helper_functions as hf


# Physical Constants
elementary_charge = 1.60217657E-19  # C
k_B = 1.3806488E-23  # m^{2} kg s^{-2} K^{-1}
hbar = 1.05457173E-34  # m^{2} kg s^{-1}
light_speed = 299792458  # ms^{-1}
epsilon_nought = 8.85418782E-12  # m^{-3} kg^{-1} s^{4} A^{2}
ELECTRON = 0
HOLE = 1

# Global Variables
global_chromophore_data = []  # Contains all of the chromophore_lists for each
#                                 morphology moiety for access by the hopping
#                                 routines
global_morphology_data = []  # Contains all of the AA_morphology_dicts for
#                                 each morphology moiety
global_time = 0  # The total simulation time that is
#                                 required/updated all over the place
global_carrier_dict = {}  # A dictionary of all carriers in the system for
#                                 the potential/recombination calculations
current_field_value = 0  # The field value calculated from the voltage
#                                 this child process was given
number_of_extractions = 0  # The total number of charges that have hopped
#                                 out of the device through the `correct'
#                                 contact
KMC_iterations = 0
# DEBUG Nothing is forbidden. Everything is permitted.
fastest_event_allowed = 1E-99  # 1E-15
slowest_event_allowed = 1E99  # 1E-9
log_file = None


class exciton:
    def __init__(self, index, global_time, initial_posn_device, parameter_dict):
        # Initialise the important variables to be used later
        self.ID = index
        self.creation_time = global_time
        self.removed_time = None
        self.T = parameter_dict["system_temperature"]
        self.lifetime_parameter = parameter_dict["exciton_lifetime"]
        # DEBUG excitons live forever
        self.recombination_time = -np.log(np.random.random()) * self.lifetime_parameter
        self.r_F = parameter_dict["forster_radius"]
        self.prefactor = parameter_dict["hopping_prefactor"]
        self.number_of_hops = 0
        # NOTE For now, we'll just inject it randomly somewhere in the system
        self.initial_device_posn = initial_posn_device
        self.current_device_posn = copy.deepcopy(initial_posn_device)
        self.initial_chromophore = global_chromophore_data.return_random_chromophore(
            initial_posn_device
        )
        self.current_chromophore = copy.deepcopy(self.initial_chromophore)
        self.calculate_behaviour()
        if parameter_dict["record_carrier_history"] is True:
            self.history = [[self.initial_device_posn, self.initial_chromophore.posn]]
        else:
            self.history = None

    def calculate_behaviour(self):
        # Set a flag to indicate whether exciton can dissociate or not
        self.can_dissociate = self.check_dissociation()
        # Dissociate the exciton immediately after creation if it would be
        # created at an interface
        if self.can_dissociate is True:
            # Determine all potential dissociation options, randomly select one,
            # plop a hole on the donor and electron on the acceptor, then remove
            # this exciton from the system by updating its removed_time
            if self.current_chromophore.species.lower() == "donor":
                h_chromo = global_chromophore_data.return_specific_chromophore(
                    self.current_device_posn, self.current_chromophore.ID
                )
                # The electron chromophore is a randomly selected chromophore of
                # the opposing type that is in range of the current one
                e_chromo_ID = np.random.choice(
                    self.current_chromophore.dissociation_neighbours
                )[0]
                e_chromo = global_chromophore_data.return_specific_chromophore(
                    self.current_device_posn, e_chromo_ID
                )
                if (self.current_device_posn not in h_chromo.occupied) and (
                    self.current_device_posn not in e_chromo.occupied
                ):
                    self.hole_chromophore = h_chromo
                    self.electron_chromophore = e_chromo
                else:
                    # hf.write_to_file(log_file, ["Debug: Cannot dissociate after all."
                    #                             " current_device_posn = "
                    #                             + repr(self.current_device_posn)
                    #                             + " h_chromo.occupied = "
                    #                             + str(h_chromo.occupied)
                    #                             + " e_chromo.occupied = "
                    #                             + str(e_chromo.occupied)])
                    self.can_dissociate = False
            elif self.current_chromophore.species.lower() == "acceptor":
                e_chromo = global_chromophore_data.return_specific_chromophore(
                    self.current_device_posn, self.current_chromophore.ID
                )
                # The hole chromophore is a randomly selected chromophore of the
                # opposing type that is in range of the current one
                h_chromo_ID = np.random.choice(
                    self.current_chromophore.dissociation_neighbours
                )[0]
                h_chromo = global_chromophore_data.return_specific_chromophore(
                    self.current_device_posn, h_chromo_ID
                )
                if (self.current_device_posn not in h_chromo.occupied) and (
                    self.current_device_posn not in e_chromo.occupied
                ):
                    self.hole_chromophore = h_chromo
                    self.electron_chromophore = e_chromo
                else:
                    # hf.write_to_file(log_file, ["Debug: Cannot dissociate after all."
                    #                             " current_device_posn = "
                    #                             + repr(self.current_device_posn)
                    #                             + " h_chromo.occupied = "
                    #                             + str(h_chromo.occupied)
                    #                             + "e_chromo.occupied = "
                    #                             + str(e_chromo.occupied)])
                    self.can_dissociate = False
            if self.can_dissociate is True:
                # Notify execute() that this exciton should not be queued up
                # again by setting self.hopTime == None
                [
                    self.destination_chromophore,
                    self.hop_time,
                    self.destination_image,
                ] = [None, None, None]
                self.removed_time = global_time
                return
        # If we're not instantaneously dissociating, check the recombination
        # time
        if global_time >= self.creation_time + self.recombination_time:
            # Exciton has run out of time and is now recombining, so set its
            # next hop_time to None
            [self.destination_chromophore, self.hop_time, self.destination_image] = [
                None,
                None,
                None,
            ]
            self.removed_time = global_time
            return
        # Otherwise, calculate the fastest hop from the current chromophore
        try:
            [
                self.destination_chromophore,
                self.hop_time,
                self.destination_image,
            ] = self.calculate_hop()
        except ValueError:
            # To remove the exciton from the system (or, more acurately, to
            # prevent the code from queueing up the exciton again because it has
            # already been popped from the main KMC queue), self.hopTime needs
            # to be None
            [self.destination_chromophore, self.hop_time, self.destination_image] = [
                None,
                None,
                None,
            ]
            self.removed_time = global_time

    def check_dissociation(self):
        # Return True if there are neighbours of the opposing electronic species
        # present, otherwise return False
        if len(self.current_chromophore.dissociation_neighbours) > 0:
            return True
        return False

    def calculate_hop(self):
        # First things first, if the current global time is > creationTime
        # + lifeTime then this exciton recombined before the hop so we need to
        # remove it from the system
        if global_time > (self.creation_time + self.recombination_time):
            return []

        # Get the chromophoreList so that we can work out which chromophores to
        # hop to
        chromophore_list = global_chromophore_data.return_chromophore_list(
            self.current_device_posn
        )
        # Determine the hop times to all possible neighbours
        hop_times = []
        for neighbour_index, transfer_integral in enumerate(
            self.current_chromophore.neighbours_TI
        ):
            # Ignore any hops with a NoneType transfer integral (usually due to
            # an orca error)
            if transfer_integral is None:
                continue
            # For the hop, we need to know the change in E_ij for the Boltzmann
            # factor, as well as the distance rij being hopped.
            # In Durham we used a prefactor of sqrt(2) here, can't remember why
            delta_E_ij = self.current_chromophore.neighbours_delta_E[neighbour_index]
            # The current posn inside the wrapped morphology is easy
            current_chromo_posn = self.current_chromophore.posn
            # The hop destination in wrapped morphology is slightly more
            # complicated
            neighbour_posn = chromophore_list[
                self.current_chromophore.neighbours[neighbour_index][0]
            ].posn
            # Need to determine the relative image between the two (recorded in
            # obtainChromophores). We will also need this to check if the
            # exciton is hopping out of this morphology cell and into an
            # adjacent one
            neighbour_relative_image = self.current_chromophore.neighbours[
                neighbour_index
            ][1]
            # Morphology shape so we can work out the actual relative positions
            morphology_shape = [
                global_morphology_data.return_AA_morphology(self.current_device_posn)[
                    key
                ]
                for key in ["lx", "ly", "lz"]
            ]
            remapped_posn = neighbour_posn + [
                neighbour_relative_image[axis] * morphology_shape[axis]
                for axis in range(3)
            ]
            rij = (
                hf.calculate_separation(
                    np.array(current_chromo_posn), np.array(remapped_posn)
                )
                * 1E-10
            )
            # Note, separations are recorded in angstroems, so spin this down to
            # metres. Additionally, all of the energies are in eV currently, so
            # convert them to J
            hop_rate = hf.calculate_FRET_hop_rate(
                self.prefactor,
                self.lifetime_parameter,
                self.r_F,
                rij,
                delta_E_ij * elementary_charge,
                self.T,
            )
            hop_time = hf.determine_event_tau(
                hop_rate,
                event_type="exciton-hop",
                slowest_event_allowed=slowest_event_allowed,
                fastest_event=fastest_event_allowed,
                maximum_attempts=100,
                log_file=log_file,
            )
            # Keep track of the destination chromophore ID, the corresponding
            # tau, and the relative image (to see if we're hopping over a
            # boundary)
            if hop_time is not None:
                hop_times.append(
                    [
                        global_chromophore_data.return_specific_chromophore(
                            self.current_device_posn,
                            self.current_chromophore.neighbours[neighbour_index][0],
                        ),
                        hop_time,
                        neighbour_relative_image,
                    ]
                )
        # Sort by ascending hop time
        hop_times.sort(key=lambda x: x[1])
        # Only want to take the quickest hop if it is NOT going to hop outside
        # the device (top or bottom - i.e. z axis is > shape[2] or < 0)
        # Other-axis periodic hops (i.e. x and y) are permitted, and will allow
        # the exciton to loop round again.
        while len(hop_times) > 0:
            # Get the hop destination cell
            new_location = list(
                np.array(self.current_device_posn) + np.array(hop_times[0][2])
            )
            # If the z axis component is > shape[2] or < 0, then forbid this hop
            # by removing it from the hopTimes list.
            if (new_location[2] >= global_chromophore_data.device_array.shape[2]) or (
                new_location[2] < 0
            ):
                hop_times.pop(0)
            else:
                break
        # Take the quickest hop
        if len(hop_times) > 0:
            return hop_times[0]
        else:
            # The exciton has not yet recombined, but there are no legal hops
            # that can be performed so the only fate for this exciton is
            # recombination through photoluminescence
            return []

    def perform_hop(self):
        destination_position = self.destination_chromophore.posn
        if self.destination_image == [0, 0, 0]:
            # Exciton is not hopping over a boundary, so can simply update its
            # current position
            self.current_chromophore = self.destination_chromophore
        else:
            # We're hopping over a boundary. Permit the hop with the
            # already-calculated hopTime, but then find the closest chromophore
            # to the destination position in the adjacent device cell.
            target_device_posn = list(
                np.array(self.current_device_posn) + np.array(self.destination_image)
            )
            new_chromo = global_chromophore_data.return_closest_chromophore_to_position(
                target_device_posn, destination_position
            )
            if (new_chromo.lower() == "top") or (new_chromo.lower() == "bottom"):
                # This exciton is hopping out of the active layer and into the
                # contacts.
                # Ensure it doesn't get queued up again
                [
                    self.destination_chromophore,
                    self.hop_time,
                    self.destination_image,
                ] = [None, None, None]
                self.removed_time = global_time
                # TODO This carrier has left the simulation volume so we now
                # need to ensure we remove it from the carrier dictionary so
                # it's not included in any of the Coulombic calculations
            elif new_chromo.lower() == "out of bounds":
                # Ignore this hop, do a different one instead
                pass
            else:
                self.current_device_posn = list(
                    np.array(self.current_device_posn)
                    + np.array(self.destination_image)
                )
                self.current_chromophore = new_chromo
        self.number_of_hops += 1
        self.can_dissociate = self.check_dissociation()
        if self.history is not None:
            self.history.append(
                [self.current_device_posn, self.current_chromophore.posn]
            )


class carrier:
    def __init__(
        self,
        index,
        global_time,
        initial_posn_device,
        initial_chromophore,
        injected_from,
        parameter_dict,
        injected_onto_site=None,
    ):
        self.ID = index
        self.creation_time = global_time
        self.removed_time = None
        self.initial_device_posn = initial_posn_device
        self.current_device_posn = initial_posn_device
        self.initial_chromophore = initial_chromophore
        self.current_chromophore = initial_chromophore
        self.injected_from = injected_from
        self.injected_onto_site = injected_onto_site
        self.prefactor = parameter_dict["hopping_prefactor"]
        self.recombining = False
        self.recombining_with = None
        self.T = parameter_dict["system_temperature"]
        self.wrapxy = parameter_dict["wrap_device_xy"]
        self.disable_coulombic = parameter_dict["disable_coulombic"]
        if self.current_chromophore.sub_species == self.initial_chromophore.sub_species:
            self.lambda_ij = self.current_chromophore.reorganisation_energy
        else:
            self.lambda_ij = (
                self.current_chromophore.reorganisation_energy
                + self.initial_chromophore.reorganisation_energy
            ) / 2
        # self.carrier_type: Set carrier type to be == 0 if Electron and == 1 if
        # Hole. This allows us to do quick arithmetic to get the signs correct
        # in the potential calculations without having to burn through a ton of
        # conditionals.
        if self.current_chromophore.species.lower() == "donor":
            self.carrier_type = HOLE
        elif self.current_chromophore.species.lower() == "acceptor":
            self.carrier_type = ELECTRON
        self.relative_permittivity = parameter_dict["relative_permittivity"]
        if parameter_dict["record_carrier_history"] is True:
            self.history = [[self.initial_device_posn, initial_chromophore.posn]]
        else:
            self.history = None
        # Set the use of Koopmans' approximation to false if the key does not
        # exist in the parameter dict
        try:
            self.use_koopmans_approximation = parameter_dict[
                "use_koopmans_approximation"
            ]
        except KeyError:
            self.use_koopmans_approximation = False
        # Are we using a simple Boltzmann penalty?
        try:
            self.use_simple_energetic_penalty = parameter_dict[
                "use_simple_energetic_penalty"
            ]
        except KeyError:
            self.use_simple_energetic_penalty = False
        # Are we applying a distance penalty beyond the transfer integral?
        try:
            self.use_VRH = parameter_dict["use_VRH"]
        except KeyError:
            self.use_VRH = False
        if self.use_VRH is True:
            self.VRH_scaling = 1.0 / self.current_chromophore.VRH_delocalisation

    def calculate_behaviour(self):
        try:
            hop_data = self.calculate_hop()
            self.destination_chromophore = hop_data[0]
            self.hop_time = hop_data[1]
            self.relative_image = hop_data[2]
            self.destination_image = hop_data[3]
            # Update the destination chromophore (if it's a chromophore) so that
            # it's marked as occupied
            if not isinstance(self.destination_chromophore, str):
                global_chromophore_data.return_specific_chromophore(
                    self.destination_image, self.destination_chromophore.ID
                ).occupied.append(self.destination_image)
        except ValueError:
            self.destination_chromophore = None
            self.hop_time = None
            self.relative_image = None
            self.destination_image = None

    def calculate_hop(self):
        # Determine the hop times to all possible neighbours
        hop_times = []
        # Obtain the reorganisation energy in J (from eV in the parameter file)
        for neighbour_index, transfer_integral in enumerate(
            self.current_chromophore.neighbours_TI
        ):
            # Ignore any hops with a NoneType transfer integral (usually due to
            # an orca error), or zero
            if (transfer_integral is None) or (transfer_integral < 1E-10):
                continue
            neighbour_chromophore = global_chromophore_data.return_specific_chromophore(
                self.current_device_posn,
                self.current_chromophore.neighbours[neighbour_index][0],
            )
            # The destination chromophore will be the actual chromophore we end
            # up on (i.e. not neighbour if we hop across a boundary)
            destination_chromophore = neighbour_chromophore
            # Need to determine the relative image between the two (recorded in
            # obtainChromophores) to check if the carrier is hopping out of this
            # morphology cell and into an adjacent one
            neighbour_relative_image = self.current_chromophore.neighbours[
                neighbour_index
            ][1]
            destination_image = list(
                np.array(self.current_device_posn) + np.array(neighbour_relative_image)
            )
            # If we're hopping out of this cell and into a new one, make sure
            # that it's in-range if we're not going to wrap it
            skip_this_neighbour = False
            # Skip this neighbour if the destination is already occupied
            if destination_image in destination_chromophore.occupied:
                continue
            elif self.wrapxy is False:
                # Skip this neighbour if we're going to hop out of the device
                # along the non-electrode axes
                if neighbour_relative_image != [0, 0, 0]:
                    for axis_no, val in enumerate(destination_image[:-1]):
                        if (
                            val >= global_chromophore_data.device_array.shape[axis_no]
                        ) or (val < 0):
                            skip_this_neighbour = True
                            break
                    # Don't need to perform this check if we've already
                    # satisfied the first skip condition
                    if skip_this_neighbour is False:
                        destination_chromophore = global_chromophore_data.return_closest_chromophore_to_position(
                            destination_image, neighbour_chromophore.posn
                        )
                        if destination_chromophore.lower() == "out of bounds":
                            continue
                    else:
                        continue
            # Make sure the destination chromophore is of the correct type and
            # is not occupied, if it is a chromophore and not a string (saying
            # 'Top' or 'Bottom' if it's leaving the device)
            if not isinstance(destination_chromophore, str):
                if (destination_image in destination_chromophore.occupied) or (
                    destination_chromophore.species != self.current_chromophore.species
                ):
                    continue
            # Otherwise, we're good to go. Calculate the hop as previously (but
            # a hop to the neighbourChromophore)
            delta_E_ij = self.calculate_delta_E(
                neighbour_chromophore,
                neighbour_relative_image,
                self.current_chromophore.neighbours_delta_E[neighbour_index],
            )
            # All of the energies (EXCEPT EIJ WHICH IS ALREADY IN J) are in eV
            # currently, so convert them to J
            if self.use_VRH is True:
                relative_image = np.array(destination_image) - np.array(
                    self.current_device_posn
                )
                destination_chromo_posn = destination_chromo.posn + (
                    np.array(relative_image)
                    * np.array([axis[1] - axis[0] for axis in self.sim_dims])
                )
                # Convert from ang to m
                chromophore_separation = (
                    hf.calculate_separation(
                        self.current_chromophore.posn, destination_chromo_posn
                    )
                    * 1E-10
                )
                hop_rate = hf.calculate_carrier_hop_rate(
                    self.lambda_ij * elementary_charge,
                    transfer_integral * elementary_charge,
                    delta_E_ij,
                    self.prefactor,
                    self.T,
                    use_VRH=True,
                    rij=chromophore_separation,
                    VRH_prefactor=self.VRH_scaling,
                    boltz_pen=self.use_simple_energetic_penalty,
                )
            else:
                hop_rate = hf.calculate_carrier_hop_rate(
                    self.lambda_ij * elementary_charge,
                    transfer_integral * elementary_charge,
                    delta_E_ij,
                    self.prefactor,
                    self.T,
                    boltz_pen=self.use_simple_energetic_penalty,
                )
            hop_time = hf.determine_event_tau(
                hop_rate,
                event_type="carrier-hop",
                slowest_event_allowed=slowest_event_allowed,
                fastest_event=fastest_event_allowed,
                maximum_attempts=100,
                log_file=log_file,
            )
            if hop_time is not None:
                # Keep track of the chromophoreID and the corresponding tau
                hop_times.append(
                    [
                        destination_chromophore,
                        hop_time,
                        neighbour_relative_image,
                        destination_image,
                        self.prefactor,
                        self.lambda_ij * elementary_charge,
                        transfer_integral * elementary_charge,
                        delta_E_ij,
                    ]
                )
        # Sort by ascending hop time
        hop_times.sort(key=lambda x: x[1])
        # Take the quickest hop
        if len(hop_times) > 0:
            if hop_times[0][1] > 1:
                hf.write_to_file(
                    log_file,
                    [
                        "---=== WARNING, selected extremely long carrier hop ===---",
                        repr(hop_times[-1]),
                    ],
                )
                hf.write_to_file(
                    log_file, [repr(hop_times[0]), repr(global_carrier_dict)]
                )
                for ID, carrier in global_carrier_dict.items():
                    hf.write_to_file(
                        log_file,
                        [
                            "".join(
                                [
                                    "ID = {:d}, current posn =".format(ID),
                                    repr(carrier.current_device_posn),
                                    ", current_chromophore_ID = {:d}".format(
                                        carrier.current_chromophore.ID
                                    ),
                                ]
                            )
                        ],
                    )
                print("NO VIABLE HOPS FOUND")
                return []
            return hop_times[0][:4]
        else:
            # We are trapped here, so just return in order to raise the
            # calculate_behaviour exception
            print("NO VIABLE HOPS FOUND")
            return []

    def perform_hop(self):
        global number_of_extractions
        global KMC_iterations
        global global_time
        # Unset the current chromophore's occupation status as we're about to
        # hop away. The destination should have already been updated
        global_chromophore_data.return_specific_chromophore(
            self.current_device_posn, self.current_chromophore.ID
        ).occupied.remove(self.current_device_posn)
        if (self.destination_chromophore.lower() == "top") or (
            self.destination_chromophore.lower() == "bottom"
        ):
            # This carrier is hopping out of the active layer and into the
            # contacts. Firstly, work out whether this is a `correct' hop
            # (i.e. hole hopping to anode or electron hopping to cathode)
            # that causes photovoltaic current.
            if self.destination_chromophore.lower() == "top":
                # Leaving through top (anode)
                if self.injected_from.lower() != "anode":
                    if self.carrier_type == HOLE:
                        number_of_extractions += 1
                    else:
                        number_of_extractions -= 1
                # Else (injected from anode), number of extractions doesn't
                # change.
            else:
                # Leaving through bottom (cathode)
                if self.injected_from.lower() != "cathode":
                    if self.carrier_type == ELECTRON:
                        number_of_extractions += 1
                    else:
                        number_of_extractions -= 1
                # Else (injected from cathode), number of extractions
                # doesn't change.
            if self.carrier_type == ELECTRON:
                hf.write_to_file(
                    log_file,
                    [
                        "".join(
                            [
                                "EVENT: Electron left out of {0:s} of device! New number of extractions: {0:d} after {1:d} iterations (global_time = {2:.2e})".format(
                                    self.destination_chromophore,
                                    number_of_extractions,
                                    KMC_iterations,
                                    global_time,
                                )
                            ]
                        )
                    ],
                )
            else:
                hf.write_to_file(
                    log_file,
                    [
                        "".join(
                            [
                                "EVENT: Hole left out of {0:s} of device! New number of extractions: {0:d} after {1:d} iterations (global_time = {2:.2e})".format(
                                    self.destination_chromophore,
                                    number_of_extractions,
                                    KMC_iterations,
                                    global_time,
                                )
                            ]
                        )
                    ],
                )
            # Now ensure that it doesn't get queued up again
            self.destination_chromophore = None
            self.hop_time = None
            self.relative_image = None
            self.destination_image = None
            self.removed_time = global_time
        else:
            self.current_device_posn = self.destination_image
            self.current_chromophore = self.destination_chromophore
        if self.history is not None:
            self.history.append(
                [self.current_device_posn, self.current_chromophore.posn]
            )

    def calculate_delta_E(
        self, destination_chromophore, neighbour_relative_image, chromo_E_ij
    ):
        # Delta_E_ij has 3 main components: 1) the energetic disorder (difference
        # in HOMO/LUMO levels), 2) the field within the device, and 3) the
        # Coulombic effect from nearby charges
        delta_E_ij = 0.0  # report this in J
        # 1) Energetic Disorder
        delta_E_ij += chromo_E_ij * elementary_charge
        # 2) Field within the device
        current_absolute_position = np.array(self.current_device_posn) * parameter_dict[
            "morphology_cell_size"
        ] + (np.array(self.current_chromophore.posn) * 1E-10)
        destination_absolute_position = (
            np.array(self.current_device_posn) + np.array(neighbour_relative_image)
        ) * parameter_dict["morphology_cell_size"] + (
            np.array(destination_chromophore.posn) * 1E-10
        )
        # Field has negative sign because device is flipped with anode at +Z and
        # Cathode at 0
        z_sep = -(destination_absolute_position[2] - current_absolute_position[2])
        charge = elementary_charge * ((2 * self.carrier_type) - 1)
        delta_E_ij += z_sep * current_field_value * charge
        if self.disable_coulombic is not True:
            # 3) Difference in Coulombic Potential at dest compared to origin
            origin_coulomb, recombine_flag, recombine_ID = calculate_coulomb(
                current_absolute_position,
                self.ID,
                self.carrier_type,
                carrier_is_recombining=self.recombining,
            )
            destination_coulomb, dummy1, dummy2 = calculate_coulomb(
                destination_absolute_position, self.ID, self.carrier_type
            )
            delta_E_ij += destination_coulomb - origin_coulomb
            if (recombine_flag is True) and (self.recombining is False):
                # Need to update the recombining details
                self.recombining = recombine_flag
                self.recombining_with = recombine_ID
        return delta_E_ij


class inject_site:
    def __init__(self, device_posn, chromophore, inject_rate, electrode):
        self.device_posn = device_posn
        self.chromophore = chromophore
        self.inject_rate = inject_rate
        self.electrode = electrode
        self.calculate_inject_time()

    def calculate_inject_time(self):
        self.inject_time = hf.determine_event_tau(
            self.inject_rate,
            event_type="".join([self.electrode, "-injection"]),
            slowest_event_allowed=slowest_event_allowed,
            fastest_event=fastest_event_allowed,
            maximum_attempts=100,
            log_file=log_file,
        )


def plot_hop_distance(distribution):
    plt.figure()
    plt.hist(distribution)
    plt.show()
    exit()


def plot_device_components(device_array):
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    colour = ["r", "g", "b", "y", "k"]
    for x_val in range(9):
        for y_val in range(9):
            for z_val in range(9):
                ax.scatter(
                    x_val,
                    y_val,
                    z_val,
                    zdir="z",
                    c=colour[device_array[x_val, y_val, z_val]],
                )
    plt.show()
    exit()


def calculate_photoinjection_rate(parameter_dict, device_shape):
    # Photoinjection rate is given by the following equation. Calculations will
    # be performed in SI always.
    # Flux is multiplied by 10 to comvert from mW/cm^{2} to W/m^{2}
    rate = (
        (parameter_dict["incident_flux"] * 10)
        * (parameter_dict["incident_wavelength"] / (hbar * 2 * np.pi * light_speed))
        * device_shape[0]
        * device_shape[1]
        * parameter_dict["morphology_cell_size"] ** 2
        * (
            1
            - np.exp(
                -100
                * parameter_dict["absorption_coefficient"]
                * device_shape[2]
                * parameter_dict["morphology_cell_size"]
            )
        )
    )
    return rate


def decrement_time(event_queue, event_time):
    global global_time
    # A function that increments the global time by whatever the time this event
    # is, and then decrements all of the remaining times in the queue.
    global_time += event_time
    event_queue = [(event[0] - event_time, event[1], event[2]) for event in event_queue]
    return event_queue


def plot_carrier_Z_profiles(all_carriers, parameter_dict, device_array, output_dir):
    carriers_to_plot = []
    # Only consider electrons and holes
    for carrier in all_carriers:
        if "r_F" not in carrier.__dict__:
            carriers_to_plot.append(carrier)
    # Determine the ylims of the zProfile plot
    z_len = 1E10 * (
        np.array(device_array.shape[2]) * parameter_dict["morphology_cell_size"]
    )
    # Now make the plots
    for carrier in carriers_to_plot:
        if (carrier.injected_from.lower() == "anode") or (
            carrier.injected_from.lower() == "cathode"
        ):
            continue
        plot_Z_profile(carrier, z_len, output_dir)


def plot_Z_profile(carrier, z_dim_size, output_dir):
    x_vals = []
    y_vals = []
    for hop_index, hop in enumerate(carrier.history):
        x_vals.append(hop_index)
        current_Z = (
            1E10
            * ((np.array(hop[0][2]) + 0.5))
            * parameter_dict["morphology_cell_size"]
        ) + np.array(hop[1][2])
        y_vals.append(current_Z)
        if carrier.carrier_type == ELECTRON:
            colour = "b"
        else:
            colour = "r"
            file_name = "{0:s}carrier_{1:05d}_Z_profile.pdf".format(
                output_dir, carrier.ID
            )
    plt.figure()
    plt.plot(x_vals, y_vals, color=colour)
    plt.ylim([0, z_dim_size])
    plt.xlabel("Hop Number (Arb. U.)")
    plt.ylabel("Z-position (Ang)")
    plt.savefig(file_name)
    print("".join(["Carrier Z-Profile saved as ", fileName]))
    plt.close()


def calculate_coulomb(
    absolute_position, self_ID, self_carrier_type, carrier_is_recombining=None
):
    global global_carrier_dict
    coulombic_potential = 0.0
    coulomb_constant = 1.0 / (
        4 * np.pi * epsilon_nought * parameter_dict["relative_permittivity"]
    )
    can_recombine_here = False
    recombining_with = None
    for carrier_ID, carrier in global_carrier_dict.items():
        # Each carrier feels coulombic effects from every other carrier in the
        # system, along with the image charges induced in the top and bottom
        # contacts, AND the image charges induced by the image charges.
        # The image charges are needed in order to keep a constant field at the
        # interface, correcting for the non-periodic Neumann-like boundary
        # conditions. Check papers from Chris Groves (Marsh, Groves and
        # Greenham2007 and beyond), Ben Lyons (2011/2012), and van der Holst
        # (2011) for more details.
        carrier_posn = np.array(carrier.current_device_posn) * parameter_dict[
            "morphology_cell_size"
        ] + (np.array(carrier.current_chromophore.posn) * 1E-10)
        if carrier_ID != self_ID:
            # Only consider the current carrier if it is not us!
            separation = hf.calculate_separation(absolute_position, carrier_posn)
            if separation == 0.0:
                print("ZERO SEPARATION BETWEEN CARRIERS, EXITING")
                if self_carrier_type == carrier.carrier_type:
                    print("ZERO SEPARATION BETWEEN LIKE CARRIERS, EXITING")
                    exit()
                exit()
            else:
                coulombic_potential += coulomb_constant * (
                    (elementary_charge * ((2 * self_carrier_type) - 1))
                    * (elementary_charge * ((2 * carrier.carrier_type) - 1))
                    / separation
                )
            # I'm also going to use this opportunity (as we're iterating over
            # all carriers in the system) to see if we're close enough to any to
            # recombine.
            if carrier_is_recombining is False:
                # Only do this if we're not currently recombining with something
                if (separation <= parameter_dict["coulomb_capture_radius"]) and (
                    self_carrier_type != carrier.carrier_type
                ):
                    # If carriers are within 1nm of each other, then assume that
                    # they are within the Coulomb capture radius and are about
                    # to recombine.
                    can_recombine_here = True
                    recombining_with = carrier_ID
                    # We don't necessarily need to update the
                    # carrier.recombining, since all of the behaviour will be
                    # determined in the main program by checking
                    # hoppingCarrier.recombining, but just in case, we can
                    # update the carrier too as well as self.
                    global_carrier_dict[carrier_ID].recombining = True
                    global_carrier_dict[carrier_ID].recombining_with = self_ID
        # TODO: NOTE: DEBUG: TURNING OFF IMAGE CHARGES TO SEE IF I CAN WORK OUT
        # WHY EVERYTHING IS FAILING
        # # Now consider the image charges and the images of the images.
        # # Direct images have the opposing charge to the actual carrier, so the
        # # potential contains (+q * -q), or -(q**2). Images of the images have
        # # the opposing charge to the image of the carrier, which is the same
        # # charge as the actual carrier. Therefore, images of images are either
        # # (+q * +q) or (-q * -q) = +(q**2)
        # # First do the top image charge
        # deviceZSize = globalChromophoreData.deviceArray.shape[2]\
        #                   * parameterDict['morphologyCellSize']
        # topImagePosn = copy.deepcopy(carrierPosn)
        # topImagePosn[2] = deviceZSize + (deviceZSize - carrierPosn[2])
        # separation = hf.calculateSeparation(absolutePosition, topImagePosn)
        # print("Separation between", selfID, "and top image charge =", separation)
        # coulombicPotential -= (elementaryCharge**2) * (coulombConstant / separation)
        # # Now do the bottom image charge
        # bottomImagePosn = copy.deepcopy(carrierPosn)
        # bottomImagePosn[2] = - carrierPosn[2]
        # separation = hf.calculateSeparation(absolutePosition, bottomImagePosn)
        # print("Separation between", selfID, "and bottom image charge =", separation)
        # coulombicPotential -= (elementaryCharge**2) * (coulombConstant / separation)
        # # Now do the top image of the bottom image charge
        # topImageOfBottomImagePosn = copy.deepcopy(carrierPosn)
        # topImageOfBottomImagePosn[2] = (2 * deviceZSize) + carrierPosn[2]
        # separation = hf.calculateSeparation(absolutePosition,
        #                                     topImageOfBottomImagePosn)
        # print("Separation between", selfID, "and top image of bottom image charge =",
        #       separation)
        # coulombicPotential += (elementaryCharge**2) * (coulombConstant / separation)
        # # And finally the bottom image of the top image charge
        # bottomImageOfTopImagePosn = copy.deepcopy(carrierPosn)
        # bottomImageOfTopImagePosn[2] = - (deviceZSize + (deviceZSize
        #                                                  - carrierPosn[2]))
        # separation = hf.calculateSeparation(absolutePosition,
        #                                     bottomImageOfTopImagePosn)
        # print("Separation between", selfID, "and bottom image of top image charge =",
        #       separation, "\n")
        # coulombicPotential += (elementaryCharge**2) * (coulombConstant / separation)
    return coulombic_potential, can_recombine_here, recombining_with


def calculate_dark_current_injections(device_array, parameter_dict):
    # NOTE Had a nightmare with Marcus hopping here, and I'm not convinced it
    # makes sense to have hops from metals be described wrt reorganisation
    # energies and transfer integrals. Instead, I'm going to use a more generic
    # MA hopping methodology here. This therefore needs a prefactor, which is
    # defined in the parameter file. This is another variable to calibrate, but
    # for now I've set it such that we get a reasonable hopping rate for
    # dark-injection hops (somewhere in the middle of our normal hopping-rate
    # distribution (i.e. ~1E12)). Get the device shape so that we know which
    # cells to inject into
    device_shape = device_array.shape
    # Calculate the important energy levels
    bandgap = parameter_dict["acceptor_LUMO"] - parameter_dict["donor_HOMO"]
    electron_inject_barrier = (
        parameter_dict["acceptor_LUMO"] - parameter_dict["cathode_work_function"]
    )
    hole_inject_barrier = (
        parameter_dict["anode_work_function"] - parameter_dict["donor_HOMO"]
    )
    valid_cathode_inj_sites = 0
    valid_anode_inj_sites = 0
    cathode_inject_rates_data = []
    anode_inject_rates_data = []
    cathode_inject_wait_times = []
    anode_inject_wait_times = []
    # Consider electron-injecting electrode (cathode) at bottom of device first:
    z_val = 0
    for x_val in range(device_shape[0]):
        for y_val in range(device_shape[1]):
            cathode_inject_chromophores = []
            AA_morphology = global_morphology_data.return_AA_morphology(
                [x_val, y_val, z_val]
            )
            morphology_chromophores = global_chromophore_data.return_chromophore_list(
                [x_val, y_val, z_val]
            )
            for chromophore in morphology_chromophores:
                # Find all chromophores that are between 5 and 10 Ang from the
                # bottom 10 Ang of the device cell
                # if (chromophore.posn[2] <= -(AA_morphology['lz'] / 2.0) + 10)\
                #    and (sum(1 for _ in filter(
                #    None.__ne__, chromophore.neighbours_delta_E)) > 0):
                if (
                    (chromophore.posn[2] <= -(AA_morphology["lz"] / 2.0) + 10)
                    and (chromophore.posn[2] >= -(AA_morphology["lz"] / 2.0) + 5)
                    and (
                        len(
                            [
                                _
                                for _ in chromophore.neighbours_TI
                                if (_ is not None) and (_ > 1E-5)
                            ]
                        )
                        > 0
                    )
                ):
                    cathode_inject_chromophores.append(
                        [
                            chromophore,
                            1E-10
                            * (chromophore.posn[2] - (-(AA_morphology["lz"] / 2.0))),
                        ]
                    )
                    valid_cathode_inj_sites += 1
            for [chromophore, separation] in cathode_inject_chromophores:
                if chromophore.species.lower() == "acceptor":
                    # Injecting an electron from the cathode (easy)
                    delta_E = elementary_charge * (electron_inject_barrier)
                else:
                    # Injecting a hole from the cathode (hard)
                    # Here, we make the assumption that the energy difference
                    # between the chromophore and literature HOMOs would be the
                    # same as the energy difference between the chromophore and
                    # literature LUMOs for injection (this is because we don't
                    # record the donor LUMO or the acceptor HOMO)
                    delta_E = elementary_charge * (bandgap - electron_inject_barrier)
                inject_rate = hf.calculate_miller_abrahams_hop_rate(
                    parameter_dict["MA_prefactor"],
                    separation,
                    parameter_dict["MA_localisation_radius"],
                    delta_E,
                    parameter_dict["system_temperature"],
                )
                cathode_inject_rates_data.append(inject_rate)
                # Create inject site object
                site = inject_site(
                    [x_val, y_val, z_val], chromophore, inject_rate, "cathode"
                )
                inject_time = site.inject_time
                cathode_inject_wait_times = push_to_queue(
                    cathode_inject_wait_times, (inject_time, "cathode-injection", site)
                )
    # Now consider hole-injecting electrode at top of device (anode):
    z_val = device_shape[2] - 1
    for x_val in range(device_shape[0]):
        for y_val in range(device_shape[1]):
            anode_inject_chromophores = []
            AA_morphology = global_morphology_data.return_AA_morphology(
                [x_val, y_val, z_val]
            )
            morphology_chromophores = global_chromophore_data.return_chromophore_list(
                [x_val, y_val, z_val]
            )
            for chromophore in morphology_chromophores:
                # Find all chromophores that are between 5 and 10 Ang from the
                # top of the device cell
                # if (chromophore.posn[2] >= (AA_morphology['lz'] / 2.0) - 10)\
                #    and (sum(1 for _ in filter(
                #    None.__ne__, chromophore.neighbours_delta_E)) > 0):
                if (
                    (chromophore.posn[2] >= (AA_morphology["lz"] / 2.0) - 10)
                    and (chromophore.posn[2] <= (AA_morphology["lz"] / 2.0) - 5)
                    and (
                        len(
                            [
                                _
                                for _ in chromophore.neighbours_TI
                                if (_ is not None) and (_ > 1E-5)
                            ]
                        )
                        > 0
                    )
                ):
                    anode_inject_chromophores.append(
                        [
                            chromophore,
                            1E-10 * ((AA_morphology["lz"] / 2.0) - chromophore.posn[2]),
                        ]
                    )
                    valid_anode_inj_sites += 1
            for [chromophore, separation] in anode_inject_chromophores:
                if chromophore.species.lower() == "acceptor":
                    # Injecting an electron from the anode (hard)
                    delta_E = elementary_charge * (bandgap - hole_inject_barrier)
                else:
                    # Injecting a hole from the anode (easy)
                    delta_E = elementary_charge * (hole_inject_barrier)
                inject_rate = hf.calculate_miller_abrahams_hop_rate(
                    parameter_dict["MA_prefactor"],
                    separation,
                    parameter_dict["MA_localisation_radius"],
                    delta_E,
                    parameter_dict["system_temperature"],
                )
                anode_inject_rates_data.append(inject_rate)
                # Create inject site object
                site = inject_site(
                    [x_val, y_val, z_val], chromophore, inject_rate, "anode"
                )
                inject_time = site.inject_time
                anode_inject_wait_times = push_to_queue(
                    anode_inject_wait_times, (inject_time, "anode-injection", site)
                )
    # plt.figure()
    # plt.hist(injectRates, bins = np.logspace(9, 12, 30))
    # plt.gca().set_xscale('log')
    # print("Valide Cathode Injection Sites =", validCathodeInjSites)
    # print("Valide Anode Injection Sites =", validAnodeInjSites)
    # print("Mean =", np.mean(injectRates), "std =", np.std(injectRates))
    # print("Max =", max(injectRates), "Min =", min(injectRates))
    # plt.show()
    # exit()
    cathode_inject_rate = np.mean(cathode_inject_rates_data)
    anode_inject_rate = np.mean(anode_inject_rates_data)
    return [
        cathode_inject_rate,
        anode_inject_rate,
        cathode_inject_wait_times,
        anode_inject_wait_times,
    ]


def get_next_dark_event(queue, electrode):
    # A function that pops the next event from a queue (cathode or anode
    # injection queues), decrements the time, requeues the inject_site with a
    # new time, pushes to the queue and returns the next event and the updated
    # queue
    next_event = heapq.heappop(queue)
    queue = decrement_time(queue, next_event[0])
    requeued_event = copy.deepcopy(next_event[2])
    requeued_event.calculate_inject_time()
    queue = push_to_queue(
        queue,
        (
            requeued_event.inject_time,
            "".join([electrode, "-injection"]),
            requeued_event,
        ),
    )
    return next_event, queue


def estimate_transfer_integral(chromophore_list):
    # Iterate over all non-periodic chromophore hopping events in the system
    TIs = []
    for chromophore in chromophore_list:
        for index, chromo in enumerate(chromophore.neighbours):
            if chromo[1] != [0, 0, 0]:
                continue
            TI = chromophore.neighbours_TI[index]
            if TI is not None:
                TIs.append(TI)
    return np.average(TIs)


def push_to_queue(queue, event):
    event = list(event)
    try:
        event[0] = np.float64(event[0])
    except:
        hf.write_to_file(
            log_file,
            [
                "Tried to push an event to the queue that does not have a"
                " np.float64 as the first element (and so can't be queued)",
                "".join(["Queue = ", repr(queue)]),
                "".join(["Event = ", repr(event)]),
                "Terminating...",
            ],
        )
        raise KeyboardInterrupt
    if (event[0] == np.float64(1E99)) and (event[1].lower() != "photo"):
        log_line_to_write = [
            "---=== TRIED TO QUEUE EVENT WITH CRAZY LONG WAIT TIME ===---"
            "".join([" Event = ", event])
        ]
        if event[2].history is not None:
            log_line_to_write.append(
                "This carrier has completed {:d} hops.".format(len(event[2].history))
            )
        hf.write_to_file(log_file, log_line_to_write)
        try:
            event[2].__dict__.pop("history")
        except KeyError:
            pass
        hf.write_to_file(log_file, ["".join(["Instance = ", repr(event[2].__dict__)])])
        try:
            hf.write_to_file(
                log_file,
                [
                    "".join(
                        [
                            "Current chromophore = ",
                            repr(event[2].current_chromophore.__dict__),
                        ]
                    ),
                    "".join(
                        [
                            "Destination chromophore = ",
                            repr(event[2].destination_chromophore.__dict__),
                        ]
                    ),
                ],
            )
        except AttributeError:
            pass
        hf.write_to_file(log_file, ["Terminating..."])
        for index, carrier in global_carrier_dict.items():
            hf.write_to_file(
                log_file,
                [
                    " ".join(
                        [
                            str(carrier.current_device_posn),
                            str(carrier.current_chromophore.posn),
                        ]
                    )
                ],
            )
        for carrier1ID, carrier1 in global_carrier_dict.items():
            for carrier2ID, carrier2 in global_carrier_dict.items():
                if carrier1ID >= carrier2ID:
                    continue
                carrier1posn = (78 * np.array(carrier1.current_device_posn)) + np.array(
                    carrier1.current_chromophore.posn
                )
                carrier2posn = (78 * np.array(carrier2.current_device_posn)) + np.array(
                    carrier2.current_chromophore.posn
                )
                hf.write_to_file(
                    log_file,
                    [
                        "{0:d} {1:d} {2:f}".format(
                            carrier1ID,
                            carrier2ID,
                            hf.calculate_separation(carrier1posn, carrier2posn),
                        )
                    ],
                )
        raise keyboard_interrupt
    event = tuple(event)
    heapq.heappush(queue, event)
    return queue


def plot_event_time_distribution(event_log, output_dir, fastest, slowest):
    plt.figure()
    plt.hist(
        [event_log[event_type] for event_type in sorted(event_log.keys())],
        bins=np.logspace(
            int(np.floor(np.log10(fastest))), int(np.ceil(np.log10(slowest))), 10
        ),
        color=["r", "g", "c", "m", "b", "y"],
        label=sorted(event_log.keys()),
        linewidth=0,
    )
    plt.legend(loc=1, prop={"size": 6})
    plt.gca().set_xscale("log")
    plt.gca().set_yscale("log")
    plt.xlabel(r"$\mathrm{\tau}$ (s)")
    plt.ylabel("Freq (Arb. U.)")
    file_name = "".join([output_dir, "event_time_dist.pdf"])
    plt.savefig(file_name)
    print("Event time distribution saved as", file_name)


def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def plot_carrier_trajectories(all_carriers, parameter_dict, device_array, output_dir):
    combinations_to_plot = {"cathode": [], "anode": []}
    pop_list = []
    for index, carrier in enumerate(all_carriers):
        # First add all the excitons
        if "r_F" in carrier.__dict__:
            # We know this is an exciton
            combinations_to_plot[carrier.ID] = [carrier]
            pop_list.append(index)
    for index in sorted(pop_list, reverse=True):
        all_carriers.pop(index)
    for carrier in all_carriers:
        # Only electrons and holes remaining
        try:
            combinations_to_plot[carrier.injected_from].append(carrier)
        except KeyError:
            pass
    for inject_source, carriers in combinations_to_plot.items():
        plot3D_trajectory(
            inject_source, carriers, parameter_dict, device_array, output_dir
        )


def plot3D_trajectory(
    inject_source, carriers_to_plot, parameter_dict, device_array, output_dir
):
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    [x_len, y_len, z_len] = 1E10 * (
        np.array(device_array.shape) * parameter_dict["morphology_cell_size"]
    )
    # The conversion is needed to get the box in ang
    # Draw boxlines
    # Varying X
    ax.plot([0, x_len], [0, 0], [0, 0], c="k", linewidth=1.0)
    ax.plot([0, x_len], [y_len, y_len], [0, 0], c="k", linewidth=1.0)
    ax.plot([0, x_len], [0, 0], [z_len, z_len], c="k", linewidth=1.0)
    ax.plot([0, x_len], [y_len, y_len], [z_len, z_len], c="k", linewidth=1.0)
    # Varying Y
    ax.plot([0, 0], [0, y_len], [0, 0], c="k", linewidth=1.0)
    ax.plot([x_len, x_len], [0, y_len], [0, 0], c="k", linewidth=1.0)
    ax.plot([0, 0], [0, y_len], [z_len, z_len], c="k", linewidth=1.0)
    ax.plot([x_len, x_len], [0, y_len], [z_len, z_len], c="k", linewidth=1.0)
    # Varying Z
    ax.plot([0, 0], [0, 0], [0, z_len], c="k", linewidth=1.0)
    ax.plot([0, 0], [y_len, y_len], [0, z_len], c="k", linewidth=1.0)
    ax.plot([x_len, x_len], [0, 0], [0, z_len], c="k", linewidth=1.0)
    ax.plot([x_len, x_len], [y_len, y_len], [0, z_len], c="k", linewidth=1.0)

    if (inject_source is "anode") or (inject_source is "cathode"):
        carrier_string = inject_source
    else:
        # Exciton
        carrier_string = "exciton_{:05d}".format(carriers_to_plot[0].ID)
    for carrier in carriers_to_plot:
        if "r_F" in carrier.__dict__:
            color = "g"
        elif carrier.carrier_type == ELECTRON:
            color = "b"
        elif carrier.carrier_type == HOLE:
            color = "r"
        try:
            for hop_index, hop in enumerate(carrier.history[:-1]):
                # Note that conversion factors are needed as hop[0] * morphCellSize is
                # in m, and hop[1] is in ang.
                # Additionally, we need to add half of the morphCellSize to whatever
                # hop[0] is as the origin is in the centre of the box
                current_posn = (
                    1E10
                    * (
                        (np.array(hop[0]) + np.array([0.5, 0.5, 0.5]))
                        * parameter_dict["morphology_cell_size"]
                    )
                ) + np.array(hop[1])
                next_hop = carrier.history[hop_index + 1]
                next_posn = (
                    1E10
                    * (
                        (np.array(next_hop[0]) + np.array([0.5, 0.5, 0.5]))
                        * parameter_dict["morphology_cell_size"]
                    )
                ) + np.array(next_hop[1])
                ax.plot(
                    [current_posn[0], next_posn[0]],
                    [current_posn[1], next_posn[1]],
                    [current_posn[2], next_posn[2]],
                    c=color,
                    linewidth=0.5,
                )
        except AttributeError:
            hf.write_to_file(
                log_file,
                [
                    "Something has gone wrong while plotting. This carrier has no"
                    " history:",
                    repr(carrier.__dict__),
                    "Continuing...",
                ],
            )
            continue
    file_name = "".join([output_dir, carrier_string, "_traj.pdf"])
    plt.savefig(file_name)
    hf.write_to_file(log_file, ["Figure saved as", file_name])
    plt.close()


def debug_check_occupied(ID, image):
    print(
        "UPDATED =", globalChromophoreData.returnSpecificChromophore(image, ID).occupied
    )


def execute(
    device_array,
    chromophore_data,
    morphology_data,
    parameter_dict,
    voltage_val,
    time_limit,
):
    # ---=== PROGRAMMER'S NOTE ===---
    # This `High Resolution' version of the code will permit chromophore-based
    # hopping through the device, rather than just approximating the
    # distribution. We have the proper resolution there, let's just use it!
    # ---=========================---
    global global_chromophore_data
    global global_morphology_data
    global global_time
    global current_field_value
    global KMC_iterations
    global fastest_event_allowed
    global slowest_event_allowed
    global log_file

    global_chromophore_data = chromophore_data
    global_morphology_data = morphology_data
    if parameter_dict["fastest_event_allowed"] is not None:
        fastest_event_allowed = parameter_dict["fastest_event_allowed"]
    if parameter_dict["slowest_event_allowed"] is not None:
        slowest_event_allowed = parameter_dict["slowest_event_allowed"]
    # Given a voltage, the field value corresponding to it = (((bandgap
    # - el_inj_barrier - ho_inj_barrier) - Voltage) / z-extent)
    current_field_value = (
        # Bandgap:
        (
            (parameter_dict["acceptor_LUMO"] - parameter_dict["donor_HOMO"])
            -
            # Electron Inject Barrier:
            (parameter_dict["acceptor_LUMO"] - parameter_dict["cathode_work_function"])
            -
            # Hole Inject Barrier:
            (parameter_dict["anode_work_function"] - parameter_dict["donor_HOMO"])
            -
            # Voltage
            voltage_val
        )
        /
        # Z-extent:
        (device_array.shape[2] * parameter_dict["morphology_cell_size"])
    )
    hf.write_to_file(
        log_file, ["Current E-field value = {:s} Vm^{-1}".format(current_field_value)]
    )
    output_figures_dir = os.path.join(
        parameter_dict["output_device_dir"],
        parameter_dict["device_morphology"],
        "figures",
        str(voltage_val),
    )

    # DEBUG
    slowest_event = 0
    fastest_event = 2E99
    event_times = []

    # Need to initalise a bunch of variables
    event_queue = []
    photoinjection_rate = calculate_photoinjection_rate(
        parameter_dict, device_array.shape
    )
    number_of_photoinjections = 0
    number_of_cathode_injections = 0
    number_of_anode_injections = 0
    exciton_index = 0
    carrier_index = 0
    output_current_converged = False
    recombining_carrier_IDs = []
    if parameter_dict["record_carrier_history"] is True:
        all_carriers = []

    # Calculate the convergence characteristics
    check_conv_every = (
        int(parameter_dict["minimum_number_of_photoinjections"] / 100) * 5
    )
    previous_check = number_of_photoinjections
    conv_global_time = []
    conv_extractions = []
    diss_exciton_disp = []
    diss_exciton_time = []
    rec_exciton_disp = []
    rec_exciton_time = []
    number_of_hops = []
    number_of_dissociations = 0
    number_of_recombinations = 0
    number_of_electrical_recombinations = 0  # The total number of charges that
    # have recombined with an injected carrier from the `correct' contact
    event_log = {
        "photo": [],
        "cathode-injection": [],
        "anode-injection": [],
        "exciton_hop": [],
        "carrier_hop": [],
        "recombine": [],
    }
    already_printed = False
    output_print_statement = False

    # As the morphology is not changing, we can calculate the dark inject rates
    # and locations at the beginning and then not have to do it again, just
    # re-queue up the injectSites as we use them by calling
    # site(calculateInjectTime)
    injection_data = calculate_dark_current_injections(device_array, parameter_dict)
    cathode_inject_rate = injection_data[0]
    anode_inject_rate = injection_data[1]
    cathode_inject_queue = injection_data[2]
    anode_inject_queue = injection_data[3]

    print("\n\n ---=== MAIN KMC LOOP START ===---\n\n")
    t0 = T.time()
    # Main KMC loop
    # Put this in a try-except to plot on keyboardInterrupt
    try:
        while True:
            KMC_iterations += 1
            t1 = T.time()
            # TODO Check whether the output current has converged, after we have
            # had sufficient photoinjections. To check convergence, create a
            # datapoint every time 5% of the parameter-file-determined
            # photoinjections have been completed. This datapoint is the
            # number_of_extractions as a function of time. When the required
            # number of photoinjections has been completed, start performing
            # linear regression on the dataset and set output_current_converged
            # when the number_of_extractions saturates.
            if number_of_photoinjections >= (previous_check + check_conv_every):
                # Regardless of how many photoinjections we have, add to the
                # convergence datasets every 5% of the minimum required
                conv_global_time.append(global_time)
                conv_extractions.append(number_of_extractions)
                if (
                    number_of_photoinjections
                    > parameter_dict["minimum_number_of_photoinjections"]
                ):
                    # If we've gone over the minimum, then perform the
                    # convergence check to see if the output current has
                    # converged.
                    # TODO Check convergence, but for now, let's just see what
                    # the curve looks like
                    break
            # Break if current has converged
            if output_current_converged is True:
                break
            # Break if less than an hour before SLURM will kill the job
            if (time_limit is not None) and (t1 > t0 + time_limit - 3600):
                hf.write_to_file(
                    log_file,
                    ["LESS THAN ONE HOUR TO GO, TERMINATING JOB AT THIS POINT"],
                )
                break
            if len(event_queue) == 0:
                # Either the simulation has just started, or everything just
                # recombined/left the device. Therefore, queue up one
                # photoinjection and one dark current injection to get us
                # rolling. These are the only available kinetic starting points
                # for the simulation. Everything else stems from these.
                # Photoinjection
                photoinjection_time = hf.determine_event_tau(
                    photoinjection_rate,
                    event_type="photoinjection",
                    slowest_event_allowed=slowest_event_allowed,
                    fastest_event=fastest_event_allowed,
                    maximum_attempts=100,
                    log_file=log_file,
                )
                event_queue = push_to_queue(
                    event_queue, (photoinjection_time, "photo", None)
                )
                number_of_photoinjections += 1
                if parameter_dict["disable_dark_injection"] is False:
                    # Dark Injection:
                    # We now split up dark injection to be separate between the
                    # anode and the cathode, to ensure detail balance. We will
                    # use the cathode_inject_rate (the mean of the distribution)
                    # and the anode_inject_rate to determine when to perform the
                    # injection. Then, we will pop the first event of the
                    # cathode_inject_wait_times and anode_inject_wait_times as
                    # required, executing it IMMEDIATELY and decrementing the
                    # time in the rest of the queue. Then, we can re-queue the
                    # inject site. This means that the wait times in
                    # cathode_inject_wait_times and anode_inject_wait_times are
                    # actually only used to determine the order of the
                    # chromophores onto which we inject from each electrode, and
                    # do not increment the global_time.
                    # The global time is incremented according to
                    # cathode_injection_time and anode_injection_time instead.
                    # First, deal with the cathode injection
                    cathode_injection_time = hf.determine_event_tau(
                        cathode_inject_rate,
                        event_type="cathode-injection",
                        slowest_event_allowed=slowest_event_allowed,
                        fastest_event=fastest_event_allowed,
                        maximum_attempts=100,
                        log_file=log_file,
                    )
                    # Sort out the cathodeInjectQueue by getting the next inject
                    # site, and requeueing it with a new time to update the
                    # queue.
                    next_cathode_event, cathode_inject_queue = get_next_dark_event(
                        cathode_inject_queue, "cathode"
                    )
                    # Now push the cathode event to the main queue
                    event_queue = push_to_queue(
                        event_queue,
                        (
                            cathode_injection_time,
                            "cathode-injection",
                            next_cathode_event[2],
                        ),
                    )

                    # Now, deal with the anode injection
                    anode_injection_time = hf.determine_event_tau(
                        anode_inject_rate,
                        event_type="anode-injection",
                        slowest_event_allowed=slowest_event_allowed,
                        fastest_event=fastest_event_allowed,
                        maximum_attempts=100,
                        log_file=log_file,
                    )
                    # Sort out the anodeInjectQueue by popping the next inject
                    # site, decrementing the queue and re-queueing the inject
                    # site
                    next_anode_event, anode_inject_queue = get_next_dark_event(
                        anode_inject_queue, "anode"
                    )
                    # Now push the anode event to the main queue
                    event_queue = push_to_queue(
                        event_queue,
                        (anode_injection_time, "anode-injection", next_anode_event[2]),
                    )

            if int(t1 - t0) % 10 == 0:
                output_print_statement = True
            else:
                output_print_statement = False
                already_printed = False
            if output_print_statement and not already_printed:
                hf.write_to_file(
                    log_file,
                    [
                        "Current runtime = {0:d} s, with {1:d} events currently in the queue and {2:d} carriers currently in the system. Currently completed {3:d} iterations and simulated {4:.2e} s.".format(
                            t1 - t0,
                            len(event_queue),
                            len(global_carrier_dict.keys()),
                            KMC_iterations,
                            global_time,
                        )
                    ],
                )
                already_printed = True

            # Now find out what the next event is
            next_event = heapq.heappop(event_queue)
            # Increment the global time and decrement all of the other times in
            # the queue
            event_queue = decrement_time(event_queue, next_event[0])

            # DEBUG
            if next_event[0] > slowest_event:
                slowest_event = next_event[0]
            if next_event[0] < fastest_event:
                fastest_event = next_event[0]
            event_times.append(next_event[0])
            event_log[next_event[1]].append(next_event[0])

            # Execute the next behaviour (being sure to recalculate rates for
            # any new particles that have appeared)
            if next_event[1].lower() == "photo":
                # Complete the event by injecting an exciton
                # First find an injection location. For now, this will just be
                # somewhere random in the system
                random_device_position = [
                    np.random.randint(0, x - 1) for x in device_array.shape
                ]
                hf.write_to_file(
                    log_file,
                    [
                        "".join(
                            [
                                "EVENT: Photoinjection #{:d} into ".format(
                                    number_of_photoinjections
                                ),
                                repr(random_device_position),
                                " (which has type ",
                                repr(device_array[tuple(random_device_position)]),
                                ") after {0:d} iterations (global_time = {1:.2e})".format(
                                    KMC_iterations, global_time
                                ),
                            ]
                        )
                    ],
                )
                injected_exciton = exciton(
                    exciton_index, global_time, random_device_position, parameter_dict
                )
                if (
                    (injected_exciton.can_dissociate is True)
                    or (injected_exciton.hop_time is None)
                    or (injected_exciton.destination_chromophore is None)
                ):
                    # Injected onto either a dissociation site or a trap site
                    if injected_exciton.can_dissociate is True:
                        hf.write_to_file(
                            log_file, ["event: exciton dissociating immediately"]
                        )
                        number_of_dissociations += 1
                        number_of_hops.append(injected_exciton.number_of_hops)
                        # Create the carrier instances, but don't yet calculate
                        # their behaviour (need to add them to the carrier list
                        # before we can calculate the energetics)
                        # Also add the carriers to the carrier dictionary for
                        # when we need to calc delta_E in the device
                        injected_electron = carrier(
                            carrier_index,
                            global_time,
                            injected_exciton.current_device_posn,
                            injected_exciton.electron_chromophore,
                            injected_exciton.ID,
                            parameter_dict,
                        )
                        global_carrier_dict[carrier_index] = injected_electron
                        carrier_index += 1
                        injected_hole = carrier(
                            carrier_index,
                            global_time,
                            injected_exciton.current_device_posn,
                            injected_exciton.hole_chromophore,
                            injected_exciton.ID,
                            parameter_dict,
                        )
                        global_carrier_dict[carrier_index] = injected_hole
                        carrier_index += 1
                        # Update the injected chromophores to be marked as
                        # occupied
                        global_chromophore_data.return_specific_chromophore(
                            injected_exciton.current_device_posn,
                            injected_exciton.electron_chromophore.ID,
                        ).occupied.append(injected_exciton.current_device_posn)
                        global_chromophore_data.return_specific_chromophore(
                            injected_exciton.current_device_posn,
                            injected_exciton.hole_chromophore.ID,
                        ).occupied.append(injected_exciton.current_device_posn)
                        hf.write_to_file(
                            log_file,
                            [
                                " ".join(
                                    [
                                        "Exciton depositing electron at",
                                        repr(injected_exciton.current_device_posn),
                                        repr(
                                            injected_exciton.electron_chromophore.posn
                                        ),
                                    ]
                                ),
                                " ".join(
                                    [
                                        "Exciton depositing hole at",
                                        repr(injected_exciton.current_device_posn),
                                        repr(injected_exciton.hole_chromophore.posn),
                                    ]
                                ),
                            ],
                        )
                        # Now determine the behaviour of both carriers, and add
                        # their next hops to the KMC queue
                        injected_electron.calculate_behaviour()
                        if (injected_electron.hop_time is not None) and (
                            injected_electron.destination_chromophore is not None
                        ):
                            event_queue = push_to_queue(
                                event_queue,
                                (
                                    injected_electron.hop_time,
                                    "carrier_hop",
                                    injected_electron,
                                ),
                            )
                        injected_hole.calculate_behaviour()
                        if (injected_hole.hop_time is not None) and (
                            injected_hole.destination_chromophore is not None
                        ):
                            event_queue = push_to_queue(
                                event_queue,
                                (injected_hole.hop_time, "carrier_hop", injected_hole),
                            )
                    else:
                        # Injected onto a site with no connections, so this
                        # exciton will eventually die
                        hf.write_to_file(
                            log_file, ["EVENT: Exciton recombining immediately"]
                        )
                        number_of_recombinations += 1
                        number_of_hops.append(injected_exciton.number_of_hops)
                else:
                    # Hopping permitted, push the exciton to the queue.
                    # Push the exciton to the queue
                    event_queue = push_to_queue(
                        event_queue,
                        (injected_exciton.hop_time, "exciton_hop", injected_exciton),
                    )
                # A photoinjection has just occured, so now queue up a new one
                photoinjection_time = hf.determine_event_tau(
                    photoinjection_rate,
                    event_type="photoinjection",
                    slowest_event_allowed=slowest_event_allowed,
                    fastest_event=fastest_event_allowed,
                    maximum_attempts=100,
                    log_file=log_file,
                )
                event_queue = push_to_queue(
                    event_queue, (photoinjection_time, "photo", None)
                )
                # Increment the exciton and photoinjection counters
                exciton_index += 1
                number_of_photoinjections += 1

            elif (next_event[1].lower() == "cathode-injection") or (
                next_event[1].lower() == "anode-injection"
            ):
                inject_site = next_event[2]
                if inject_site.electrode.lower() == "cathode":
                    number_of_injections = number_of_cathode_injections
                else:
                    number_of_injections = number_of_anode_injections
                inject_chromophore = global_chromophore_data.return_specific_chromophore(
                    inject_site.device_posn, inject_site.chromophore.ID
                )
                if inject_site.device_posn not in inject_chromophore.occupied:
                    hf.write_to_file(
                        log_file,
                        [
                            "".join(
                                [
                                    "EVENT: Dark current injection from the {0:s} #{1:d} into ".format(
                                        inject_site.electrode, number_of_injections
                                    ),
                                    repr(inject_site.device_posn),
                                    " (which has type ",
                                    repr(device_array[tuple(inject_site.device_posn)]),
                                    ") chromophore number {0:d} after {1:d} iterations (global_time = {2:.2e})",
                                    format(
                                        inject_site.chromophore.ID,
                                        KMC_iterations,
                                        global_time,
                                    ),
                                ]
                            )
                        ],
                    )
                    # Inject the carrier
                    injected_carrier = carrier(
                        carrier_index,
                        global_time,
                        inject_site.device_posn,
                        inject_site.chromophore,
                        inject_site.electrode,
                        parameter_dict,
                        injected_onto_site=inject_site,
                    )
                    global_carrier_dict[carrier_index] = injected_carrier
                    carrier_index += 1
                    # Update the chromophore occupation
                    inject_chromophore.occupied.append(inject_site.device_posn)
                    # Determine the injected carrier's next hop and queue it
                    injected_carrier.calculate_behaviour()
                    if (injected_carrier.hop_time is not None) and (
                        injected_carrier.hop_time > 1
                    ):
                        hf.write_to_file(
                            log_file,
                            ["DARK INJECTION LED TO ELECTRON WITH CRAZY HOPTIME"],
                        )
                        for carrier_from_list in global_carrier_dict.values():
                            hf.write_to_file(
                                log_file,
                                [
                                    " ".join(
                                        [
                                            repr(carrier_from_list.current_device_posn),
                                            repr(
                                                carrier_from_list.current_chromophore.posn
                                            ),
                                        ]
                                    )
                                ],
                            )
                        hf.write_to_file(
                            log_file, [str(injected_carrier.current_chromophore.ID)]
                        )
                        hf.write_to_file(log_file, [repr(injected_carrier.__dict__)])
                        exit()
                    if (injected_carrier.hop_time is not None) and (
                        injected_carrier.destination_chromophore is not None
                    ):
                        event_queue = push_to_queue(
                            event_queue,
                            (
                                injected_carrier.hop_time,
                                "carrier_hop",
                                injected_carrier,
                            ),
                        )
                # Now determine the next DC event and queue it
                if inject_site.electrode.lower() == "cathode":
                    next_cathode_event, cathode_inject_queue = get_next_dark_event(
                        cathode_inject_queue, "cathode"
                    )
                    event_queue = push_to_queue(
                        event_queue,
                        (
                            cathode_injection_time,
                            "cathode-injection",
                            next_cathode_event[2],
                        ),
                    )
                    number_of_cathode_injections += 1
                if inject_site.electrode.lower() == "anode":
                    next_anode_event, anode_inject_queue = get_next_dark_event(
                        anode_inject_queue, "anode"
                    )
                    event_queue = push_to_queue(
                        event_queue,
                        (anode_injection_time, "anode-injection", next_anode_event[2]),
                    )
                    number_of_anode_injections += 1

            elif next_event[1].lower() == "exciton_hop":
                hopping_exciton = next_event[2]
                # There is a sporadic (rare) bug that causes excitons to
                # sometimes get queued up to hop  even though they have already
                # recombined. The following is similar to the check we do before
                # the carrier hop, and should fix that.
                if hopping_exciton.removed_time is not None:
                    # Exciton has already been removed from the system through
                    # recombination
                    continue
                hopping_exciton.perform_hop()
                if hopping_exciton.removed_time is None:
                    # As long as this exciton hasn't just been removed,
                    # recalculate its behaviour
                    hopping_exciton.calculate_behaviour()
                # At this point, dissociate the exciton or remove it from the
                # system if
                # needed.
                if (
                    (hopping_exciton.can_dissociate is True)
                    or (hopping_exciton.hop_time is None)
                    or (hopping_exciton.destination_chromophore is None)
                ):
                    # Exciton needs to be removed. As we've already popped it
                    # from the queue, we just need to not queue it up again.
                    if parameter_dict["record_carrier_history"] is True:
                        all_carriers.append(hopping_exciton)

                    if hopping_exciton.can_dissociate is True:
                        hf.write_to_file(
                            log_file,
                            [
                                "EVENT: Exciton dissociating after {0:d} iterations (global_time = {1:.2e})".format(
                                    KMC_iterations, global_time
                                )
                            ],
                        )
                        number_of_dissociations += 1
                        number_of_hops.append(injected_exciton.number_of_hops)
                        # Create the carrier instances, but don't yet calculate
                        # their behaviour (need to add them to the carrier list
                        # before we can calculate the energetics)
                        # Also add the carriers to the carrier dictionary for
                        # when we need to calc delta_E in the device. Plop an
                        # electron down on the electron_chromophore of the
                        # dissociating exciton, and a hole on the
                        # hole_chromophore
                        injected_electron = carrier(
                            carrier_index,
                            global_time,
                            hopping_exciton.current_device_posn,
                            hopping_exciton.electron_chromophore,
                            hopping_exciton.ID,
                            parameter_dict,
                        )
                        global_carrier_dict[carrier_index] = injected_electron
                        carrier_index += 1
                        injected_hole = carrier(
                            carrier_index,
                            global_time,
                            hopping_exciton.current_device_posn,
                            hopping_exciton.hole_chromophore,
                            hopping_exciton.ID,
                            parameter_dict,
                        )
                        global_carrier_dict[carrier_index] = injected_hole
                        carrier_index += 1
                        # Update the injected chromophores to be marked as
                        # occupied
                        global_chromophore_data.return_specific_chromophore(
                            hopping_exciton.current_device_posn,
                            hopping_exciton.electron_chromophore.ID,
                        ).occupied.append(hopping_exciton.current_device_posn)
                        global_chromophore_data.return_specific_chromophore(
                            hopping_exciton.current_device_posn,
                            hopping_exciton.hole_chromophore.ID,
                        ).occupied.append(hopping_exciton.current_device_posn)
                        # Add to the allCarriers list for plotting
                        if parameter_dict["record_carrier_history"] is True:
                            all_carriers += [injected_electron, injected_hole]
                        # Now determine the behaviour of both carriers, and add
                        # their next hops to the KMC queue
                        injected_electron.calculate_behaviour()
                        if (injected_electron.hop_time is not None) and (
                            injected_electron.destination_chromophore is not None
                        ):
                            event_queue = push_to_queue(
                                event_queue,
                                (
                                    injected_electron.hop_time,
                                    "carrier_hop",
                                    injected_electron,
                                ),
                            )
                        injected_hole.calculate_behaviour()
                        if (injected_hole.hop_time is not None) and (
                            injected_hole.destination_chromophore is not None
                        ):
                            event_queue = push_to_queue(
                                event_queue,
                                (injected_hole.hop_time, "carrier_hop", injected_hole),
                            )
                    else:
                        hf.write_to_file(
                            log_file,
                            [
                                "EVENT: Exciton #{0:d} recombining after {1:d} iterations (global_time = {2:.2e})".format(
                                    hopping_exciton.ID, KMC_iterations, global_time
                                )
                            ],
                        )
                        number_of_recombinations += 1
                        number_of_hops.append(injected_exciton.number_of_hops)
                    # DEBUG
                    # Calculate the initial position and final positions and
                    # append the excitonDisp with the separation
                    initial_pos = np.array(
                        hopping_exciton.initial_device_posn
                    ) * parameter_dict["morphology_cell_size"] + (
                        np.array(hopping_exciton.initial_chromophore.posn) * 1E-10
                    )
                    final_pos = np.array(
                        hopping_exciton.current_device_posn
                    ) * parameter_dict["morphology_cell_size"] + (
                        np.array(hopping_exciton.current_chromophore.posn) * 1E-10
                    )
                    if hopping_exciton.can_dissociate is True:
                        diss_exciton_disp.append(
                            hf.calculate_separation(initial_pos, final_pos) / 1E-9
                        )
                        diss_exciton_time.append(hopping_exciton.recombination_time)
                    else:
                        rec_exciton_disp.append(
                            hf.calculate_separation(initial_pos, final_pos) / 1E-9
                        )
                        rec_exciton_time.append(hopping_exciton.recombination_time)
                    # END DEBUG
                else:
                    event_queue = push_to_queue(
                        event_queue,
                        (hopping_exciton.hop_time, "exciton_hop", injected_exciton),
                    )

            elif next_event[1].lower() == "carrier_hop":
                hopping_carrier = next_event[2]
                # Check that the carrier is still in the carrier dictionary. If
                # it's not, then it has already been removed from the system and
                # we can safely ignore it
                if hopping_carrier.ID not in global_carrier_dict.keys():
                    continue
                hopping_carrier.perform_hop()
                if hopping_carrier.removed_time is None:
                    # As long as this carrier hasn't just been removed,
                    # recalculate its behaviour
                    hopping_carrier.calculate_behaviour()
                if (hopping_carrier.hop_time is None) or (
                    hopping_carrier.destination_chromophore is None
                ):
                    # Carrier is either trapped with no eligible hops or has
                    # just been extracted
                    if hopping_carrier.removed_time is not None:
                        # Carrier has been extracted, so remove it from the
                        # carriers dictionary
                        global_carrier_dict.pop(hopping_carrier.ID)
                    # Else: Carrier is trapped, so we'll just leave it there so
                    # it affects the Coulombic landscape
                    if parameter_dict["record_carrier_history"] is True:
                        all_carriers.append(hopping_carrier)
                else:
                    # Normal carrier hop, so requeue this carrier
                    event_queue = push_to_queue(
                        event_queue,
                        (hopping_carrier.hop_time, "carrier_hop", hopping_carrier),
                    )
                # Check if we're eligible to recombine
                if (
                    (hopping_carrier.recombining is True)
                    and (hopping_carrier.ID not in recombining_carrier_IDs)
                    and (
                        hopping_carrier.recombining_with not in recombining_carrier_IDs
                    )
                ):
                    recombining_carrier_IDs.append(hopping_carrier.ID)
                    recombining_carrier_IDs.append(hopping_carrier.recombining_with)
                    recombination_time = hf.determine_event_tau(
                        parameter_dict["recombination_rate"],
                        event_type="carrier-recombination",
                        slowest_event_allowed=slowest_event_allowed,
                        fastest_event=fastest_event_allowed,
                        maximum_attempts=100,
                        log_file=log_file,
                    )
                    event_queue = push_to_queue(
                        event_queue, (recombination_time, "recombine", hopping_carrier)
                    )
                # # If this carrier was injected, requeue this injection site in
                # # the dark_inject_queue to allow another injection. Only
                # # re-queueing up the dark inject site after the carrier has
                # # hopped prevents double injection onto the same chromophore,
                # # by waiting until the carrier hops away first.
                # # TODO, what if the carrier hops to a chromophore that was
                # # about to be injected into? We don't check for this, so could
                # # still get double occupancy!
                # if hopping_carrier.injected_onto_site is not None:
                #     new_dark_inject = copy.deepcopy(inject_site)
                #     new_dark_inject.calculate_inject_time()
                #     heapq.heappush(dark_inject_queue, (new_dark_inject.inject_time,
                #                                        'dark', new_dark_inject))

            elif next_event[1].lower() == "recombine":
                hf.write_to_file(
                    log_file,
                    [
                        "EVENT: Carrier recombination check after {0:d} iterations (global_time = {1:.2e})".format(
                            KMC_iterations, global_time
                        )
                    ],
                )
                # A recombination event is about to occur. At this point, we
                # should check if the carrier and its recombination partner are
                # still in range.
                delta_electrical_recombinations = 0
                try:
                    carrier1 = next_event[2]
                    carrier1posn = np.array(
                        carrier1.current_device_posn
                    ) * parameter_dict["morphology_cell_size"] + (
                        np.array(carrier1.current_chromophore.posn) * 1E-10
                    )
                    carrier2 = global_carrier_dict[carrier1.recombining_with]
                    carrier2posn = np.array(
                        carrier2.current_device_posn
                    ) * parameter_dict["morphology_cell_size"] + (
                        np.array(carrier2.current_chromophore.posn) * 1E-10
                    )
                    separation = hf.calculate_separation(carrier2posn, carrier1posn)
                    recombining_carrier_IDs.remove(carrier2.ID)
                    # Calculate the increment to the
                    # number_of_electrical_recombinations counter (this will be
                    # used to calculate the current density later)
                    # Note that optical recombinations do not contribute to
                    # photocurrent
                    for recombining_carrier in [carrier1, carrier2]:
                        if (
                            (recombining_carrier.carrier_type == HOLE)
                            and (recombining_carrier.injected_from.lower() == "anode")
                        ) or (
                            (recombining_carrier.carrier_type == ELECTRON)
                            and (recombining_carrier.injected_from.lower() == "cathode")
                        ):
                            delta_electrical_recombinations += 1
                        elif (
                            (recombining_carrier.carrier_type == ELECTRON)
                            and (recombining_carrier.injected_from.lower() == "anode")
                        ) or (
                            (recombining_carrier.carrier_type == HOLE)
                            and (recombining_carrier.injected_from.lower() == "cathode")
                        ):
                            delta_electrical_recombinations -= 1
                except (ValueError, KeyError):
                    # The second carrier is missing from the simulation (already
                    # extracted), so set the separation to be large
                    separation = 1E99
                recombining_carrier_IDs.remove(carrier1.ID)
                if separation <= parameter_dict["coulomb_capture_radius"]:
                    hf.write_to_file(
                        log_file,
                        [
                            "{0:f} <= {1:f}".format(
                                separation, parameter_dict["coulomb_capture_radius"]
                            )
                        ],
                    )
                    hf.write_to_file(
                        log_file,
                        [
                            "EVENT: Carrier recombination succeeded after {0:d} iterations (global_time = {1:.2e})".format(
                                KMC_iterations, global_time
                            )
                        ],
                    )
                    # Carriers are in range, so recombine them
                    carrier1.removed_time = global_time
                    carrier2.removed_time = global_time
                    number_of_recombinations += 1
                    number_of_electrical_recombinations += (
                        delta_electrical_recombinations
                    )
                    try:
                        global_carrier_dict.pop(carrier1.ID)
                    except KeyError:
                        # Carrier has already been removed (extracted at contacts while
                        # waiting for recombination)
                        pass
                    try:
                        global_carrier_dict.pop(carrier2.ID)
                    except KeyError:
                        # Carrier has already been removed (extracted at contacts while
                        # waiting for recombination)
                        pass
                    if parameter_dict["record_carrier_history"] is True:
                        all_carriers += [carrier1, carrier2]
                else:
                    hf.write_to_file(
                        log_file,
                        [
                            "{0:f} > {1:f}".format(
                                separation, parameter_dict["coulomb_capture_radius"]
                            )
                        ],
                    )
                    hf.write_to_file(
                        log_file,
                        [
                            "EVENT: Carrier recombination failed after {0:d} iterations (global_time = {1:.2e})".format(
                                KMC_iterations, global_time
                            )
                        ],
                    )
                    # Carriers are no longer in range, so the recombination fails.
                    # Update their recombination flags
                    carrier1.recombining = False
                    carrier1.recombining_with = None
                    if separation != 1E99:
                        # If recombining carrier was already extracted, then skip this
                        carrier2.recombining = False
                        carrier2.recombining_with = None
            else:
                print(event_queue)
                raise SystemError("Next event in queue is unknown")

    except KeyboardInterrupt:
        time = T.time() - t0
        hf.write_to_file(
            log_file,
            [
                "Kill command recieved...",
                "Plotting output graphs before terminating...",
                " ".join(["---=== Results from CPU rank", sys.argv[2], "===---"]),
                "Run terminated after {0:d} iterations (global_time = {1:.2e}) after {2:.2f} seconds.".format(
                    KMC_iterations, global_time, time
                ),
                "Number of photoinjections = {:d}".format(number_of_photoinjections),
                "Number of cathode injections = {:d}".format(
                    number_of_cathode_injections
                ),
                "Number of anode injections = {:d}".format(number_of_anode_injections),
                "Number of dissociations = {:d}".format(number_of_dissociations),
                "Number of recombinations = {:d}".format(number_of_recombinations),
                "Number of extractions = {:d}".format(number_of_extractions),
            ],
        )
        # print("DURING THIS RUN:")
        # print("Slowest Event Considered =", slowestEvent)
        # print("Fastest Event Considered =", fastestEvent)
        plot_event_time_distribution(
            event_log, output_figures_dir, fastest_event, slowest_event
        )
        if parameter_dict["record_carrier_history"] is True:
            plot_carrier_Z_profiles(
                all_carriers, parameter_dict, device_array, output_figures_dir
            )
            plot_carrier_trajectories(
                all_carriers, parameter_dict, device_array, output_figures_dir
            )
        exit()
    # THE FOLLOWING CODE WAS USED TO GENERATE THE ANALYSIS DATA USED IN THE
    # MATTYSUMMARIES EXCITON DYNAMICS STUFF
    # print("Plotting Exciton Characteristics")
    # plt.figure()
    # plt.hist(dissExcitonDisp, bins=15)
    # plt.title('Displacement until Dissociation')
    # plt.xlabel('Displacement, nm')
    # plt.ylabel('Frequency, Arb. U.')
    # plt.savefig(outputFiguresDir + 'dissExcitonDisp.pdf')
    # plt.clf()

    # plt.hist(np.array(dissExcitonTime) * 1E9, bins=15
    # plt.title('Time until Dissociation')
    # plt.xlabel('Time, ns')
    # plt.ylabel('Frequency, Arb. U.')
    # plt.savefig(outputFiguresDir + 'dissExcitonTime.pdf')
    # plt.clf()

    # plt.hist(recExcitonDisp, bins=15)
    # plt.title('Displacement until Recombination')
    # plt.xlabel('Displacement, nm')
    # plt.ylabel('Frequency, Arb. U.')
    # plt.savefig(outputFiguresDir + 'recExcitonDisp.pdf')
    # plt.clf()

    # plt.hist(np.array(recExcitonTime) * 1E10, bins=15)
    # plt.title('Time until Recombination')
    # plt.xlabel('Time, ns')
    # plt.ylabel('Frequency, Arb. U.')
    # plt.savefig(outputFiguresDir + 'recExcitonTime.pdf')
    # plt.clf()

    # plt.hist(numberOfHops, bins=15)
    # plt.title('numberOfHops')
    # plt.savefig(outputFiguresDir + 'numberOfExcitonHops.pdf')
    # plt.clf()

    # print("XDE =", numberOfDissociations
    #       / parameterDict['minimumNumberOfPhotoinjections'])
    # print("Mean/SD DissDisp =", np.mean(dissExcitonDisp), "+-",
    #       np.std(dissExcitonDisp)/np.sqrt(len(dissExcitonDisp)))
    # print("Mean/SD DissTime =", np.mean(dissExcitonTime), "+-",
    #       np.std(dissExcitonTime)/np.sqrt(len(dissExcitonTime)))
    # print("Mean/SD RecDisp =", np.mean(recExcitonDisp), "+-",
    #       np.std(recExcitonDisp)/np.sqrt(len(recExcitonDisp)))
    # print("Mean/SD RecTime =", np.mean(recExcitonTime), "+-",
    #       np.std(recExcitonTime)/np.sqrt(len(recExcitonTime)))

    # plotConnections(excitonPath, chromophoreData.returnChromophoreList([1,8,6]),
    #                 morphologyData.returnAAMorphology([1,8,6]), outputFiguresDir)

    # print("Plotting the MSD of the exciton")
    # plt.figure()
    # MSD = list(np.array(excitonDisp)**2)
    # plt.plot(excitonTime, MSD, c='b')
    # fit = np.polyfit(excitonTime, MSD, 1)
    # print("Fit =", fit)
    # print("Diffusion Coefficient =", fit[0])
    # xFit = np.linspace(0, max(excitonTime), 100)
    # yFit = [(xVal*fit[0] + fit[1]) for xVal in xFit]
    # plt.plot(xFit, yFit, c='b')
    # plt.savefig(outputFiguresDir + 'excitonMSD.pdf')

    # print("Plotting the extraction data as a function of simulation time")
    # plt.figure()
    # plt.plot(convGlobalTime, convExtractions)
    # plt.savefig(outputFiguresDir + 'convergence.pdf')
    # return
    time = T.time() - t0
    hf.write_to_file(
        log_file,
        [
            "SIMULATION COMPLETED",
            "Plotting output graphs...",
            " ".join(["---=== Results from CPU rank", sys.argv[2], "===---"]),
            "Simulation voltage = {:.2f}".format(voltage_val),
            "Run completed after {0:d} iterations (global_time = {1:.2e}) after {2:.2f} seconds.".format(
                KMC_iterations, global_time, time
            ),
            "Number of photoinjections = {:d}".format(number_of_photoinjections),
            "Number of cathode injections = {:d}".format(number_of_cathode_injections),
            "Number of anode injections = {:d}".format(number_of_anode_injections),
            "Number of dissociations = {:d}".format(number_of_dissociations),
            "Number of recombinations = {:d}".format(number_of_recombinations),
            "Number of extractions = {:d}".format(number_of_extractions),
        ],
    )
    plot_event_time_distribution(
        event_log, output_figures_dir, fastest_event, slowest_event
    )
    if parameter_dict["record_carrier_history"] is True:
        plot_carrier_Z_profiles(
            all_carriers, parameter_dict, device_array, output_figures_dir
        )
        plot_carrier_trajectories(
            all_carriers, parameter_dict, device_array, output_figures_dir
        )


def slurm_time_in_S(slurm_time):
    # Expects slurmTime in HH:MM:SS
    split_time = slurm_time.split(":")
    time_in_S = (
        int(split_time[0]) * 60 ** 2 + int(split_time[1]) * 60 + int(split_time[2])
    )
    return time_in_S


def main():
    global log_file

    KMC_directory = sys.argv[1]
    CPU_rank = int(sys.argv[2])
    np.random.seed(int(sys.argv[3]))
    overwrite = False
    try:
        overwrite = bool(sys.argv[4])
    except:
        pass

    # Get the time limit
    time_limit = os.getenv("SLURM_TIMELIMIT", None)
    if time_limit is not None:
        time_limit = slurm_time_in_S(time_limit)

    jobs_file_name = os.path.join(
        KMC_directory, "KMC_data_{:02d}.pickle".format(CPU_rank)
    )
    device_data_file_name = KMC_directory.replace("/KMC", "/code/device_data.pickle")

    with open(device_data_file_name, "rb") as pickle_file:
        [device_array, chromophore_data, morphology_data, parameter_dict] = pickle.load(
            pickle_file
        )
    with open(jobs_file_name, "rb") as pickle_file:
        jobs_to_run = pickle.load(pickle_file)
    if parameter_dict["output_log_to_stdout"] is True:
        print("Redirecting log to standard out.")
        log_file = "stdout"
    else:
        log_file = os.path.join(KMC_directory, "KMC_log_{:02d}.log".format(CPU_rank))
        # Reset the log file
        with open(log_file, "wb+") as log_file_handle:
            pass
    hf.write_to_file(
        log_file,
        [
            "".join(
                ["Found {:d} jobs to run:".format(len(jobs_to_run)), repr(jobs_to_run)]
            )
        ],
    )
    hf.write_to_file(
        log_file, ["Using random number seed {:d}".format(int(sys.argv[3]))]
    )
    if parameter_dict["disable_coulombic"] is True:
        hf.write_to_file(log_file, ["COULOMBIC INTERACTIONS DISABLED"])
    if parameter_dict["disable_dark_injection"] is True:
        hf.write_to_file(log_file, ["DARK CURRENT INJECTION (FROM CONTACTS) DISABLED"])
    # hf.write_to_file(log_file, ['Found ' + str(len(jobs_to_run)) + ' jobs to run.'])
    # Set the affinities for this current process to make sure it's maximising
    # available CPU usage
    current_PID = os.getpid()
    # try:
    #     affinity_job = sp.Popen(['taskset', '-pc', str(CPU_rank), str(current_PID)],
    #                             stdin=sp.PIPE, stdout=sp.PIPE,
    #                             stderr=sp.PIPE).communicate()
    #     # hf.write_to_file(log_file, affinity_job[0].split('\n'))
    #     # hf.write_to_file(log_file, affinity_job[1].split('\n'))
    # except OSError:
    #     hf.writeToFile(log_file, ["Taskset command not found, skipping setting of"
    #                               " processor affinity..."])
    #     # hf.write_to_file(log_file, ["Taskset command not found, skipping setting"
    #                                   " of processor affinity..."])
    # Begin the simulation
    for voltage_val in jobs_to_run:
        execute(
            device_array,
            chromophore_data,
            morphology_data,
            parameter_dict,
            voltage_val,
            time_limit,
        )


if __name__ == "__main__":
    main()
