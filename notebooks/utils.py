# This file contains some utility functions for the notebooks

import h5py
from glob import glob
import pykep as pk
import numpy as np
from tqdm import tqdm
import pandas as pd


def load_particle_data_from_hdf5(filename, starting_t, dt, stepsize=1, verbose=False):
    """
    Loads an hdf5 file and returns the data.
    """
    print("Loading hdf5 from: ", filename)
    # Load hdf5 file
    data = h5py.File(filename)

    if verbose:
        print("Found keys", list(data.keys()))

    # Determine the iterations at which output was written
    iterations_idx_str = list(data["ParticleData"].keys())
    iterations_idx = []
    for it in iterations_idx_str:
        if it != "ConstantProperties":
            iterations_idx.append(int(it))
    iterations_idx.sort()
    max_iterations = max(iterations_idx)
    if verbose:
        print("Found a total of", max_iterations, " iterations.")

    # Find out the simulation runtime in days
    end_t = pk.epoch(starting_t.mjd2000 + max_iterations * dt * pk.SEC2DAY)
    total_days = end_t.mjd - starting_t.mjd
    if verbose:
        print("Simulation ran for a total of", total_days, " days.")

    iterations_idx = np.asarray(iterations_idx)
    if verbose:
        print("Total # of indices = ", len(iterations_idx))
    iterations_idx = iterations_idx[iterations_idx % 1000 == 0]
    if verbose:
        print("After cleaning # of indices = ", len(iterations_idx))

    rs, vs, ids = [], [], []  # will hold r,v, ids for whole simulation
    iterations_idx = iterations_idx[::stepsize]
    for idx in tqdm(iterations_idx, disable=not verbose):

        # Load velocities
        v_x = data["ParticleData"][str(idx)]["Particles"]["Velocities"][:]["x"]
        v_y = data["ParticleData"][str(idx)]["Particles"]["Velocities"][:]["y"]
        v_z = data["ParticleData"][str(idx)]["Particles"]["Velocities"][:]["z"]
        v = np.vstack([v_x, v_y, v_z]).transpose()

        # Load positions
        r_x = data["ParticleData"][str(idx)]["Particles"]["Positions"][:]["x"]
        r_y = data["ParticleData"][str(idx)]["Particles"]["Positions"][:]["y"]
        r_z = data["ParticleData"][str(idx)]["Particles"]["Positions"][:]["z"]
        r = np.vstack([r_x, r_y, r_z]).transpose()

        # Get IDs
        ID = np.array(data["ParticleData"][str(idx)]["Particles"]["IDs"])

        # Store in single large list
        rs.append(r * 1000.0)  # convert to m
        vs.append(v * 1000.0)  # convert to m
        ids.append(ID)

    # Some code breaks if constant properties change
    assert len(data["ParticleData"]["ConstantProperties"][0]) == 6
    constant_props = data["ParticleData"]["ConstantProperties"][()]

    return (
        iterations_idx,
        rs,
        vs,
        ids,
        constant_props,
        end_t,
        max_iterations,
    )


def load_conjunctions_and_evasions(filename, constant_props, verbose=False):
    """
    Loads the conjunction and evasion data from the hdf5 file and returns the data.
    """

    # Load hdf5 file
    data = h5py.File(filename)

    def get_particle_size(ID):
        for p in constant_props:
            if p[0] == ID:
                return p[3]

    def get_particle_status(ID):
        for p in constant_props:
            if p[0] == ID:
                return p[5]

    # Convert conjunctions to a pandas dataframe for convenience
    conj = pd.DataFrame(
        columns=[
            "P1",
            "P2",
            "Iteration",
            "SquaredDistance",
            "Size_P1[m]",
            "Size_P2[m]",
            "$\kappa$",
            "Status_P1",
            "Status_P2",
        ]
    )

    try:
        collision_keys = data["CollisionData"].keys()
    except Exception as e:
        collision_keys = None
        if verbose:
            print("No collisions found in ", filename)

    if collision_keys is not None:
        for idx, it in tqdm(
            enumerate(collision_keys), total=len(collision_keys), disable=not verbose
        ):
            iteration = int(it)
            collisions = data["CollisionData"][it]["Collisions"]
            for collision in collisions:
                size_1 = get_particle_size(collision[0])
                size_2 = get_particle_size(collision[1])
                status_1 = get_particle_status(collision[0])
                status_2 = get_particle_status(collision[1])
                conj = conj.append(
                    {
                        "P1": collision[0],
                        "P2": collision[1],
                        "Iteration": iteration,
                        "SquaredDistance": collision[2],
                        "Size_P1[m]": size_1,
                        "Size_P2[m]": size_2,
                        "$\kappa$": np.sqrt(collision[2])
                        / ((size_1 + size_2) / 1000.0),
                        "Status_P1": status_1,
                        "Status_P2": status_2,
                    },
                    ignore_index=True,
                )

    # Convert evasions to a pandas dataframe for convenience
    evasions_df = pd.DataFrame(
        columns=[
            "P1",
            "P2",
            "Iteration",
            "SquaredDistance",
            "Size_P1[m]",
            "Size_P2[m]",
            "$\kappa$",
            "Status_P1",
            "Status_P2",
        ]
    )

    try:
        evasion_keys = data["EvasionData"].keys()
    except Exception as e:
        if verbose:
            print("No evasions found in ", filename)
        evasion_keys = None

    if evasion_keys is not None:
        for idx, it in tqdm(
            enumerate(evasion_keys), total=len(evasion_keys), disable=not verbose
        ):
            iteration = int(it)
            evasions = data["EvasionData"][it]["Evasions"]
            for evasion in evasions:
                size_1 = get_particle_size(evasion[0])
                size_2 = get_particle_size(evasion[1])
                status_1 = get_particle_status(evasion[0])
                status_2 = get_particle_status(evasion[1])
                evasions_df = evasions_df.append(
                    {
                        "P1": evasion[0],
                        "P2": evasion[1],
                        "Iteration": iteration,
                        "SquaredDistance": evasion[2],
                        "Size_P1[m]": size_1,
                        "Size_P2[m]": size_2,
                        "$\kappa$": np.sqrt(evasion[2]) / ((size_1 + size_2) / 1000.0),
                        "Status_P1": status_1,
                        "Status_P2": status_2,
                    },
                    ignore_index=True,
                )

    return (
        conj,
        evasions_df,
    )


def load_hdf5s_from_mpi_run(path_to_results_folder, starting_t, dt, stepsize=1):
    """
    Loads all hdf5 files in a folder and returns the data.
    """
    print("Loading hdf5s from MPI run")

    # Get all hdf5 files in folder
    runs = glob(path_to_results_folder + "/*.h5_rank_*")

    # Init vars to hold the data
    rs, vs, ids = [], [], []
    end_t, max_iterations, iterations_idx = None, None, None
    constant_props = None

    conj = pd.DataFrame(
        columns=[
            "P1",
            "P2",
            "Iteration",
            "SquaredDistance",
            "Size_P1[m]",
            "Size_P2[m]",
            "$\kappa$",
            "Status_P1",
            "Status_P2",
        ]
    )

    evasions = pd.DataFrame(
        columns=[
            "P1",
            "P2",
            "Iteration",
            "SquaredDistance",
            "Size_P1[m]",
            "Size_P2[m]",
            "$\kappa$",
            "Status_P1",
            "Status_P2",
        ]
    )

    # Loop over all hdf5 files and load data,
    # we need two loops because we need some of the data for
    # the rest of it
    for run in runs:
        (
            iterations_idx_new,
            rs_i,
            vs_i,
            ids_i,
            constant_props_i,
            end_t_new,
            max_iterations_new,
        ) = load_particle_data_from_hdf5(run, starting_t, dt, stepsize)

        # Check data among ranks is consistent
        if end_t is None:
            iterations_idx = iterations_idx_new
            end_t = end_t_new
            max_iterations = max_iterations_new
        else:
            assert end_t.mjd == end_t_new.mjd
            assert max_iterations == max_iterations_new
            assert (iterations_idx == iterations_idx_new).all()

        if constant_props is None:
            constant_props = constant_props_i
            rs = rs_i
            vs = vs_i
            ids = ids_i
        else:
            constant_props = np.concatenate((constant_props, constant_props_i), axis=0)
            for it in range(len(iterations_idx)):
                rs[it] = np.concatenate((rs[it], rs_i[it]), axis=0)
                vs[it] = np.concatenate((vs[it], vs_i[it]), axis=0)
                ids[it] = np.concatenate((ids[it], ids_i[it]), axis=0)

    for run in runs:
        (conj_i, evasions_i,) = load_conjunctions_and_evasions(run, constant_props)

        conj = pd.concat([conj, conj_i])
        evasions = pd.concat([evasions, evasions_i])

    print("Found a total of", max_iterations, " iterations.")
    total_days = end_t.mjd - starting_t.mjd
    print("Simulation ran for a total of", total_days, " days.")
    print("Total # of indices (available timesteps) = ", len(iterations_idx))
    print("Total # of particles = ", len(ids[0]))
    print("Total # of evasions = ", len(evasions))
    print("Total # of conjunctions = ", len(conj))
    return (
        iterations_idx,
        rs,
        vs,
        ids,
        conj,
        evasions,
        constant_props,
        end_t,
        max_iterations,
    )

