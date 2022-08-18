# This file contains some utility functions for the notebooks

import h5py
from glob import glob
import pykep as pk
import numpy as np
from tqdm import tqdm
import pandas as pd


def load_hdf5(filename, starting_t, dt, stepsize=1, verbose=False):
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
    except:
        if verbose:
            print("No collisions found in ", filename)

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
    except Exception as e:
        if verbose:
            print("No evasions found in ", filename)

    # Some code breaks if constant properties change
    assert len(data["ParticleData"]["ConstantProperties"][0]) == 6
    constant_props = data["ParticleData"]["ConstantProperties"][()]

    return (
        iterations_idx,
        rs,
        vs,
        ids,
        conj,
        evasions_df,
        constant_props,
        end_t,
        max_iterations,
    )


def load_hdf5s_from_mpi_run(path_to_results_folder, starting_t, dt, stepsize=1):
    """
    Loads all hdf5 files in a folder and returns the data.
    """
    print("Loading hdf5s from MPI run")
    runs = glob(path_to_results_folder + "/*.h5_rank_*")
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

    for run in runs:
        (
            iterations_idx_new,
            rs_i,
            vs_i,
            ids_i,
            conj_i,
            evasions_i,
            constant_props_i,
            end_t_new,
            max_iterations_new,
        ) = load_hdf5(run, starting_t, dt, stepsize)

        # Check data among ranks is consistent
        if end_t is None:
            iterations_idx = iterations_idx_new
            end_t = end_t_new
            max_iterations = max_iterations_new
        else:
            assert end_t.mjd == end_t_new.mjd
            assert max_iterations == max_iterations_new
            assert (iterations_idx == iterations_idx_new).all()

        rs = rs + rs_i
        vs = vs + vs_i
        ids = ids + ids_i

        if constant_props is None:
            constant_props = constant_props_i
        else:
            constant_props = np.concatenate((constant_props, constant_props_i), axis=0)

        conj = pd.concat([conj, conj_i])
        evasions = pd.concat([evasions, evasions_i])

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

