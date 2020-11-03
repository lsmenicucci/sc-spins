# import packages
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl
import time
import logging
import sys
import threading
import shortuuid
from os import path
from mpi4py import MPI
from mpi4py.futures import MPIPoolExecutor

# local packages
from spintronics import spintronics
from scripts import configuration_plot
from scripts.geometry import initialize as initialize_geometry
from scripts.mpi_logger import MPIFileHandler, get_logger
from scripts.interactive import InteractiveView
from scripts.io import SpintronicsSnapshoter, TaskIO, read_input_tasks, get_task_filepath

# mpl.use("TkAgg")
mpl.use("Agg")

logger = get_logger('main')


class SpintronicsDynamics(threading.Thread):
    def __init__(self, spintronics):
        threading.Thread.__init__(self)
        self.spintronics = spintronics

    def run(self):
        results = spintronics.metropolis(int(1e5), 0.25)


plt.rcParams['savefig.dpi'] = 500

kekulene_vortex_path = np.array(range(32, 49), dtype=np.int8)


def run_task(task):
    initialize_geometry(
        spintronics, task['geometry']['type'], task['geometry']['parameters'])

    task_data_filepath = get_task_filepath(task["output_folder"])
    camera = SpintronicsSnapshoter(spintronics, task_data_filepath, task)

    # Themarlize
    logger.info(f'Thermalizing the system')
    t_start = time.time()

    results = spintronics.metropolis(
        int(task['parameters']['mc_steps']), task['parameters']["beta"])

    logger.info(f"Thermalization done in {time.time() - t_start:7.2f} seconds")

    # Measure
    logger.info(f'Measuring the system')
    t_start = time.time()

    vortexes = 0
    for i in range(int(task['parameters']['measures'])):
        results = spintronics.metropolis(
            int(task['parameters']['mc_inter_steps']), task['parameters']["beta"])

        has_vortex = spintronics.has_vortex(kekulene_vortex_path) != 0
        if(has_vortex):
            vortexes += 1

        camera.snapshot(beta=task['parameters']["beta"], has_vortex=has_vortex)

    logger.info(f"Measurement done in {time.time() - t_start:7.2f} seconds")

    # Gather results
    measures_zip = zip(["mean_energy", "mean_mag_x", "mean_mag_y", "mean_mag_z"],
                       spintronics.metropolis_measures)
    metropolis_measures = {n: v for (n, v) in measures_zip}
    results = {
        **task['parameters'],
        **metropolis_measures,
        "vortex_density": vortexes/task['parameters']['measures'],
        "filepath": task_data_filepath
    }

    camera.update_attrs(**results)
    camera.close()

    return {
        "geometry": task['geometry'],
        **results
    }


if __name__ == '__main__':
    tasks = read_input_tasks('./run.json')
    resumes = TaskIO(tasks[0]['output_folder'], 'resumes.h5')

    with MPIPoolExecutor(max_workers=2) as executor:
        res = executor.map(run_task, tasks)
        resumes.save(list(res))
