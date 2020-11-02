# import packages
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl
import time
import logging
import sys
import threading
from mpi4py import MPI
from mpi4py.futures import MPIPoolExecutor

# local packages
from spintronics import spintronics
from scripts import configuration_plot
from scripts.geometry import initialize as initialize_geometry
from scripts.mpi_logger import MPIFileHandler, get_logger
from scripts.interactive import InteractiveView
from scripts.io import SpintronicsSnapshoter, TaskIO, read_input_tasks

mpl.use("TkAgg")

logger = get_logger('main')


class SpintronicsDynamics(threading.Thread):
    def __init__(self, spintronics):
        threading.Thread.__init__(self)
        self.spintronics = spintronics

    def run(self):
        logger.info(f'Thermalizing the system')
        t_start = time.time()

        spintronics.metropolis(int(1e5), beta)

        logger.info(
            f"Thermalization done in {time.time() - t_start:7.2f} seconds")

        logger.info(
            f"Energy before integration = {spintronics.total_energy():12.6f}")
        #spintronics.integrate(10000, 0.05)
        logger.info(
            f"Energy after integration = {spintronics.total_energy():12.6f}")


plt.rcParams['savefig.dpi'] = 500

L = 10
a = 1.0
dipolar_cut = 15.0
J = -1.0
D = -1.0
beta = 1 / 1.0


def run_task(task):
    initialize_geometry(
        spintronics, task['geometry']['type'], task['geometry']['parameters'])

    camera = SpintronicsSnapshoter(spintronics, 'test.h5', task)
    camera.snapshot()

    logger.info(f'Measuring the system')
    t_start = time.time()

    results = spintronics.metropolis(int(1e3), beta)
    camera.snapshot(beta=beta)

    logger.info(f"Measurement done in {time.time() - t_start:7.2f} seconds")
    logger.debug(f"Measurement output is: {results}")
    camera.close()

    return {}


if __name__ == '__main__':
    task_data = read_input_tasks('./run.json')
    tasks = []
    for group in task_data["task_groups"]:
        for task in group['tasks']:
            taskGeometryParameters = {**group['geometry']["parameters"], **task['geometric_parms']} if('geometric_parms' in task) else group['geometry']["parameters"]
            tasks.append({
                "parameters": task["parameters"],
                "geometry": {**group['geometry'], "parameters": taskGeometryParameters}
            })

    run_task(tasks[0])
    sys.exit()
    with MPIPoolExecutor(max_workers=2) as executor:
        res = executor.submit(run_task, tasks[0])
        res.result()
