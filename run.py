# import packages
import numpy as np
import datetime
import time
import logging
import sys
import threading
import shortuuid
from os import path
from mpi4py import MPI
from mpi4py.futures import MPIPoolExecutor
import matplotlib.pyplot as plt 
import matplotlib as mpl

# local packages
from spintronics import spintronics
from scripts import configuration_plot
from scripts.geometry import initialize as initialize_geometry
from scripts.mpi_logger import MPIFileHandler, get_logger
from scripts.interactive import InteractiveView
from scripts.io import SpintronicsSnapshoter, TaskIO, read_input_tasks, get_task_filepath

mpl.use("TkAgg")

logger = get_logger('main')


class SpintronicsDynamics(threading.Thread):
    def __init__(self, spintronics):
        threading.Thread.__init__(self)
        self.spintronics = spintronics

    def run(self):
        spintronics.metropolis(int(1e4), 0.25) 

        logger.info("Integrating")

        T, dt = 20.0, 1e-3
        spintronics.integrate(T, dt)

        logger.info("Saving log plot")

        print(len(spintronics.energy_log[:]))

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(np.linspace(4*dt, T, len(spintronics.energy_log[:])), spintronics.energy_log[:], '-')
        ax.grid(True)

        ax.set_xlabel("t")
        ax.set_ylabel("Energy")
        
        plt.show()



kekulene_vortex_path = np.array(range(32, 49), dtype=np.int8)


def run_task(task):
    initialize_geometry(
        spintronics, task['geometry']['type'], task['geometry']['parameters'])

    task_data_filepath = get_task_filepath(task["output_folder"])
    camera = SpintronicsSnapshoter(spintronics, task_data_filepath, task)

    # Themarlize
    logger.info(f'Thermalizing the system')
    t_start = time.time()

    beta = task['parameters']["beta"] if 'beta' in task['parameters'] else 1 / \
        task['parameters']["T"]

    results = spintronics.metropolis(
        int(task['parameters']['mc_steps']), beta)

    logger.info(f"Thermalization done in {time.time() - t_start:7.2f} seconds")

    # Measure
    logger.info(f'Measuring the system')
    t_start = time.time()

    vortexes = 0
    for i in range(int(task['parameters']['measures'])):
        results = spintronics.metropolis(
            int(task['parameters']['mc_inter_steps']), beta)

        has_vortex = spintronics.has_vortex(kekulene_vortex_path) != 0
        if(has_vortex):
            vortexes += 1

        camera.snapshot(beta=beta, has_vortex=has_vortex)

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


def run_interactive():
    # initialize geometry first
    sim_parms = {
        "W": 30,
        "H": 10,
        "layers": 1,
        "a": 1.0,
        "dipolar_cut": 10.0,
        "J": -1.0,
        "D": -0.2,
    }

    # sim_parms = {
    #     "L": 10,
    #     "a": 1.0,
    #     "dipolar_cut": 1.1,
    #     "J": -1.0,
    #     "D": 0.0
    # }

    initialize_geometry(spintronics, 'rectangle', sim_parms)


    thread_sim = SpintronicsDynamics(spintronics)
    #thread_int = InteractiveView(spintronics)

    thread_sim.run()
    #thread_int.start()


if __name__ == '__main__':
    run_interactive()
    
    sys.exit()

    t_start = time.time()

    input_data = read_input_tasks('./run.json')
    tasks = input_data['tasks']
    pool_size = input_data['pool_size']

    print(tasks[0])

    # resumes = TaskIO(tasks[0]['output_folder'], 'resumes.h5')

    # logger.debug(f"Starting {pool_size} mpi workers")

    # with MPIPoolExecutor(max_workers=pool_size) as executor:
    #     res = executor.map(run_task, tasks)
    #     resumes.save(list(res))
    # logger.info(f"All tasks done. Duration: {str(datetime.timedelta(seconds=time.time() - t_start))}")
