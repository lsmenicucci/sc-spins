# import packages
import time
import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt
import logging
from mpi4py import MPI
import sys

# local packages
from spintronics import spintronics
from scripts import geometry as spintronicGeometries
from scripts.mpi_logger import MPIFileHandler

# Get MPI communicator
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Start logging
logger = logging.getLogger("rank[%i]"%comm.rank)
logger.setLevel(logging.DEBUG)

mh = MPIFileHandler("spintronics.log")
formatter = logging.Formatter('%(asctime)s:%(name)s:%(levelname)s:%(message)s')
mh.setFormatter(formatter)

stdout_h = logging.StreamHandler(sys.stdout)
stdout_h.setFormatter(formatter)

logger.addHandler(mh)
logger.addHandler(stdout_h)
logger.debug('Starting worker')


plt.rcParams['savefig.dpi'] = 500

L = 10
a = 1.0
dipolar_cut = 1.1
J = -1.0
D = -1.0
beta = 1/1.0

spintronicGeometries.initialize_disk(spintronics, 7, 1, a, dipolar_cut, J, D)
#spintronicGeometries.initialize_line(spintronics, 2, a, dipolar_cut, J, D)


logger.info(f'Thermalizing the system')
t_start = time.time()

spintronics.metropolis(int(1e4), beta)

logger.info(f"Thermalization done in {time.time() - t_start:7.2f} seconds")

logger.info(f"Energy before integration = {spintronics.total_energy():12.6f}")
spintronics.integrate(10, 0.005)
logger.info(f"Energy after integration = {spintronics.total_energy():12.6f}")

# logger.info(f'Measuring the system')
# t_start = time.time()
# results = spintronics.metropolis(int(1e5), beta)

# logger.info(f"Measurement done in {time.time() - t_start:7.2f} seconds")
# logger.debug(f"Measurement output is: {results}")

# print(f"- Mean energy: {results[0]:12.6f}")
# print(f"- Mean mag_x:  {results[1]:12.6f}")
# print(f"- Mean mag_y:  {results[2]:12.6f}")
# print(f"- Mean mag_z:  {results[3]:12.6f}")

# print("Plotting")
# fig, ax = plot_current_config_2d()

# fig.suptitle(fr'$D = {D}, Dcut = {dipolar_cut:4.2f}, \beta = {beta}$')
# plt.savefig(f'test - {rank}.png')