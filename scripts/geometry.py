# import packages
import logging
import numpy as np 
from mpi4py import MPI

# import local modules
from scripts.mpi_logger import MPIFileHandler

# Get MPI communicator
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Start logging
logger = logging.getLogger(f"rank[{comm.rank}](geometry)")
logger.setLevel(logging.DEBUG)

mh = MPIFileHandler("spintronics.log")
formatter = logging.Formatter('%(asctime)s:%(name)s:%(levelname)s:%(message)s')
mh.setFormatter(formatter)

logger.addHandler(mh)
logger.debug('Starting geometry module')


def get_sites_in_range(spintronics, center_site, cut_range, include_border = False):
	n = len(spintronics.rx)

	site_x = spintronics.rx[center_site]
	site_y = spintronics.ry[center_site]	
	site_z = spintronics.rz[center_site]	

	def dist(i): 
		return np.sqrt( 
			(site_x - spintronics.rx[i])**2 + 
			(site_y - spintronics.ry[i])**2 + 
			(site_z - spintronics.rz[i])**2)

	if(include_border):
		return list(filter(lambda s: s != center_site and dist(s) <= cut_range , range(n))) 
	else:
		return list(filter(lambda s: s != center_site and dist(s) < cut_range, range(n))) 

def initialize_general_geometry(spintronics, points, a, dipolar_cut, J, D):
	logger.debug(f'Initializin general geometry with {len(points)} points ')
	n = len(points)

	# Set the sites positions
	spintronics.rx = list(map(lambda p: p[0], points))
	spintronics.ry = list(map(lambda p: p[1], points))
	spintronics.rz = list(map(lambda p: p[2], points))
	
	# Initialize spin vectors
	spintronics.sx = np.ones(n)
	spintronics.sy = np.ones(n)
	spintronics.sz = np.ones(n)	

	# get neibs
	exc_neibs = list(map(lambda s: get_sites_in_range(spintronics, s, a, include_border = True), range(n)))
	dip_neibs = list(map(lambda s: get_sites_in_range(spintronics, s, dipolar_cut), range(n)))

	max_exc_neibs_count = max(map(lambda neibs: len(neibs), exc_neibs))
	max_exc_dip_count = max(map(lambda neibs: len(neibs), dip_neibs))

	# Initalize exchange interaction vectors
	spintronics.v_exc = J * np.ones([3, max_exc_neibs_count, n])	
	spintronics.v_interacts_exc = -1 * np.ones([max_exc_neibs_count, n])	
	spintronics.v_interacts_exc_count = -1 * np.ones(n)

	# Initialize dipolar interaction vectors
	spintronics.v_dip = D * np.ones([3, max_exc_dip_count, n])
	spintronics.v_interacts_dip = -1 * np.ones([max_exc_dip_count, n])
	spintronics.v_interacts_dip_count = -1 * np.ones(n)

	logger.debug('Calculating neibs')

	# Load neib data
	for i in range(n):
		exc_neib_size = len(exc_neibs[i])
		dip_neib_size = len(dip_neibs[i])

		spintronics.v_interacts_exc_count[i] = exc_neib_size
		for exc_j in range(exc_neib_size):
			spintronics.v_interacts_exc[exc_j, i] = exc_neibs[i][exc_j] + 1

		spintronics.v_interacts_dip_count[i] = dip_neib_size
		for dip_j in range(dip_neib_size):
			spintronics.v_interacts_dip[dip_j, i] = dip_neibs[i][dip_j] + 1

	logger.debug('Done setting geometry')

def initialize_disk(spintronics, r, layers, a, dipolar_cut, J, D):
	logger.info(f'Initializing disk geometry with {layers} layers and radius = {r}')

	slimit = r*int(np.round(float(r)/2.0))
	points = [ (i*a, j*a, k*a) 
		for k in range(layers)
		for j in range(-slimit, slimit + 1) 
		for i in range(-slimit, slimit + 1) 
			if i**2 + j**2 < r**2  ]

	initialize_general_geometry(spintronics, points, a, dipolar_cut, J, D)

def initialize_square_layers(spintronics, L, layers, a, dipolar_cut, J, D):
	logger.info(f'Initializing square geometry with {layers} layers and side = {L}')

	slimit = int(np.round(float(L)/2.0))
	points = [ (i*a, j*a, k*a) 
		for k in range(layers)
		for j in range(-slimit, slimit + 1) 
		for i in range(-slimit, slimit + 1) ]

	initialize_general_geometry(spintronics, points, a, dipolar_cut, J, D)
