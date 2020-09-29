# import packages
import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt
from spintronics import *

def get_sites_in_range(center_site, cut_range, include_border = False):
	n = len(spintronics.rx)

	site_x = spintronics.rx[center_site]
	site_y = spintronics.ry[center_site]	

	def dist(i): return np.sqrt( (site_x - spintronics.rx[i])**2 + (site_y - spintronics.ry[i])**2 )

	if(include_border):
		return list(filter(lambda s: s != center_site and dist(s) <= cut_range , range(n))) 
	else:
		return list(filter(lambda s: s != center_site and dist(s) < cut_range, range(n))) 

def initialize_periodic_square(n, a, J, D):
	# Initialize vectors
	spintronics.sx = np.ones(n*n)
	spintronics.sy = np.ones(n*n)
	spintronics.sz = np.ones(n*n)

	spintronics.rx = np.zeros(n*n)
	spintronics.ry = np.zeros(n*n)

	# Calculate exchange potential
	spintronics.v_exc = J * np.ones([3, 4, n*n])	
	spintronics.v_interacts_exc = -1 * np.ones([4, n*n])	
	v_interacts_exc_count = -1 * np.ones(n*n)

	# Calculate dipolar potential
	spintronics.v_dip = D * np.ones([3, 4, n*n])
	spintronics.v_interacts_dip = -1 * np.ones([4, n*n])
	spintronics.v_interacts_dip_count = -1 * np.ones(n*n)

	for i in range(n*n):
		# Set neibs
		spintronics.v_interacts_exc[0, i] = int(i / n) * n + (i + 1) % n
		spintronics.v_interacts_exc[1, i] = int(i / n) * n + (i - 1) % n
		spintronics.v_interacts_exc[2, i] = (i + n) % (n**2)
		spintronics.v_interacts_exc[3, i] = (i - n) % (n**2)

		# Set position
		spintronics.rx[i] = (i % n) * a 
		spintronics.ry[i] = int(i / n) * a 

		# Set dipolar interaction range
		spintronics.v_interacts_dip_count[i] = 0
		v_interacts_exc_count[i] = 4


	spintronics.v_interacts_exc = spintronics.v_interacts_exc + 1

def initialize_disk(r, a, dipolar_cut, J, D):
	slimit = r*int(np.round(float(r)/2.0))
	points = [ (i, j) for i in range(-slimit, slimit + 1) for j in range(-slimit, slimit + 1) if i**2 + j**2 < r**2  ]
	
	n = len(points)

	# Set the sites positions
	spintronics.rx = list(map(lambda p: p[0]*a, points))
	spintronics.ry = list(map(lambda p: p[1]*a, points))
	
	# Initialize spin vectors
	spintronics.sx = np.ones(n)
	spintronics.sy = np.ones(n)
	spintronics.sz = np.ones(n)	

	# get neibs
	exc_neibs = list(map(lambda s: get_sites_in_range(s, a, include_border = True), range(n)))
	dip_neibs = list(map(lambda s: get_sites_in_range(s, dipolar_cut), range(n)))

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

def plot_current_setup(number_sites = False,selected_spin = None):
	if(spintronics.check_parameters(0) == 0):
		print("Can't plot, invalid parameters on module's globals")

	with plt.style.context('bmh'):
		plt.rcParams['savefig.dpi'] = 500
		plt.axes().set_aspect('equal')
		plt.plot(spintronics.rx, spintronics.ry, 'bo')	

		if (selected_spin != None):
			plt.plot([spintronics.rx[selected_spin]], [spintronics.ry[selected_spin]], 'ro')	

			# Plot exc neibs
			neib_size = spintronics.v_interacts_exc_count[selected_spin]
			neibs = spintronics.v_interacts_exc

			print(neibs[:, selected_spin])
			
			exc_neibs_x = list(map(lambda n: spintronics.rx[n - 1], neibs[:neib_size, selected_spin]))
			exc_neibs_y = list(map(lambda n: spintronics.ry[n - 1], neibs[:neib_size, selected_spin]))

			print(exc_neibs_x)
			print(exc_neibs_y)
			plt.plot(exc_neibs_x, exc_neibs_y, 'rx')

			# Plot exc neibs
			neib_size = spintronics.v_interacts_dip_count[selected_spin]
			neibs = spintronics.v_interacts_dip
			
			dip_neibs_x = list(map(lambda n: spintronics.rx[n - 1], neibs[:neib_size, selected_spin]))
			dip_neibs_y = list(map(lambda n: spintronics.ry[n - 1], neibs[:neib_size, selected_spin]))
			plt.plot(dip_neibs_x, dip_neibs_y, 'w+')

		if(number_sites):
			n = len(spintronics.sx)
			for i in range(n):
				plt.text( spintronics.rx[i], spintronics.ry[i] + 0.2, 
					f"{i}",
	 				verticalalignment='bottom', horizontalalignment='center',
					fontsize = 5)

		plt.savefig('test.png')

def plot_current_config():
	with plt.style.context('bmh'):
		plt.rcParams['savefig.dpi'] = 500
		plt.axes().set_aspect('equal')


		smooth_sz = interpolate.LinearNDInterpolator(list(zip(spintronics.rx, spintronics.ry)), spintronics.sx)

		L = np.arange(min(spintronics.rx)*1.5, max(spintronics.rx)*1.5, 0.25)
		XX, YY = np.meshgrid(L,L)
		#plt.pcolor(XX, YY, smooth_sz(XX, YY), shading="nearest", vmin=-1, vmax=1, cmap="PiYG") 
		#plt.colorbar()

		plt.quiver(spintronics.rx, spintronics.ry, spintronics.sx, spintronics.sy)	

		plt.savefig('test.png')

r = 5
a = 1.0
dipolar_cut = 10.0
J = -1.0
D = 0.1

print("Initializing system topology")
initialize_disk(r, a, dipolar_cut, J, D)

print("Thermalizing")
mean_energ = spintronics.metropolis(int(1e5), 1/.7)
print(f"Mean energy: {mean_energ}")
print(f"Mean Mz: {np.mean(spintronics.sz)}")

print("Plotting")
plot_current_config()

