# import packages
import matplotlib.pyplot as plt

plt.rcParams['savefig.dpi'] = 500

def plot_current_config_3d(spintronics):
	with plt.style.context('classic'):
		fig = plt.figure()

		ax = fig.add_subplot(111, projection='3d')

		ax.scatter(spintronics.rx, spintronics.ry, spintronics.rz, c='b', marker="o")
		ax.quiver(spintronics.rx, spintronics.ry, spintronics.rz, spintronics.sx, spintronics.sy, spintronics.sz)
		ax.set_zlim([-2, 2])

		ax.set_xlabel(r'$x$')
		ax.set_ylabel(r'$y$')
		ax.set_zlabel(r'$z$')

		return fig, ax

def plot_current_config_2d(spintronics):
	with plt.style.context('bmh'):
		fig = plt.figure()

		ax = fig.add_subplot(111)
		ax.set_aspect('equal')

		smooth_sz = interpolate.LinearNDInterpolator(list(zip(spintronics.rx, spintronics.ry)), spintronics.sx)

		L = np.arange(min(spintronics.rx)*1.5, max(spintronics.rx)*1.5, 0.25)
		XX, YY = np.meshgrid(L,L)
		colorbar = ax.pcolor(XX, YY, smooth_sz(XX, YY), shading="nearest", vmin=-1, vmax=1, cmap="PiYG") 
		fig.colorbar(colorbar)

		ax.quiver(spintronics.rx, spintronics.ry, spintronics.sx, spintronics.sy)	

		return fig, ax


def plot_current_setup(spintronics, number_sites = False, selected_spin = None):
	if(spintronics.check_parameters(0) == 0):
		print("Can't plot, invalid parameters on module's globals")

	with plt.style.context('bmh'):
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')

		ax.scatter(spintronics.rx, spintronics.ry, spintronics.rz, c='b', marker=".")	

		if (selected_spin != None):
			ax.scatter([spintronics.rx[selected_spin]], [spintronics.ry[selected_spin]], [spintronics.rz[selected_spin]], 'ro')	

			# Plot exc neibs
			neib_size = spintronics.v_interacts_dip_count[selected_spin]
			neibs = spintronics.v_interacts_dip
			
			print(neibs[:, selected_spin])

			dip_neibs_x = list(map(lambda n: spintronics.rx[n - 1], neibs[:neib_size, selected_spin]))
			dip_neibs_y = list(map(lambda n: spintronics.ry[n - 1], neibs[:neib_size, selected_spin]))
			dip_neibs_z = list(map(lambda n: spintronics.rz[n - 1], neibs[:neib_size, selected_spin]))
			ax.scatter(dip_neibs_x, dip_neibs_y, dip_neibs_z, c='g', marker='o', zorder = 1)

			# Plot exc neibs
			neib_size = spintronics.v_interacts_exc_count[selected_spin]
			neibs = spintronics.v_interacts_exc

			print(neibs[:, selected_spin])
			
			exc_neibs_x = list(map(lambda n: spintronics.rx[n - 1], neibs[:neib_size, selected_spin]))
			exc_neibs_y = list(map(lambda n: spintronics.ry[n - 1], neibs[:neib_size, selected_spin]))
			exc_neibs_z = list(map(lambda n: spintronics.rz[n - 1], neibs[:neib_size, selected_spin]))

			ax.scatter(exc_neibs_x, exc_neibs_y, exc_neibs_z, c='r', marker="*")

		if(number_sites):
			n = len(spintronics.sx)
			for i in range(n):
				ax.text( spintronics.rx[i], spintronics.ry[i] + 0.2, spintronics.rz[i],
					f"{i + 1}",
	 				verticalalignment='bottom', horizontalalignment='center',
					fontsize = 5)

		return fig, ax
