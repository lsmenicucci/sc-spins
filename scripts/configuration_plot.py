# import packages
import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np

plt.rcParams['savefig.dpi'] = 500


def plot_current_config_3d(spintronics):
    with plt.style.context('classic'):
        fig = plt.figure()

        ax = fig.add_subplot(111, projection='3d')

        ax.scatter(spintronics.rx, spintronics.ry,
                   spintronics.rz, c='b', marker="o")
        ax.quiver(spintronics.rx, spintronics.ry, spintronics.rz,
                  spintronics.sx, spintronics.sy, spintronics.sz)
        ax.set_zlim([-2, 2])

        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        ax.set_zlabel(r'$z$')

        return fig, ax


def plot_current_config_2d(spintronics, heat_map=False):
    with plt.style.context('bmh'):
        fig = plt.figure()

        ax = fig.add_subplot(111)
        ax.set_aspect('equal')

        if(heat_map == True):
            smooth_sz = interpolate.LinearNDInterpolator(
                list(zip(spintronics.rx, spintronics.ry)), spintronics.sx)

            L = np.arange(min(spintronics.rx)*1.1,
                          max(spintronics.rx)*1.1, 0.25)
            XX, YY = np.meshgrid(L, L)
            colorbar = ax.pcolor(XX, YY, smooth_sz(
                XX, YY), shading="nearest", vmin=-1, vmax=1, cmap="PiYG")
            fig.colorbar(colorbar)

        ax.quiver(spintronics.rx, spintronics.ry,
                  spintronics.sx, spintronics.sy)

        ax.plot(spintronics.rx, spintronics.ry, 'ko', mec="1.0")

        return fig, ax


def plot_current_setup(spintronics, mode_3d=True, number_sites=False, selected_spin=None):
    if(spintronics.check_parameters(0) == 0):
        print("Can't plot, invalid parameters on module's globals")

    with plt.style.context('bmh'):
        fig = plt.figure()
        ax = None
        if(mode_3d):
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(spintronics.rx, spintronics.ry,
                       spintronics.rz, c='b', marker=".")
        else:
            ax = fig.add_subplot(111)
            ax.set_aspect('equal')
            ax.scatter(spintronics.rx, spintronics.ry, c='b', marker=".")

        if (selected_spin != None):
            if(mode_3d):
                ax.plot([spintronics.rx[selected_spin]], [spintronics.ry[selected_spin]], [
                        spintronics.rz[selected_spin]], 'ro')
            else:
                ax.plot([spintronics.rx[selected_spin]], [
                        spintronics.ry[selected_spin]], 'ro')

            # Plot exc neibs
            neib_size = spintronics.v_interacts_dip_count[selected_spin]
            neibs = spintronics.v_interacts_dip

            dip_neibs_x = list(
                map(lambda n: spintronics.rx[n - 1], neibs[:neib_size, selected_spin]))
            dip_neibs_y = list(
                map(lambda n: spintronics.ry[n - 1], neibs[:neib_size, selected_spin]))
            dip_neibs_z = list(
                map(lambda n: spintronics.rz[n - 1], neibs[:neib_size, selected_spin]))

            if(mode_3d):
                ax.scatter(dip_neibs_x, dip_neibs_y, dip_neibs_z,
                           c='g', marker='o', zorder=1)
            else:
                ax.scatter(dip_neibs_x, dip_neibs_y,
                           c='g', marker='o', zorder=1)

            # Plot exc neibs
            neib_size = spintronics.v_interacts_exc_count[selected_spin]
            neibs = spintronics.v_interacts_exc

            exc_neibs_x = list(
                map(lambda n: spintronics.rx[n - 1], neibs[:neib_size, selected_spin]))
            exc_neibs_y = list(
                map(lambda n: spintronics.ry[n - 1], neibs[:neib_size, selected_spin]))
            exc_neibs_z = list(
                map(lambda n: spintronics.rz[n - 1], neibs[:neib_size, selected_spin]))

            if(mode_3d):
                ax.scatter(exc_neibs_x, exc_neibs_y,
                           exc_neibs_z, c='r', marker="*")
            else:
                ax.scatter(exc_neibs_x, exc_neibs_y, c='r', marker="*")

        if(number_sites):
            n = len(spintronics.sx)
            for i in range(n):
                positions = [spintronics.rx[i], spintronics.ry[i] + 0.2, spintronics.rz[i]
                             ] if (mode_3d) else [spintronics.rx[i], spintronics.ry[i] + 0.2]
                ax.text(*positions,
                        f"{i + 1}",
                        verticalalignment='bottom', horizontalalignment='center',
                        fontsize=5)
        return fig, ax
