from contextlib import ContextDecorator
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d
import mpl_toolkits.mplot3d.axis3d
from mpl_toolkits.mplot3d.axis3d import Axis
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.text import TextPath, Text
from matplotlib.transforms import Affine2D
from matplotlib.patches import Rectangle, Circle, PathPatch
from matplotlib.font_manager import FontProperties

# Patch for enabling propper limits on 3D plots
if not hasattr(Axis, "_get_coord_info_old"):
    def _get_coord_info_new(self, renderer):
        mins, maxs, centers, deltas, tc, highs = self._get_coord_info_old(
            renderer)
        mins += deltas / 4
        maxs -= deltas / 4
        return mins, maxs, centers, deltas, tc, highs
    Axis._get_coord_info_old = Axis._get_coord_info
    Axis._get_coord_info = _get_coord_info_new


def text3d(ax, xyz, s, zdir="z", size=None, angle=0, usetex=False, **kwargs):

    x, y, z = xyz
    if zdir == "y":
        xy1, z1 = (x, z), y
    elif zdir == "x":
        xy1, z1 = (y, z), x
    else:
        xy1, z1 = (x, y), z

    text_path = TextPath((0, 0), s, size=size, usetex=usetex)
    trans = Affine2D().rotate(angle).translate(xy1[0], xy1[1])

    p1 = PathPatch(trans.transform_path(text_path), antialiased = True, lw = 0)
    ax.add_patch(p1)
    art3d.pathpatch_2d_to_3d(p1, z=z1, zdir=zdir)


def decorate_3d_axis(axis_obj):
    axis_obj.pane.fill = False
    axis_obj.pane.set_alpha(0)
    axis_obj.pane.set_edgecolor('k')
    axis_obj.pane.set_closed(True)
    axis_obj._axinfo['axisline']['linewidth'] = 0.0
    axis_obj._axinfo['tick']['linewidth'][True] = 0.0
    axis_obj.line.set_linewidth(0)


def missing_direction(one, other):
    return [x for x in range(3) if x != one and x != other][0]


def plot_rect(ax, direction, zdiretion, center, width, height, *args, **kwargs):
    edges = [center.copy() for k in range(4)]

    other_direction = missing_direction(direction, zdiretion)

    # width
    edges[2][other_direction] += width
    edges[3][other_direction] += width

    # heigt
    edges[1][direction] += height
    edges[2][direction] += height

    return ax.add_collection3d(Poly3DCollection([edges], *args, **kwargs))


def plot_manual_frame(ax, lw=0.5, tick_size=2.0, tick_thickness=0.8):
    lw = 1e-2 * lw * mpl.rcParams['axes.linewidth']
    tick_size = 1e-2 * tick_size * mpl.rcParams['axes.linewidth']
    tick_thickness = 1e-2 * tick_thickness * mpl.rcParams['axes.linewidth']

    # Plot the spines
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    zlim = ax.get_zlim()

    # calculate the plot dimensions
    dimensions = [
        xlim[1] - xlim[0],
        ylim[1] - ylim[0],
        zlim[1] - zlim[0]
    ]

    # set facings
    facings = [2, 2, 0]
    centers = [
        np.zeros(3),
        np.array(
            [dimensions[0], 0, 0]),
        np.zeros(3)]

    spine_zip = zip(dimensions, facings, centers)

    for i, (dim, zidx, center) in enumerate(spine_zip):
        other_direction = missing_direction(i, zidx)
        spine_width = dimensions[other_direction] * lw
        spine_height = dim
        tick_width = dim * tick_size

        # plot_rect(ax, edges, shade=False, edgecolor='none',
        #           antialiased=True, color='k', clip_on=False)
        plot_rect(ax, i, zidx, center, spine_width,
                  spine_height, color='k')

    # Plot the ticks
    xticks = ax.get_xticks()
    yticks = ax.get_yticks()
    zticks = ax.get_zticks()

    tick_zip = zip(dimensions, facings, [xticks, yticks, zticks], centers)

    for i, (axis_dim, zidx, ticks, center) in enumerate(tick_zip):
        other_direction = missing_direction(i, zidx)

        tick_width = axis_dim * tick_thickness
        tick_height = dimensions[other_direction] * tick_size

        inout_factor = 1 if center[other_direction] > 0.5 * axis_dim else -1

        for tick in ticks:
            tick_center = center.copy()
            tick_center[i] += tick

            plot_rect(ax, other_direction, zidx, tick_center,
                      tick_width, inout_factor * tick_height, color='k')

            # plot tick label
            tick_label_center = [0, 0, 0]
            tick_label_center[zidx] = tick_center[zidx]
            tick_label_center[i] = tick + inout_factor * tick_height * 0.25

            if(inout_factor == 1):
                tick_label_center[other_direction] = dimensions[other_direction]

            tick_label_center[other_direction] += inout_factor * dimensions[other_direction] * tick_width
            tick_dir = ['x', 'y', 'z'][zidx]
            angle = np.pi/2 if i == 0 else 0

            text3d(ax, tick_label_center, f"{tick:2.0f}".replace(' ', ''), 
                zdir=tick_dir,
                angle=angle, 
                size= tick_width * dimensions[i], fc='k')

        
    return


class Custom3DPlot(ContextDecorator):
    def __init__(self, fig, subplot_args=()):
        if(fig == None):
            self.fig = plt.figure()
        else:
            self.fig = fig

        self.ax = fig.add_subplot(
            *subplot_args, projection='3d', proj_type='ortho')

    def __enter__(self):
        return self.ax

    def __exit__(self, *exc):
        ax.set_xmargin(0)
        ax.set_ymargin(0)
        ax.set_zmargin(0)

        for axis_obj in [ax.w_xaxis, ax.w_yaxis, ax.w_zaxis]:
            decorate_3d_axis(axis_obj)

        ax.w_zaxis._axinfo['juggled'] = (1, 2, 1)

        plot_manual_frame(ax)


if __name__ == '__main__':
    mpl.use("TkAgg")

    line = np.linspace(0, 1, 40)
    XX, YY = np.meshgrid(line, line)
    ZZ = np.sin(XX**2 + YY**2)

    fig = plt.figure(figsize=(10, 10))

    with Custom3DPlot(fig, (1, 1, 1)) as ax:
        # ax.plot_surface(XX, YY, ZZ)

        ax.set_zlim([0, 1])
        ax.set_xlim([0, 10])
        ax.set_ylim([0, 5])

        ax.grid(True, wich='major')

    plt.show()
