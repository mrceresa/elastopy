import numpy as np
from elastopy import dynamics
from elastopy import gmsh
from elastopy import plotter
from elastopy import data

mesh_name = 'patch'

model = gmsh.Parse(mesh_name)

material = data.Collect()
material.E[9] = 1000
material.nu[9] = 0.2
material.dnsty[9] = 1.


def b_force(x1, x2, t=1):
    return np.array([0.0,
                     -1.0])


def trac_bc(x1, x2, t=1):
    return {('line', 3): [-1., 0.0],
            ('line', 1): [1., 0.0]}


def displ_bc(x1, x2):
    return {('node', 0): [0.0, 0.0],
            ('node', 1): ['free', 0.0]}

t_int = 0.1
dt = 0.1/300

frm, fig = dynamics.solver(model, material, b_force,
                           trac_bc, displ_bc, t_int, dt,
                           plot_u=True, vmin=0., vmax=60., magf=100)

ani = plotter.anime(frm, fig, t_int)
ani.save('elastodyna_displ4.gif', writer='imagemagick', bitrate=300)
plotter.show()
