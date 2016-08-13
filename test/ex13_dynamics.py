import numpy as np
from elastopy import dynamics
from elastopy import gmsh
from elastopy import plotter
from elastopy import data

mesh_name = 'patch'

model = gmsh.Parse(mesh_name)

material = data.Collect()
s = list(model.surf.keys())

material.cndtvt[s[0]] = 1.
material.spcfht[s[0]] = 1.
material.dnsty[s[0]] = 0.1
material.E[s[0]] = 1000
material.nu[s[0]] = 0.2
material.dnsty[s[0]] = 1.
material.thrmexp[s[0]] = 1e-6


def b_force(x1, x2, t=1):
    return np.array([0.0,
                     0.0])


def trac_bc(x1, x2, t=1):
    return {('line', 3): [0., 0.0],
            ('line', 1): [0., 0.0]}


def displ_bc(x1, x2):
    return {('node', 0): [0.0, 0.0],
            ('node', 1): ['free', 0.0]}


def b_heat(x1, x2, t=1):
    return t*10000


def flux_bc(x1, x2, t=1):
    return {1: 0.0,
            2: 0.0,
            0: 0.0}


def temp_bc(x1, x2, t=1):
    return {3: 1.}

t_int = 30
dt = t_int/300

U, SIG = dynamics.solver(model, material, b_force,
                         trac_bc, displ_bc, t_int, dt,
                         b_heat=b_heat, flux_bc=flux_bc, temp_bc=temp_bc)

plotter.stress_animation(SIG, model, t_int, dt, spmax=True,
                         name="stre9.gif", vmin=-.5, vmax=.5,
                         interval=100)

plotter.displ_animation(U, model, t_int, dt, name='displ9.gif', magf=500,
                        interval=100)
