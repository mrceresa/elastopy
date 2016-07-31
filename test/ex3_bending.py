import numpy as np
from elastopy import statics
from elastopy import gmsh
from elastopy import plotter
from elastopy import data

meshName = 'patch5'

model = gmsh.Parse(meshName)

material = data.Collect()
s = list(model.surf.keys())

material.E[s[0]] = 1000
material.nu[s[0]] = 0.3


def b_force(x1, x2, t=1):
    return np.array([
        0.0,
        0.0,
    ])

I = 0.01
M = 10.0
h = 0.5


def trac_bc(x1, x2, t=1):
    return {
        ('line', 1): [M*(x2-h)/I, 0.0]
    }


def displ_bc(x1, x2):
    return {
        ('line', 3): [0.0, 0.0]

    }

U, SIG = statics.solver(model, material, b_force,
                        trac_bc, displ_bc)

plotter.stresses(model, SIG, spmin=True)
plotter.stresses(model, SIG, spmax=True)
plotter.stresses(model, SIG, s11=True)
plotter.model(model, ele=True, edges_label=True)
plotter.model_deformed(model, U, magf=0.1)
plotter.show()
