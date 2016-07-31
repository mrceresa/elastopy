import numpy as np
from elastopy import statics
from elastopy import gmsh
from elastopy import plotter
from elastopy import data

mesh_name = 'patch'

model = gmsh.Parse(mesh_name)


material = data.Collect()
s = list(model.surf.keys())

material.E[s[0]] = 1000
material.nu[s[0]] = 0.3


def body_forces(x1, x2, t=1):
    return np.array([0.0, 0.0])


def traction_imposed(x1, x2, t=1):
    return {}


def displacement_imposed(x1, x2):
    return {
        ('nodes', 0, 3, 7): [0.0, 0.0],
        ('nodes', 4, 6): [0.5, 0.0],
        ('nodes', 1, 5, 2): [1.0, 0.0],
    }


U, SIG = statics.solver(model, material, body_forces,
                        traction_imposed, displacement_imposed)

plotter.model(model)
plotter.model_deformed(model, U, magf=0.1)

plotter.show()
