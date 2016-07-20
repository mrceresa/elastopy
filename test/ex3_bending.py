import numpy as np
from elastopy import elasticity2d
from elastopy import gmsh
from elastopy import plotter

meshName = 'patch5'

mesh = gmsh.Parse(meshName)

material = {'E-nu': [1000.0, 0.3]}


def body_forces(x1, x2):
    return np.array([
        0.0,
        0.0,
    ])

I = 0.01
M = 10.0
h = 0.5


def traction_imposed(x1, x2):
    return {
        ('line', 1): [M*(x2-h)/I, 0.0]
    }


def displacement_imposed(x1, x2):
    return {
        ('line', 3): [0.0, 0.0]

    }

U, sNode = elasticity2d.solver(mesh, material, body_forces,
                               traction_imposed, displacement_imposed)

plotter.stress(mesh, sNode, spmin=True)
plotter.stress(mesh, sNode, spmax=True)
plotter.stress(mesh, sNode, s11=True)
plotter.model(mesh, ele=True, edges_label=True)
plotter.model_deformed(mesh, U, magf=0.1)
plotter.show()
