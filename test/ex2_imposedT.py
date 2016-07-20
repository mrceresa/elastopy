import numpy as np
from elastopy import elasticity2d
from elastopy import gmsh
from elastopy import plotter

meshName = 'patch'

mesh = gmsh.Parse(meshName)

material = {'E-nu': [1000.0, 0.3]}


def body_forces(x1, x2):
    return np.array([
        0.0,
        0.0,
    ])


def traction_imposed(x1, x2):
    return {
        ('line', 3): [-1.0, 0.0],
        ('line', 1): [1.0, 0.0]
    }


def displacement_imposed(x1, x2):
    return {
        ('node', 0): [0.0, 0.0],
        ('node', 1): ['free', 0.0]

    }

U, sNode = elasticity2d.solver(mesh, material, body_forces,
                               traction_imposed, displacement_imposed)

plotter.model(mesh, ele=True, nodes_label=True,
              ele_label=True, edges_label=True)
plotter.model_deformed(mesh, U, magf=100, ele=True)

plotter.show()
