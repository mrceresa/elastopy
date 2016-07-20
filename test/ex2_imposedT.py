import numpy as np
from elastopy import elasticity2d
from elastopy import gmsh
from elastopy import plotter
from elastopy import data

mesh_name = 'patch'

model = gmsh.Parse(mesh_name)

material = data.Collect()

material.E[9] = 1000
material.nu[9] = 0.3


def body_forces(x1, x2):
    return np.array([
        0.0,
        0.0])


def traction_imposed(x1, x2):
    return {
        ('line', 3): [-1.0, 0.0],
        ('line', 1): [1.0, 0.0]}


def displacement_imposed(x1, x2):
    return {
        ('node', 0): [0.0, 0.0],
        ('node', 1): ['free', 0.0]}

U, sNode = elasticity2d.solver(model, material, body_forces,
                               traction_imposed, displacement_imposed)

# plotter.model(model, ele=True, nodes_label=True,
#               ele_label=True, edges_label=True)
# plotter.model_deformed(model, U, magf=100, ele=True)

# plotter.show()
