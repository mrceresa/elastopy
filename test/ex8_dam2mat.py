import numpy as np
from elastopy import elasticity2d
from elastopy import gmsh
from elastopy import plotter


meshName = 'barragem3'

mesh = gmsh.Parse(meshName)

material = {'Ef-nu': [30e9, 0.17], 'Eb-nu': [15e9, 0.17]}


def body_forces(x1, x2):
    return np.array([
        0.0,
        (x2 >= 70.0)*-24e3+(x2 < 70.0)*-30e3,
    ])

gma = 9820.0


def traction_imposed(x1, x2):
    return {
        ('line', 7): [(x2 <= 186.13)*(gma*(186.13-x2)), 0.0],
        ('line', 8): [(17.07*gma - 116.13*gma)*(x2-70.0)/(169.06-70.0) +
                      116.13*gma,
                      -((17.07*gma - 116.13*gma)*(x2-70.0)/(169.06-70.0) +
                        116.13*gma)],
        ('line', 9): [0.0, - 116.13*gma]
    }


def displacement_imposed(x1, x2):
    return {
        ('line', 0): [0.0, 0.0],
        ('line', 10): [0.0, 'free'],
        ('line', 1): [0.0, 'free']
    }

U, sNode = elasticity2d.solver(mesh, material, body_forces,
                               traction_imposed, displacement_imposed)

plotter.stress(mesh, sNode, spmin=True, ftr=1000)
plotter.stress(mesh, sNode, spmax=True, ftr=1000)
plotter.model(mesh, edges_label=True, surf_label=True)
plotter.model(mesh, ele=True)
plotter.model_deformed(mesh, U, magf=1000)
plotter.show()
