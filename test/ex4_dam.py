import numpy as np
from elastopy import elasticity2d
from elastopy import gmsh
from elastopy import plotter

meshName = 'barragem2'

material = {'E-nu': [34473.786*1e6, 0.1]}

mesh = gmsh.Parse(meshName)


def body_forces(x1, x2): 
    return np.array([
        0.0,
        -24357.0,
    ])

gma = 9820.0*10


def traction_imposed(x1, x2):
    return {
        ('line', 7): [(x2 <= 186.13)*(-gma*(x2-186.13)), 0.0],
        ('line', 8): [(17.07*gma - 186.13*gma)*x2/169.06 + 186.13*gma,
                      -((17.07*gma - 186.13*gma)*x2/169.06 + 186.13*gma)],
        ('line', 9): [0.0, - 186.13*gma]
    }


def displacement_imposed(x1, x2):
    return {
        ('line', 0): [0.0, 0.0],
        ('line', 1): [0.0, 0.0],
        ('line', 10): [0.0, 0.0]
    }

U, sNode = elasticity2d.solver(mesh, material, body_forces,
                               traction_imposed, displacement_imposed)

# plotter.model(mesh)
# plotter.model_deformed(mesh, U, magf=0.1)
plotter.stress(mesh, sNode, spmin=True)
plotter.stress(mesh, sNode, spmax=True)
plotter.show()
