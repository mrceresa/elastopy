import elastopy.element2dof as element2dof
import elastopy.assemble2dof as assemble2dof
import elastopy.boundaryconditions2dof as boundaryconditions2dof
import elastopy.processing as processing
from elastopy import stiffness
from elastopy import load
from elastopy import traction
from elastopy import stress
import numpy as np


def solver(model, material, b_force, trac_bc, displ_bc, EPS0=0, t=1):

    K = stiffness.K_matrix(model, material, t)

    Pb = load.Pb_vector(model, b_force, t)

    Pt = traction.Pt_vector(model, trac_bc, t)

    Pe = load.Pe_vector(model, material, EPS0, t)

    P = Pb + Pt + Pe

    Km, Pm = boundaryconditions2dof.dirichlet(K, P, model, displ_bc)

    U = np.linalg.solve(Km, Pm)

    SIG = stress.recovery(model, material, U, EPS0)

    return U, SIG
