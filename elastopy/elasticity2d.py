from scipy.sparse.linalg import spsolve
from scipy import sparse
import elastopy.element2dof as element2dof
import elastopy.assemble2dof as assemble2dof
import elastopy.boundaryconditions2dof as boundaryconditions2dof
import elastopy.processing as processing
from elastopy import stiffness
from elastopy import load
from elastopy import traction
import numpy as np


def solver(model, material, b_force, trac,
           displacement_imposed, **kwargs):

    K = stiffness.K_matrix(model, material)

    Pb = load.Pb_vector(model, b_force)

    Pt = traction.Pt_vector(model, trac)

    return None, None
    # ele.b_force(b_force)

    # ele.initial_strain(kwargs)

    # P0q = assemble2dof.globalVector(ele.P0q, model)

    # P0t = boundaryconditions2dof.neumann(model, trac)

    # P0e = assemble2dof.globalVector(ele.P0e, model)

    # P0 = P0q + P0t + P0e

    # Km, P0m = boundaryconditions2dof.dirichlet(K, P0, model,
    #                                            displacement_imposed)

    # Ks = sparse.csc_matrix(Km)

    # U = spsolve(Ks, P0m)

    # sNode, sEle, eEle = processing.stress_recovery_simple2(model, U, matDic,
    #                                                        ele.e0)

    # return U, sNode
