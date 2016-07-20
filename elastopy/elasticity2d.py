from scipy.sparse.linalg import spsolve
from scipy import sparse
import elastopy.element2dof as element2dof
import elastopy.assemble2dof as assemble2dof
import elastopy.boundaryconditions2dof as boundaryconditions2dof
import elastopy.processing as processing
from elastopy import stiffness


def solver(model, material, body_forces, traction_imposed,
           displacement_imposed, **kwargs):

    K = stiffness.K_matrix(model, material)

    return None, None
    # ele.body_forces(body_forces)

    # ele.initial_strain(kwargs)

    # P0q = assemble2dof.globalVector(ele.P0q, model)

    # P0t = boundaryconditions2dof.neumann(model, traction_imposed)

    # P0e = assemble2dof.globalVector(ele.P0e, model)

    # P0 = P0q + P0t + P0e

    # Km, P0m = boundaryconditions2dof.dirichlet(K, P0, model,
    #                                            displacement_imposed)

    # Ks = sparse.csc_matrix(Km)

    # U = spsolve(Ks, P0m)

    # sNode, sEle, eEle = processing.stress_recovery_simple2(model, U, matDic,
    #                                                        ele.e0)

    # return U, sNode
