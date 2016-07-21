import numpy as np


def Pb_vector(model, b_force, t=1):
    """Build the load vector for the internal heat source

    """
    Pb = np.zeros(model.ndof)

    for e, conn in enumerate(model.CONN):
        xyz = model.XYZ[conn]
        dof = model.DOF[e]
        pb = pb_vector(model, xyz, b_force, t)

        id = dof

        Pb[id] += pb

    return Pb


def pb_vector(model, xyz, b_force, t=1):
    """Build the element vector for the internal heat source

    """
    gauss_points = model.chi / np.sqrt(3.0)

    pb = np.zeros(8)
    for gp in gauss_points:
        model.basis_function(gp)
        model.jacobian(xyz)
        dJ = model.detJac
        x1, x2 = model.mapping(xyz)

        pb[0] += model.phi[0]*b_force(x1, x2, t)[0]*dJ
        pb[1] += model.phi[0]*b_force(x1, x2, t)[1]*dJ
        pb[2] += model.phi[1]*b_force(x1, x2, t)[0]*dJ
        pb[3] += model.phi[1]*b_force(x1, x2, t)[1]*dJ
        pb[4] += model.phi[2]*b_force(x1, x2, t)[0]*dJ
        pb[5] += model.phi[2]*b_force(x1, x2, t)[1]*dJ
        pb[6] += model.phi[3]*b_force(x1, x2, t)[0]*dJ
        pb[7] += model.phi[3]*b_force(x1, x2, t)[1]*dJ

    return pb
