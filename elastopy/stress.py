import numpy as np


def recovery(model, material, U, EPS0=0):
    """Recovery stress at nodes from displacement

    """
    sig = np.zeros(3)
    SIG = np.zeros((model.nn, 3))

    for e, conn in enumerate(model.CONN):
        surf = model.surf_of_ele[e]
        dof = model.DOF[e]

        try:
            E = material.E[surf]
            nu = material.nu[surf]
        except:
            print('Surface {} has no property assigned!'
                  'Default values were used!'.format(surf))
            E = 1.0
            nu = 0.1

        C = c_matrix(E, nu)

        u = U[dof]

        for n, xez in enumerate(model.chi):
            model.basis_function(xez/np.sqrt(3))

            # number of elements sharing a node
            num_ele_shrg = (model.CONN == conn[n]).sum()

            dp_xi = model.dphi_xi
            B = np.array([
                [dp_xi[0, 0], 0, dp_xi[0, 1], 0, dp_xi[0, 2], 0,
                 dp_xi[0, 3], 0],
                [0, dp_xi[1, 0], 0, dp_xi[1, 1], 0, dp_xi[1, 2], 0,
                 dp_xi[1, 3]],
                [dp_xi[1, 0], dp_xi[0, 0], dp_xi[1, 1], dp_xi[0, 1],
                 dp_xi[1, 2], dp_xi[0, 2], dp_xi[1, 3], dp_xi[0, 3]]])

            # sig = [sig_11 sig_22 sig_12] for each n node
            sig = C @ (B @ u - EPS0)

            # dof 1 degree of freedom per node
            d = dof[2*n]/2

            SIG[d, :] += sig / num_ele_shrg

    return SIG


def c_matrix(E, nu):
    """Build the element constitutive matrix

    """
    C = np.zeros((3, 3))
    C[0, 0] = 1.0
    C[1, 1] = 1.0
    C[1, 0] = nu
    C[0, 1] = nu
    C[2, 2] = (1.0 - nu)/2.0
    C = (E/(1.0-nu**2.0))*C

    return C
