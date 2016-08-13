from diffuspy import steadystate
import numpy as np


def temperature(model, material, temp_bc, flux_bc, b_heat, T0, t):
    """ build the initial strain due change in temeperature


    Return:
    EPS0: [eps0_11, eps0_22, eps0_12] for each element

    """
    EPS0 = np.zeros((model.ne, 3))
    np.set_printoptions(precision=3)
    T = steadystate.solver(model, material, b_heat, flux_bc,
                           temp_bc, t)

    if np.size(T0) == 1:
        T0 = np.zeros(model.nn)

    for e, conn in enumerate(model.CONN):
        surf = model.surf_of_ele[e]
        # 1 dof per node
        dof = np.int_(model.DOF[e][::2]/2)

        try:
            alpha = material.thrmexp[surf]
        except:
            print('Surface {} has no property assigned!'
                  'Default values were used!'.format(surf))
            alpha = 1e-6

        # averageof temperature on each element node
        dT = np.average(T[dof] - T0[dof])

        EPS0[e, :] = dT * alpha * np.array([1, 1, 0])

    return EPS0, T
