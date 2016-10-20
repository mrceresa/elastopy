import numpy as np


def thermal_strain(model, material, T, T0,):
    """Compute thermal strain

    """
    EPS0 = np.zeros((model.ne, 3))
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

    return EPS0
