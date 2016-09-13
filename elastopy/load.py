import numpy as np
from elastopy.element import Element


def Pb_vector(model, b_force, t=1):
    """Build the load vector due internal body force

    """
    Pb = np.zeros(model.ndof)

    try:
        b_force(1, 1, 1)
    except:
        print("No body force applied!")

    for eid, conn in enumerate(model.CONN):
        element = Element(eid, model)
        pb = element.constructor.body_load(b_force, t)
        Pb[element.id_v] += pb

    return Pb


def Pe_vector(model, material, EPS0, t=1):
    """Build the load vector due initial strain EPS0

    """
    Pe = np.zeros(model.ndof)

    # correct if EPS0 is not specified
    if np.size(EPS0) == 1:
        EPS0 = np.zeros((model.ne, 3))

    for eid, conn in enumerate(model.CONN):
        element = Element(eid, model)
        eps0 = EPS0[eid]
        pe = element.constructor.initial_strain_load(material, eps0, t)
        Pe[element.id_v] += pe

    return Pe
