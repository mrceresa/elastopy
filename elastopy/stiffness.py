"""Assemble the global stiffness matrix

"""
import numpy as np
from elastopy.element import Element


def K_matrix(model, material, t=1):
    """Build the global stiffness matrix

    """
    K = np.zeros((model.ndof, model.ndof))

    for eid, conn in enumerate(model.CONN):
        element = Element(eid, model)
        k = element.constructor.stiffness(material, t)
        K[element.id] += k

    return K
