"""Creates an element object with basic attributes

"""
import numpy as np
from elastopy.elements.quad4 import Quad4


class Element(object):
    """Build an Element object with a constructor for a specific element type

    """
    def __init__(self, eid, model):
        self.type = model.TYPE[eid]
        self.conn = model.CONN[eid]
        self.xyz = model.XYZ[self.conn]
        self.dof = model.DOF[eid]
        self.surf = model.surf_of_ele[eid]

        self.id = np.ix_(self.dof, self.dof)
        self.id_v = self.dof

        # add new elements in here according to its gmsh type!
        if self.type == 3:
            # passing element into quad4 specialization
            self.constructor = Quad4(self)
        else:
            raise Exception('Element not implemented yet!')
