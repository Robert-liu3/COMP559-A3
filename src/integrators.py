# YOUR NAME AND STUDENT NUMBER HERE
# Make any changes you need to this file!

import numpy as np

class phase_space_integrator:
    def __init__(self,name,method):
        self.name = name
        self.method = method
    
    @staticmethod
    def phase_space_unwrapper(p,system):
        x = p[0:system.n]
        v = p[system.n:]
        a = system.derivs(x,v)
        return np.concatenate((v,a))

    def step(self,x0,v0,h,system):
        p0 = np.concatenate((x0,v0))
        f = lambda p : self.phase_space_unwrapper(p,system)
        p1 = self.method(p0,h,f)
        x1 = p1[0:system.n]
        v1 = p1[system.n:]
        return x1,v1

def forward_euler(p0, h, derivs): 
    # TODO don't just return the inputs
    dp = derivs(p0) #is this how your suppose to use derivs? how do you use it?
    p1 = p0 + h * dp 
    self.x = p1 #am I suppose to update self.x?

    #how do I use these two functions?
    # self.ps_verts.update_point_positions(self.x) 
	# self.ps_curves.update_node_positions(self.x)

    return p0

def midpoint(p0, h, derivs):
    # TODO don't just return the inputs
    return p0

def modified_midpoint(p0, h, derivs):
    # TODO don't just return the inputs
    return p0

def rk4(p0, h, derivs):
    # TODO don't just return the inputs
    return p0

class symplectic_euler:
    name = "Symplectic Euler"
    @staticmethod
    def step(x0,v0, h, system):
        # TODO: don't just return the inputs
        return x0, v0

class backward_euler:
    name = "Backward Euler"
    @staticmethod
    def step(x0, v0, h, system):
        # TODO: don't just return the inputs
        return x0, v0

integration_methods = [
    symplectic_euler(),
    phase_space_integrator("Forward Euler",forward_euler),
    phase_space_integrator("midpoint Method",midpoint),
    phase_space_integrator("RK4",modified_midpoint),
    phase_space_integrator("RK4",rk4),
    backward_euler(),
]