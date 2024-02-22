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
    p1 = p0 + h * derivs(p0) 

    return p1

def midpoint(p0, h, derivs):
    # TODO don't just return the inputs
    p1 = p0 + h * derivs(p0+(1/2)*h*derivs(p0))
    
    return p1

def modified_midpoint(p0, h, derivs):
    # TODO don't just return the inputs
    p1 = p0 + h * derivs(p0+(2/3)*h*derivs(p0))

    return p1

def rk4(p0, h, derivs):
    # TODO don't just return the inputs
    k1 = derivs(p0)
    k2 = derivs(p0+(h/2)*k1)
    k3 = derivs(p0+(h/2)*k2)
    k4 = derivs(p0+h*k3)
    p1 = p0 + (h/6)*(k1+2*k2+2*k3+k4)
    return p1

class symplectic_euler:
    name = "Symplectic Euler"
    @staticmethod
    def step(x0,v0, h, system):
        # TODO: don't just return the inputs
        v1 = v0 + h * system.derivs(x0, v0)
        x1 = x0 + h * v1 
        return x1, v1

class backward_euler:
    name = "Backward Euler"
    @staticmethod
    def step(x0, v0, h, system):
        mass_matrix = system.M
        ext_forces = system.compute_forces(x0, v0).flatten()
        v0_og = v0.copy()

        damping = system.damping*np.eye(len(mass_matrix))
        stiffness_matrix = system.compute_stiffness_matrix(x0)

        v0_flatten = v0.flatten()
        b = mass_matrix@v0_flatten + h*ext_forces #what would my fsi be
        A = mass_matrix - h*damping - (h**2)*stiffness_matrix

        pinned = []
        for i in system.pinned:
            pinned.append(2*i)
            pinned.append(2*i+1)
            
        b = np.delete(b, pinned, axis=0)
        # b = np.delete(b, pinned, axis=0)
        A = np.delete(A, pinned, axis=0)
        A = np.delete(A, pinned, axis=1)
        # A = np.delete(A, pinned, axis=0)
        # A = np.delete(A, pinned, axis=1)


        v1 = np.linalg.solve(A, b)
        for i in system.pinned:
            v1 = np.insert(v1, i, 0, axis=0)
            v1 = np.insert(v1, i+1, 0, axis=0)
        v1 = v1.reshape(v0_og.shape)
        
        x1 = x0 + h*v1
        return x1, v1

integration_methods = [
    symplectic_euler(),
    phase_space_integrator("Forward Euler",forward_euler),
    phase_space_integrator("midpoint Method",midpoint),
    phase_space_integrator("RK4",modified_midpoint),
    phase_space_integrator("RK4",rk4),
    backward_euler(),
]