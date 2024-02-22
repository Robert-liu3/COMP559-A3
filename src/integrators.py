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
    # self.x = p1 am I suppose to update self.x?

    #how do I use these two functions?
    # self.ps_verts.update_point_positions(self.x) 
	# self.ps_curves.update_node_positions(self.x)

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
        # TODO: don't just return the inputs
        #linealize and write a linear system
        #use conjugate gradient to solve the linear system
        #make the system the cloff
        mass_matrix = system.M
        ext_forces = system.compute_forces(x0, v0)
        
        # deleted_rows = {}

        # for i in system.pinned:
        #     deleted_rows[i] = {'x0': x0[i], 'v0': v0[i], 'mass_matrix': mass_matrix[i]}
        #     v0 = np.delete(v0, i, axis=0)
        #     mass_matrix = np.delete(mass_matrix, i, axis=0)
        #     ext_forces = np.delete(ext_forces, i, axis=0)
        
        v0_flatten = v0.flatten()
        damping = -system.damping*v0_flatten
        print("this is v0: ", v0_flatten)
        stiffness_matrix = system.compute_stiffness_matrix(x0)

        for i in system.pinned:
            ext_forces[i] = 0
            damping[i] = 0
            v0_flatten[i] = 0

        numerator = mass_matrix@v0_flatten + h*ext_forces.flatten() #what would my fsi be
        print("this is the stiffness matrix: ", stiffness_matrix)
        print("this is the mass matrix: ", mass_matrix)
        print("this is h: ", h)
        print("this is damping: ", damping)
        denominator = mass_matrix - h*damping - (h**2)*stiffness_matrix
        print("this is denominator: ", denominator)

        v1 = np.linalg.solve(denominator, numerator).reshape(v0.shape)

        print("this is v1: ", v1)
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