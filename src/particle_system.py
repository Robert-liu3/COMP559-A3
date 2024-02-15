# YOUR NAME AND STUDENT NUMBER HERE
# Make any changes you need to this file!

import numpy as np
import json
import polyscope as ps

class particle_system:
	def __init__(self, file_name):
		self.name = file_name
		data = json.load(open("..\\data\\" + file_name, 'r'))
		self.x = np.array(data['vertices'],dtype='f4')
		self.edges = np.array(data['edges'],dtype='i4') 
		self.pinned = np.array(data['pinned'])	# list of pinned particles
		self.n = len(self.x)
		if 'velocities' in data:
			self.v = np.array(data['velocities'],dtype='f4')
		else:
			self.v = np.zeros_like(self.x)
		
		self.x0 = self.x.copy() # save initial state
		self.v0 = self.v.copy()	# save initial state
		# ideally should load mass from file instead of assuming identity
		self.mass = np.ones_like(self.x) # per particle mass
		self.M = np.eye( self.n * 2	) # mass matrix
		
		# edges define springs, compute rest lengths (otherwise could load from file)
		self.rest_length = np.linalg.norm(self.x[self.edges[:,0]] - self.x[self.edges[:,1]], axis=1)

		self.ps_verts = ps.register_point_cloud("verts", self.x, radius=0.01)
		self.ps_curves = ps.register_curve_network("edges", self.x, self.edges)
		self.colors = np.repeat([[0.5,0.5,1]], self.n, axis=0)
		self.colors[self.pinned] = [1,0,0]
		self.ps_verts.add_color_quantity("pinned", self.colors,enabled=True)
		
		self.gravity = -9.8
		self.h = 0.01
		self.stiffness = 100
		self.damping = 0
		self.elapsed = 0	
		
		# TODO: add any other setup code you need here!
		
	def compute_stiffness_matrix(self,x):
		# TODO: you probably want to implement this (unless you are using CG to solve Ax=b without assembly)
		pass

	def compute_damping_matrix(self,x):
		# TODO: you probably want to implement this (unless you are using CG to solve Ax=b without assembly)
		pass

	def compute_forces(self, x, v):
		# TODO: compute and return forces (gravity, spring, damping)
		a_minus_b = x[self.edges[:,0]] - x[self.edges[:,1]]
		a_minus_b_normal = np.linalg.norm(a_minus_b)

		#is stiffness the material stiffness?
		spring_forces = self.stiffness * (a_minus_b_normal - self.rest_length) * a_minus_b / a_minus_b_normal

		#i need to calculate the relative velocity of both particles?
		damping_forces = -self.damping * (v[self.edges[:,0]] - v[self.edges[:,1]])

		gravity_forces = self.gravity * self.mass
		gravity_force_vector = np.array([0, gravity_forces])


		return spring_forces + gravity_force_vector + damping_forces
	
	def derivs(self, x, v):
		f = self.compute_forces(x, v)
		f[self.pinned] = np.zeros((len(self.pinned), 2))
		return np.divide(f, self.mass)

	def advance_time( self, method ):	
		self.x, self.v = method.step(self.x, self.v, self.h, self)        
		self.elapsed += self.h
		self.ps_verts.update_point_positions(self.x)
		self.ps_curves.update_node_positions(self.x)

	def reset(self):
		self.x = self.x0.copy()
		self.v = self.v0.copy()
		self.elapsed = 0