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

		m = self.n
		stiffness_matrix = np.zeros((2*m, 2*m))

		for i in range(len(self.edges)):
			if i in self.pinned:
				continue
			k = self.stiffness
			lo = self.rest_length[i]
			fa_da = self.gradientFunc(x[self.edges[i,0]], x[self.edges[i,1]], k, lo)
			fa_db = -fa_da #0,0+1 1,1+1
			fb_da = -fa_da #1,1+1 0,1+1
			fb_db = fa_da #1,1+1 1,1+1

			index_a = self.edges[i,0]
			index_a_first = 2*index_a
			index_a_second = 2*(index_a + 1)

			index_b = self.edges[i,1]
			index_b_first = 2*index_b
			index_b_second = 2*(index_b + 1)

			stiffness_matrix[index_a_first:index_a_second, index_a_first:index_a_second] += fa_da
			stiffness_matrix[index_a_first:index_a_second, index_b_first:index_b_second] += fa_db
			stiffness_matrix[index_b_first:index_b_second, index_a_first:index_a_second] += fb_da
			stiffness_matrix[index_b_first:index_b_second, index_b_first:index_b_second] += fb_db

		return stiffness_matrix
	
	def gradientFunc(self, A, B, k, l):
		assert isinstance(A, np.ndarray)
		dim = A.shape
		assert len(dim) == 1
		A_rows = dim[0]
		assert isinstance(B, np.ndarray)
		dim = B.shape
		assert len(dim) == 1
		B_rows = dim[0]
		if isinstance(k, np.ndarray):
			dim = k.shape
			assert dim == (1, )
		if isinstance(l, np.ndarray):
			dim = l.shape
			assert dim == (1, )
		assert B_rows == A_rows

		t_0 = (A - B)
		t_1 = np.linalg.norm(t_0)
		t_2 = (k * (t_1 - l))
		T_3 = np.outer(t_0, t_0)
		t_4 = (t_2 / t_1)
		gradient = -((((k / (t_1 ** 2)) * T_3) - ((t_2 / (t_1 ** 3)) * T_3)) + (t_4 * np.eye(B_rows, A_rows)))

		return gradient

	def compute_damping_matrix(self,x):
		# TODO: you probably want to implement this (unless you are using CG to solve Ax=b without assembly)
		# Initialize an empty matrix
		damping_matrix = np.zeros((self.n, self.n))

		# For each spring in the system
		for i in range(len(self.edges)):
			# Get the indices of the particles at the ends of the spring
			a, b = self.edges[i]

			# Compute the damping coefficient for this spring
			d = self.damping[i]

			# Update the damping matrix
			damping_matrix[a, a] += d
			damping_matrix[b, b] += d
			damping_matrix[a, b] -= d
			damping_matrix[b, a] -= d

		return damping_matrix
		pass

	def compute_forces(self, x, v):
		# TODO: compute and return forces (gravity, spring, damping)
		force_particle = np.zeros((len(self.x), 2))

		for i in range(len(self.edges)):
			# a_minus_b = x[self.edges[i,0]] - x[self.edges[i,1]]
			# a_minus_b_normal = np.linalg.norm(a_minus_b)
			# displacement = a_minus_b_normal - self.rest_length[i]
			# spring_forces = self.stiffness * displacement * a_minus_b / a_minus_b_normal

			spring_forces = -self.stiffness * (np.linalg.norm(x[self.edges[i,0]] - x[self.edges[i,1]]) - self.rest_length[i]) * (x[self.edges[i,0]] - x[self.edges[i,1]]) / np.linalg.norm(x[self.edges[i,0]] - x[self.edges[i,1]])

			force_particle[self.edges[i,0]] += spring_forces
			# Initialize an empty matrix

			force_particle[self.edges[i,1]] -= spring_forces

		damping_forces = -self.damping*v
		gravity_force_vector = np.zeros((len(self.x), 2))
		gravity_force_vector[:, 1] = self.gravity * self.mass[:,1]
		force_particle += gravity_force_vector + damping_forces



		# a_minus_b = x[self.edges[:,0]] - x[self.edges[:,1]]
		# # print("rest length: ", self.rest_length)
		# a_minus_b_normal = np.linalg.norm(a_minus_b, axis=1)
		# # print("a_minus_b_normal: ", a_minus_b_normal)
		# a_minus_b_normal_reshaped = a_minus_b_normal.reshape(-1, 1)

		# displacement = (a_minus_b_normal - self.rest_length).reshape(-1, 1)
		# spring_forces = -self.stiffness * displacement * a_minus_b / a_minus_b_normal_reshaped
		# print("spring forces: ", spring_forces)

		# #i need to calculate the relative velocity of both particles?
		# damping_forces = -self.damping * v

		# print("damping forces: ", damping_forces)

		# # print("this is the mass: ", self.mass)
		# gravity_forces = self.gravity * self.mass
		# print("gravity forces: ", gravity_forces)
		# gravity_force_vector = np.zeros((len(self.x), 2))
		# gravity_force_vector[:, 1] = gravity_forces[:,1]
		# # print("gravity force vector: ", gravity_force_vector)

		# total_edge_forces = np.zeros_like(gravity_force_vector)
		# np.add.at(total_edge_forces, self.edges[:, 0], spring_forces + damping_forces)
		# np.subtract.at(total_edge_forces, self.edges[:, 1], spring_forces + damping_forces)

		# print("damping matrix", self.compute_stiffness_matrix(x))

		# force_particle = total_edge_forces + gravity_force_vector

		print("force particle: ", force_particle)

		return force_particle
	
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