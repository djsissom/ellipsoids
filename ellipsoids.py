#!/usr/bin/env python

import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import bgc2
from ipdb import set_trace



def main():
	'''
	ellipsoids.py
	
	This program reads in halo and particle data from bgc2 files and halo shape
	data from the Rockstar's ASCII output.  Halo shape information generated
	from Rockstar is used to fit halos with appropriately-rotated ellipsoidal
	shells, and the half-mass radius is found from the ellipsoidal radius of
	the (n/2)th particle.

	The halo IDs and resulting half-mass radii are saved to disk to be added as
	an additional column in the master halo database.  Optionally, plots are
	generated to demonstrate the ellipsoid halo fitting.

	Requires bgc2.py as a dependency.
	'''

	opts, args = get_args(sys.argv[1:])
	output_file, bgc2_files, ascii_files = parse_args(opts, args)

	for (ascii_file, bgc2_file) in zip(ascii_files, bgc2_files):
		#  read in halo ID, shape data, etc. from Rockstar output
		ascii_header, ascii_data = read_files(ascii_file, header_line=0)
		#  read in bgc2 files and make arrays of halo and particle data
		bgc2_header, halos, particles = bgc2.read_bgc2_numpy(bgc2_file)

		#  find array to sort halos by number of particles to work from the biggest down
		halo_indicies = np.argsort(halos.npart)
		halo_indices = halo_indices[::-1]

		#  loop through halos to work on one halo and particle list at a time
		for iteration, halo_index in enumerate(halo_indices):
			#  exit loop if max iteration reached
			if (max_iteration > 0) and (iteration >= max_iteration):
				break

			#  get data for current halo
			halo = halos[halo_index]
			halo_particles = particles[halo_index]
			ascii_halo = ascii_data[ascii_data.id == halo.id]

			#  check for id match and duplicate halos
			if len(ascii_halo) != 1:
				print "Error:  found %d ASCII halo matches for halo ID %d." % (len(ascii_halo), ascii_halo.id[0])
				continue

			#  skip halos with fewer than specified number of particles
			if (npart_threshold > 0) and (halos.npart < npart_threshold):
				print "Skipping remaining halos with fewer than %d particles." % npart_threshold
				break

			#  convert Mpc to kpc for halo and particle positions
			for pos in halo.x, halo.y, halo.z, halo_particles.x, halo_particles.y, halo_particles.z:
				pos[...] = pos * dist_scale

			#  make particle positions relative to halo center
			for particle_pos, halo_pos in zip([halo_particles.x, halo_particles.y, halo_particles.z], [halo.x, halo.y, halo.z]):
				particle_pos[...] = particle_pos - halo_pos

			#  find half-mass radius
			if method == 'sphere':
				r_half_mass = get_sphere_half_mass_r(halo_particles)
			if method == 'ellipsoid':
				r_half_mass = get_ellipsoid_half_mass_r(ascii_halo, halo_particles)


	#  save results to file

	#  make plots

	print 'Finished.'
	return



def read_files(files, header_line=None, comment_char='#', rec_array=False):
	header = None
	data = None
	if type(files) == str:
		files = [files]

	if header_line != None:
		with open(files[0], 'r') as fd:
			for line in range(header_line):
				fd.readline()
			header = fd.readline()
		if header[0] != comment_char:
			print "Header must start with a '%s'." % comment_char
			sys.exit(4)
		header = header[1:]
		header = header.split()

	for file in files:
		print "Reading file%s..." % (file)
		if data == None:
			if rec_array:
				data = np.genfromtxt(file, dtype=None, comments=comment_char, names=header, deletechars='[]/|')
				data = data.view(np.recarray)
			else:
				data = np.genfromtxt(file, dtype=None, comments=comment_char)
		else:
			if rec_array:
				data = np.append(data, np.genfromtxt(file, dtype=None, comments=comment_char, names=header, deletechars='[]/|'), axis=0)
				data = data.view(np.recarray)
			else:
				data = np.append(data, np.genfromtxt(file, dtype=None, comments=comment_char), axis=0)

	if header_line == None:
		return data
	else:
		return header, data



def get_sphere_half_mass_r(particles):
	#  convert cartesian coordinates to radii
	r = np.sqrt((particles.x)**2 + (particles.y)**2 + (particles.z)**2)

	#  find (n/2)th particle(s) and coresponding half-mass radius
	sort_indices = np.argsort(r)
	if len(r) % 2 != 0:
		#  if odd number of particles, simply find the radius of the middle particle
		half_index = sort_indices[len(sort_indices) / 2]
		r_half_mass = r[half_index]
	elif len(r) % 2 == 0:
		#  if even number of particles, find two middle particles and radius halfway between them
		half_index1 = sort_indices[len(sort_indices) / 2 - 1]
		half_index2 = sort_indices[len(sort_indices) / 2]
		r_half_mass = (r[half_index1] + r[half_index2]) / 2.

	return r_half_mass



def get_ellipsoid_half_mass_r(ascii_halo, particles):
	#  get rotation angles from A vector
	theta_z = atan(ascii_halo.Ay, ascii_halo.Ax)
	theta_x = atan(ascii_halo.Az, ascii_halo.Ay)

	#  form a diagonal matrix of the inverse-squared axis ratios
	ratios = np.diag([1.0, 1.0 / (ascii_halo.b_to_a)**2, 1.0 / (ascii_halo.c_to_a)**2])

	#  rotate axis ratio matrix about z-axis  -->  X_rot = R^T * X * R
	ratios = z_rotation_matrix(theta_z).T.dot(ratios).dot(z_rotation_matrix(theta_z)

	#  rotate axis ratio matrix about x-axis  -->  X_rot = R^T * X * R
	ratios = x_rotation_matrix(theta_x).T.dot(ratios).dot(x_rotation_matrix(theta_x)

	#  convert particle cartesian coordinates to "ellipsoidal" radii
	pos = particles[['x', 'y', 'z']]
	r_ell = pos.T.dot(ratios).dot(pos)			# <-- probably need to use np.tensordot() here

	#  find (n/2)th particle(s) and coresponding half-mass radius
	sort_indices = np.argsort(r)
	if len(r) % 2 != 0:
		#  if odd number of particles, simply find the radius of the middle particle
		half_index = sort_indices[len(sort_indices) / 2]
		r_half_mass = r[half_index]
	elif len(r) % 2 == 0:
		#  if even number of particles, find two middle particles and radius halfway between them
		half_index1 = sort_indices[len(sort_indices) / 2 - 1]
		half_index2 = sort_indices[len(sort_indices) / 2]
		r_half_mass = (r[half_index1] + r[half_index2]) / 2.

	return r_half_mass



def x_rotation_matrix(theta):
	return np.array([[1., 0.,           0.], \
	                 [0., cos(theta), -(sin(theta))], \
	                 [0., sin(theta),   cos(theta)]])



def y_rotation_matrix(theta):
	return np.array([[  cos(theta),  0., sin(theta)], \
	                 [  0.,          1., 0.], \
	                 [-(sin(theta)), 0., cos(theta)]])



def z_rotation_matrix(theta):
	return np.array([[cos(theta), -(sin(theta)), 0.], \
	                 [sin(theta),   cos(theta),  0.], \
	                 [0.,           0.,          1.]])



def sin(theta):
	return np.sin(theta)



def cos(theta):
	return np.cos(theta)



def tan(theta):
	return np.tan(theta)



def atan(x1, x2):
	return np.arctan2(x1, x2)



def add_white_to_colormap(orig_map, num):
	from matplotlib import cm
	temp_cmap = cm.get_cmap(orig_map, num)
	vals = temp_cmap(np.arange(num))
	nfade = num / 7
	vals[:nfade,0] = np.linspace(1., vals[nfade-1,0], nfade)
	vals[:nfade,1] = np.linspace(1., vals[nfade-1,1], nfade)
	vals[:nfade,2] = np.linspace(1., vals[nfade-1,2], nfade)
	newcmap = mpl.colors.LinearSegmentedColormap.from_list("custom_1", vals)
	return newcmap



colormap = add_white_to_colormap('rainbow', 30)
plot_dest_type = 'paper'
if plot_dest_type == 'paper':
	mpl.rcParams['font.family'] = 'serif'
	mpl.rcParams['font.size'] = 16
	mpl.rcParams['axes.linewidth'] = 3
	mpl.rcParams['lines.linewidth'] = 4
	mpl.rcParams['patch.linewidth'] = 4
	mpl.rcParams['xtick.major.width'] = 3
	mpl.rcParams['ytick.major.width'] = 3
	mpl.rcParams['xtick.major.size'] = 8
	mpl.rcParams['ytick.major.size'] = 8
	mpl.rcParams['xtick.minor.width'] = 2
	mpl.rcParams['ytick.minor.width'] = 2
	mpl.rcParams['xtick.minor.size'] = 4
	mpl.rcParams['ytick.minor.size'] = 4



max_iteration					# maximum number of halos to work on
npart_threshold = 100			# minimum number of particles per halo
dist_scale = 1.e3				# convert Mpc to kpc
method = 'sphere'				# use spherical shells for finding half-mass radius
#method = 'ellipsoid'			# use ellipsoidal shells for finding half-mass radius



if __name__ == '__main__':
	main()
