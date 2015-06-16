#!/usr/bin/env python

import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import bgc2
from matplotlib.patches import Circle
from scipy.ndimage.filters import gaussian_filter
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

	#opts, args = get_args(sys.argv[1:])
	#output_file, bgc2_files, ascii_files = parse_args(opts, args)
	ascii_files = [sys.argv[1]]
	bgc2_files = [sys.argv[2]]

	for (ascii_file, bgc2_file) in zip(ascii_files, bgc2_files):
		#  read in halo ID, shape data, etc. from Rockstar output
		ascii_header, ascii_data = read_files(ascii_file, header_line=0, rec_array=True)
		#  read in bgc2 files and make arrays of halo and particle data
		bgc2_header, halos, particles = bgc2.read_bgc2_numpy(bgc2_file)

		#  find array to sort halos by number of particles to work from the biggest down
		halo_indices = np.argsort(halos.npart)
		halo_indices = halo_indices[::-1]
		if start_halo > 0:
			halo_indices = halo_indices[start_halo:]

		#  loop through halos to work on one halo and particle list at a time
		for iteration, halo_index in enumerate(halo_indices):
			#  exit loop if max iteration reached
			if (max_iteration > 0) and (iteration >= max_iteration):
				break

			#  get data for current halo
			halo = np.array(halos[halo_index]).view(np.recarray)
			halo_particles = particles[halo_index]
			ascii_halo = ascii_data[ascii_data.id == halo.id]

			#  check for id match and duplicate halos
			if len(ascii_halo) != 1:
				print "Error:  found %d ASCII halo matches for halo ID %d." % (len(ascii_halo), ascii_halo.id[0])
				continue

			#  skip halos with fewer than specified number of particles
			if (npart_threshold > 0) and (halo.npart < npart_threshold):
				print "Skipping remaining halos with fewer than %d particles." % npart_threshold
				break

			#  convert Mpc to kpc for halo and particle positions
			print "Converting units to kpc..."
			halo.radius = halo.radius * dist_scale
			for pos in halo.x, halo.y, halo.z, halo_particles.x, halo_particles.y, halo_particles.z:
				pos[...] = pos * dist_scale

			#  make particle positions relative to halo center
			print "Making particle positions relative to halo center..."
			print halo_particles[0].x
			for particle_pos, halo_pos in zip([halo_particles.x, halo_particles.y, halo_particles.z], [halo.x, halo.y, halo.z]):
				particle_pos[...] = particle_pos - halo_pos
			print halo_particles[0].x

			#  convert particle cartesian coordinates to (spherical or ellipsoidal) radii
			print "Converting particle positions to spherical radii..."
			r_sphere = np.sqrt((halo_particles.x)**2 + (halo_particles.y)**2 + (halo_particles.z)**2)
			if method == 'sphere':
				r = r_sphere
			elif method == 'ellipsoid':
				print "Rotating eigenvalue matrix of axis ratios..."
				ratios = get_rotated_ratios_matrix(ascii_halo)
				print "Converting particle positions to ellipsoidal radii..."
				r = get_ellipsoid_r(ratios, np.column_stack((halo_particles.x, halo_particles.y, halo_particles.z)))

			#  find half-mass radius
			r_half_mass_sphere, r_half_mass_ell = get_half_mass_r(r, r_sphere)

			#  debug
			print ''
			print 'b_to_a: ', ascii_halo.b_to_a[0]
			print 'c_to_a: ', ascii_halo.c_to_a[0]
			print ''
			print 'r_half_mass_sphere: ', r_half_mass_sphere
			print 'r_half_mass_ell: ', r_half_mass_ell
			print 'r_ell.max(): ', r.max()
			print 'r_sphere.max(): ', r_sphere.max()
			print 'mean r_sphere: ', np.average(r_sphere)
			print 'r_vir (ascii): ', ascii_halo.rvir[0]
			print 'r_vir (bgc2): ', halo.radius
			print ''
			print 'len(particles): ', len(halo_particles)
			print 'npart (bgc2): ', halo.npart
			print 'npart (ascii): ', ascii_halo.num_p[0]
			print ''
			print 'halo id (bgc2): ', halo.id
			print 'halo id (ascii): ', ascii_halo.id[0]
			print ''
			print 'r_max / r_half_mass (sphere): ', r_sphere.max() / r_half_mass_sphere
			print 'r_max / r_half_mass (ell): ', r.max() / r_half_mass_ell
			print ''

			#  save result to array for later output to file
			#  !! todo -- add this

			#  make plots
			if generate_testing_plots or generate_paper_plots:
				make_plots(halo, ascii_halo, halo_particles, ratios, r, r_half_mass_ell)

		#  save results to file
		#  !! todo -- add this

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



def get_rotated_ratios_matrix(ascii_halo):
	#  get rotation angles from A vector (and make them numbers instead of np arrays)
	theta_z = atan(ascii_halo.Ay, ascii_halo.Ax)[0]
	theta_x = atan(ascii_halo.Az, ascii_halo.Ay)[0]

	#  form a diagonal matrix of the inverse-squared axis ratios
	ratios = np.diag([1.0, 1.0 / (ascii_halo.b_to_a)**2, 1.0 / (ascii_halo.c_to_a)**2])

	#  rotate axis ratio matrix about z-axis  -->  X_rot = R^T * X * R
	ratios = z_rotation_matrix(theta_z).T.dot(ratios).dot(z_rotation_matrix(theta_z))

	#  rotate axis ratio matrix about x-axis  -->  X_rot = R^T * X * R
	ratios = x_rotation_matrix(theta_x).T.dot(ratios).dot(x_rotation_matrix(theta_x))

	return ratios



def get_ellipsoid_r(ratios, pos):
	#  convert particle cartesian coordinates to "ellipsoidal" radii
	#  using einstein notation here (sorry), since it's much faster and np.dot/np.tensordot do not behave as expected
	#  we want the element-wise dot product  ->  `pos.dot(ratios).dot(pos.T)`, but for each particle individually
	#  `temp = ratios.dot(pos.T)` works as expected, but `pos.dot(temp)` yields a (N, N) array, but we want (N, )
	#  np.einsum follows einstein notation of summation over repeated indices

	#  temporarily store results for `ratios <dot> pos`
	#  should be equivalent to `ratios.dot(pos.T)`
	tempdot = np.einsum('ij, jk -> ik', ratios, pos.T, order='C')

	#  now find `pos <dot> tempdot` for each row in pos
	#  should be equivalent to `pos[:,0]*tempdot[0,:] + pos[:,1]*tempdot[1,:] + pos[:,2]*tempdot[2,:]`
	r_ell = np.einsum('ij, ji -> i', pos, tempdot, order='C')

	#  take square root and get rid of the extra array dimension to make 1-D
	r_ell = np.sqrt(r_ell)
	r_ell = np.squeeze(r_ell)

	return r_ell



def get_half_mass_r(r, r_real=None):
	#  set r_real to r if using spherical radii and not already done
	if r_real == None:
		r_real = r

	#  find (n/2)th particle(s) and corresponding half-mass radius
	sort_indices = np.argsort(r)
	if len(r) % 2 != 0:
		#  if odd number of particles, simply find the radius of the middle particle
		half_index = sort_indices[len(sort_indices) / 2]
		r_half_mass_sphere = r_real[half_index]
		r_half_mass_ell = r[half_index]
	elif len(r) % 2 == 0:
		#  if even number of particles, find two middle particles and radius halfway between them
		half_index1 = sort_indices[len(sort_indices) / 2 - 1]
		half_index2 = sort_indices[len(sort_indices) / 2]
		r_half_mass_sphere = (r_real[half_index1] + r_real[half_index2]) / 2.
		r_half_mass_ell = (r[half_index1] + r[half_index2]) / 2.

	return r_half_mass_sphere, r_half_mass_ell



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



def make_plots(halo, ascii_halo, particles, ratios, r, r_half_mass):
	if generate_testing_plots:
		ellipse_ring_fig = make_ellipse_ring_plot(particles, r, r_half_mass)
		projectinos_fig  = make_projections_plot(particles, r_half_mass, halo.radius)
	if generate_paper_plots:
		pass

	return 0



def make_ellipse_ring_plot(particles, r, r_half_mass):
	fig = plt.figure(figsize = (9.0, 6.0))
	ax = fig.add_subplot(111, aspect='equal')

	ax = draw_ellipse_ring(ax, particles, r, r_half_mass)

	fig.tight_layout()
	plot_name = "%s%s%s" % (plot_base, 'ellipse_ring', plot_ext)
	plt.savefig(plot_name, bbox_inches='tight')

	return fig



def draw_ellipse_ring(ax, particles, r, r_half_mass):
	print len(particles)
	mask = (np.abs(particles.z) <= particles.z.max() * 0.05)
	particles = particles[mask]
	r = r[mask]
	print len(particles)

	mask = (np.abs(r - r_half_mass) / r_half_mass <= 0.05)
	#mask = (r <= r_half_mass)
	particles = particles[mask]
	print len(particles)

	ax.plot(particles.x, particles.y, linestyle='', marker='.', markersize=2, markeredgecolor='blue')
	ax.add_patch(Circle((0., 0.), r_half_mass, fc="None", ec="black", lw=1))

	return ax




def make_projections_plot(particles, r_half_mass, r_vir):
	fig = plt.figure(figsize = (9.0, 6.0))
	ax = fig.add_subplot(111, aspect='equal')

	plot_lim = np.max((particles.x.max(), particles.y.max()))
	ax = draw_projection(ax, particles.x, particles.y, plot_lim + plot_lim*0.1)
	ax.add_patch(Circle((0., 0.), r_half_mass, fc="None", ec="black", lw=1))
	ax.add_patch(Circle((0., 0.), r_vir, fc="None", ec="black", lw=1))

	fig.tight_layout()
	plot_name = "%s%s%s" % (plot_base, 'projections', plot_ext)
	plt.savefig(plot_name, bbox_inches='tight')

	return 0



def draw_projection(ax, x, y, plot_lim, hx = None, hy = None, r = None):
	limits = [[-plot_lim, plot_lim], [-plot_lim, plot_lim]]
	z, xedges, yedges = np.histogram2d(x, y, bins=npixels, range=limits)
	if log_scale_projections:
		z[z<1.0] = 0.5
		plot_norm = mpl.colors.LogNorm(vmin = 1, vmax = z.max(), clip=True)
	else:
		plot_norm = None
	if extra_smoothing:
		z = gaussian_filter(z, smoothing_radius)
	im = ax.imshow(z.T, extent=(-plot_lim, plot_lim, -plot_lim, plot_lim), \
			interpolation='gaussian', origin='lower', cmap=colormap, norm=plot_norm)
	ax.locator_params(nbins=6)
	if draw_circle and hx != None and hy != None and r != None:
		ax.add_patch(Circle((hx, hy), r, fc="None", ec="black", lw=1))
	if draw_contours:
		x_midpoints = (xedges[:-1] + xedges[1:]) / 2.0
		y_midpoints = (yedges[:-1] + yedges[1:]) / 2.0
		X, Y = np.meshgrid(x_midpoints, y_midpoints)
		ax.contour(X, Y, z.T, 2, colors='black', linewidths=4)
		ax.contour(X, Y, z.T, 2, colors='white', linewidths=2)
	if label_colorbar:
		if log_scale_projections:
			log_format = mpl.ticker.LogFormatterMathtext(10, labelOnlyBase=False)
			ax.cax.colorbar(im, format=log_format)
		else:
			ax.cax.colorbar(im)
	elif 0:
		bar = ax.cax.colorbar(im, ticks=[])
		bar.ax.set_yticklabels([])
		#plt.setp(bar.ax.get_yticklabels(), visible=False)

	return ax



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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	user-settable control parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
start_halo = 0					# first halo to analyze
max_iteration = 1				# number of halos to analyze
#max_iteration = None			# number of halos to analyze
npart_threshold = 100			# minimum number of particles per halo
dist_scale = 1.e3				# convert Mpc to kpc
#method = 'sphere'				# use spherical shells for finding half-mass radius
method = 'ellipsoid'			# use ellipsoidal shells for finding half-mass radius
generate_testing_plots = True	# output plots for testing purposes
generate_paper_plots = False	# output plots for use in a paper
plot_base = 'plots/'			# string to prepend to plot file path
plot_ext = '.eps'				# string to append to plot file path
npixels = 250
smoothing_radius = 0.9
label_colorbar = False
draw_circle = True
draw_contours = False
log_scale_projections = True
extra_smoothing = True
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__ == '__main__':
	main()
