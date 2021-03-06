"""
Copyright (C) 2021 Gebri Mishtaku

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses.
"""

import time
import os
import json
import numpy as np
from math import inf
from astropy.io import fits
from scipy.signal import fftconvolve
from .utils import *
from .kernel import Kernel



class CenterFinder:


	def __init__(
		self, 
		galaxy_file: str, 
		wtd: bool, 
		params_file: str, 
		save: bool, 
		printout: bool,
		kernel_radius: float = 110., 
		kernel_type: str = 'step', 
		kernel_args: list = [], 
		vote_threshold: float = -inf,
		factorize: bool = False,
		):

		self.kernel_radius = kernel_radius
		self.kernel_type = kernel_type
		self.kernel_args = kernel_args
		self.show_kernel = False
		self.vote_threshold = vote_threshold
		self.weighted = wtd
		self.factorize_randoms = factorize

		self.printout = printout
		self.filename = remove_ext(galaxy_file)
		self.loadname = 'data/' + galaxy_file
		self.save = save
		self.savename = f'out_cf_{self.filename}/'
		if self.save:
			try:
				os.mkdir(self.savename)
			except FileExistsError:
				pass

		# loads galaxy data arrays
		if not self.weighted:
			self.G_ra, self.G_dec, self.G_redshift = load_data(self.loadname)
			self.G_weights = np.ones(len(self.G_ra), dtype=float)
		else:
			self.G_ra, self.G_dec, self.G_redshift, self.G_weights \
				= load_data_weighted(self.loadname)

		# gets cosmology and other hyperparameters
		self.cosmology, self.grid_spacing, self.rounding_precision \
			= load_hyperparameters(params_file)
		# calculates lookup tables for fast conversion from r to z and vice versa
		self.LUT_radii, self.LUT_redshifts = interpolate_r_z(self.G_redshift.min(), 
			self.G_redshift.max(), self.cosmology)

		self.G_radii = self.LUT_radii(self.G_redshift)

		self.randoms_grid: np.ndarray = None
		self.density_grid: np.ndarray = None
		self.density_grid_edges = None
		self.kernel: Kernel = None
		self.background_grid: np.ndarray = None
		self.centers_grid: np.ndarray = None


	def __str__(self):
		return 'CenterFinder object\n'\
				f'Galaxy data file: {self.filename}\n'\
				f'Kernel radius: {self.kernel_radius}\n'\
				f'Vote threshold: {self.vote_threshold}\n'\
				f'RA range: [{self.G_ra.min()}, {self.G_ra.max()}]\n'\
				f'DEC range: [{self.G_dec.min()}, {self.G_dec.max()}]\n'\
				f'Z range: [{self.G_redshift.min()}, {self.G_redshift.max()}]'


	def set_kernel_radius(self, kr: float):
		self.kernel_radius = kr


	def get_kernel_r_idx_units(self):
		return self.kernel.kernel_r_idx_units


	def set_kernel_type(self, kt: str, args = None):
		self.kernel_type = kt
		self.kernel_args = args


	def set_show_kernel(self, sk: bool):
		self.show_kernel = sk


	def set_vote_threshold(self, vt: float):
		self.vote_threshold = vt


	def set_density_grid(self, dg):
		self.density_grid = dg


	def get_density_grid(self):
		return self.density_grid


	def set_density_grid_edges(self, dge):
		self.density_grid_edges = dge


	def get_density_grid_edges(self):
		return self.density_grid_edges


	def set_randoms_grid(self, rg):
		self.randoms_grid = rg


	def get_randoms_grid(self):
		return self.randoms_grid


	def set_kernel(self, k):
		self.kernel = k


	def get_kernel(self):
		return self.kernel


	def set_centers_grid(self, cg):
		self.centers_grid = cg


	def get_centers_grid(self):
		return self.centers_grid


	def set_background_grid(self, bg):
		self.background_grid = bg


	def get_background_grid(self):
		return self.background_grid


	def make_histo_grid(self):
		xyzs = sky2cartesian(self.G_ra, self.G_dec, self.G_redshift, self.LUT_radii) # galaxy x, y and z coords
		self.galaxies_cartesian = np.array(xyzs).T  # each galaxy is represented by (x, y, z)

		# for cf1 and cf0 types in conker
		if not self.density_grid_edges:
			# gets the 3d histogram (density_grid) and the grid bin coordintes in cartesian (grid_edges)
			bin_counts_3d = np.array([np.ceil((xyzs[i].max() - xyzs[i].min()) / self.grid_spacing) 
										for i in range(len(xyzs))], dtype=int)
			# histograms the data points in real space with given weights
			self.density_grid, self.density_grid_edges = np.histogramdd(
				self.galaxies_cartesian, 
				bins=bin_counts_3d, 
				weights=self.G_weights
				)

		# for cfR, cf1 (cf1 w/ x-corr, diff bounds only) types in conker
		else:
			self.density_grid, self.density_grid_edges = np.histogramdd(
				self.galaxies_cartesian,
				bins=self.density_grid_edges,
				weights=self.G_weights
				)

		if self.printout:
			print('Histogramming completed successfully...')
			print('Histogram grid shape:', self.density_grid.shape)

		return np.sum(self.density_grid)


	def make_randoms_grid(self):
		if not self.randoms_grid:
			if self.factorize_randoms:
				self.randoms_grid = self._project_and_sample(
					self.density_grid, 
					self.density_grid_edges
					)
			else:
				self.randoms_grid = self.density_grid


	def make_density_grid(self, dencon, overden):
		# density contrast: subtracts randoms grid (non-conv background)
		if dencon[0]:
			self.density_grid -= self.randoms_grid

			# keep or discard negative valued weights
			# dencon[1] set to True means keep negative weights
			if not dencon[1]:
				if self.printout:
					print('Discarding all negative weights in density grid...')
				self.density_grid[self.density_grid < 0.] = 0.
			if self.printout:
				print('Background subtraction completed successfully...')
		
		# calculates avg density of all nonempty grid cells of the 
		# weighted density field and subtracts it from the density field
		elif overden:
			denavg = np.average(self.density_grid[self.density_grid!=0])
			self.density_grid[self.density_grid!=0] -= denavg
			self.density_grid[self.density_grid!=0] /= denavg
			if self.printout:
				print('Overdensity calculation completed successfully...')

		if self.printout:
			print('Minimum and maximum values of density field grid cells:\n',
				'[{}, {}]'.format(self.density_grid.min(), self.density_grid.max()))		


	def _convolve_density_grid(self):
		"""
		Convolves density grid and sets it as centers grid (signal grid).
		"""
	
		# makes the kernel for scanning over the density grid
		self.kernel = Kernel(self.kernel_type, self.kernel_radius, self.grid_spacing,
			self.printout, self.show_kernel, *self.kernel_args)

		# this scans the kernel over the whole volume of the galaxy density grid
		# calculates the tensor inner product of the two at each step
		# and finally stores this value as the number of voters per that bin in the centers grid
		self.centers_grid = np.round(
			fftconvolve(self.density_grid, self.kernel.get_grid(), mode='same'),
			decimals = self.rounding_precision)
		
		if self.printout:
			print('Convolution of density grid completed successfully...')
			print('Signal grid shape:', self.centers_grid.shape)
			print('Maximum value per single bin in signal grid W:', self.centers_grid.max())
			print('Minimum value per single bin in signal grid W:', self.centers_grid.min())

		# save whole grid without a vote cut
		if self.save:
			np.save(self.savename + f'convolved_density_grid_r_{self.kernel_radius}'
				f'_t_{self.vote_threshold}.npy', self.centers_grid)



	def _convolve_randoms_grid(self):
		"""
		Convolves randoms grid and sets it as background grid.
		"""

		# needed for background in conker
		self.background_grid = np.round(
			fftconvolve(self.randoms_grid, self.kernel.get_grid(), mode='same'),
			decimals = self.rounding_precision)
		
		if self.printout:
			print('Convolution of randoms grid completed successfully...')
			print('Background grid shape:', self.background_grid.shape)
			print('Maximum value per single bin in background grid B:', self.background_grid.max())
			print('Minimum value per single bin in background grid B:', self.background_grid.min())

		# save whole grid without a vote cut
		if self.save:
			np.save(self.savename + f'background_grid_r_{self.kernel_radius}'
				f'_t_{self.vote_threshold}.npy', self.background_grid)



	def _project_and_sample(self, grid: np.ndarray, grid_edges: list) -> np.ndarray:
		
		# TODO: these are unnecessary, remove and reference self's attribute
		bin_centers_edges_xs, bin_centers_edges_ys, bin_centers_edges_zs = \
			np.array([(grid_edges[i][:-1] + grid_edges[i][1:]) / 2 for i in range(len(grid_edges))])

		# TODO: remove, unnecessary
		# if self.save:
		# 	np.save(self.savename + '_xbins.npy', bin_centers_edges_xs)
		# 	np.save(self.savename + '_ybins.npy', bin_centers_edges_ys)
		# 	np.save(self.savename + '_zbins.npy', bin_centers_edges_zs)

		bin_centers_xs, bin_centers_ys, bin_centers_zs = np.array([(x, y, z) 
			for x in bin_centers_edges_xs 
			for y in bin_centers_edges_ys 
			for z in bin_centers_edges_zs]).T
		del bin_centers_edges_xs, bin_centers_edges_ys, bin_centers_edges_zs
		if self.printout:
			print('Number of bin centers in cartesian coordinates:', len(bin_centers_xs))
		"""
		Why can we be sure that it is okay to interpolate the radii 
		and redshift values for these bin centers coordinates?
		Because we know that the range of values of the bin centers 
		is exactly in between the min and the max of the grid bin 
		edges x, y, z. The radii come from the 3d euclidian distance, 
		which preserves this relationship (convex function of x,y,z), 
		and thus it is fine to use the beforehand-calculated interpolation 
		lookup table to find the redshifts from the radii.
		"""
		bin_centers_ra, bin_centers_dec, _, bin_centers_radii = \
			cartesian2sky(bin_centers_xs, 
							bin_centers_ys, 
							bin_centers_zs, 
							self.LUT_redshifts, 
							self.G_ra.min(), 
							self.G_ra.max())
		del bin_centers_xs, bin_centers_ys, bin_centers_zs
		if self.printout:
			print('Number of bin centers in sky coordinates:', len(bin_centers_ra))

		# sum total of weights
		N_tot = np.sum(grid)
		if self.printout:
			print('Sum total of weights:', N_tot)

		start = time.time()
		# get volume adjustment grid, the differentials in sky coordinate dimensions 
		# and the number of bins in each dimension
		vol_adjust_ratio_grid, d_r, d_alpha, d_delta, N_bins_r, N_bins_alpha, N_bins_delta = \
			self._volume_adjustment(bin_centers_radii, bin_centers_ra, bin_centers_dec, grid.shape)
		end = time.time()
		print(f'\nvolume_adjustment took time: {end-start} seconds\n')

		start = time.time()
		# alpha-delta and z counts
		N_bins_x, N_bins_y, N_bins_z = grid.shape[0], grid.shape[1], grid.shape[2]
		sky_coords_grid_shape = (N_bins_x, N_bins_y, N_bins_z, 3)  # need to store a triple at each grid bin
		sky_coords_grid = np.array(list(zip(bin_centers_ra, bin_centers_dec, 
			bin_centers_radii))).reshape(sky_coords_grid_shape)
		if self.printout:
			print('Shape of grid containing sky coordinates of observed grid\'s bin centers:', 
				sky_coords_grid.shape)
		end = time.time()
		print(f'\ncreating sky_coords_grid took time: {end-start} seconds\n')

		start = time.time()
		# getting some variables ready for the projection step
		alpha_min = bin_centers_ra.min()
		d_alpha = np.rad2deg(d_alpha)
		delta_min = bin_centers_dec.min()
		d_delta = np.rad2deg(d_delta)
		r_min = bin_centers_radii.min()
		del bin_centers_ra, bin_centers_dec, bin_centers_radii

		# vectorial computation of the sky indices
		sky_coords_grid[:, :, :, 0] = (sky_coords_grid[:, :, :, 0] - alpha_min) // d_alpha
		sky_coords_grid[:, :, :, 1] = (sky_coords_grid[:, :, :, 1] - delta_min) // d_delta
		sky_coords_grid[:, :, :, 2] = (sky_coords_grid[:, :, :, 2] - r_min) // d_r
		sky_coords_grid = sky_coords_grid.astype(int)

		# TODO: the condition here should be >= rather than ==
		# the following fixes any indices that lie beyond the outer 
		# walls of the sky grid by pulling them to the wall
		sky_coords_grid[:, :, :, 0][sky_coords_grid[:, :, :, 0] == N_bins_alpha] = N_bins_alpha - 1
		sky_coords_grid[:, :, :, 1][sky_coords_grid[:, :, :, 1] == N_bins_delta] = N_bins_delta - 1
		sky_coords_grid[:, :, :, 2][sky_coords_grid[:, :, :, 2] == N_bins_r] = N_bins_r - 1

		end = time.time()
		print(f'\nmodifying sky_coords_grid took time: {end-start} seconds\n')

		start = time.time()
		alpha_delta_grid, r_grid = self._alpha_delta_r_projections_from_grid(grid, 
			N_bins_x, N_bins_y, N_bins_z, sky_coords_grid, N_bins_alpha, N_bins_delta, N_bins_r)
		end = time.time()
		print(f'\nalpha_delta_r_projections took time: {end-start} seconds\n')
		if self.printout:
			print('Shape of alpha-delta grid:', alpha_delta_grid.shape)
			print('Shape of r grid:', r_grid.shape)
			print('Maximum value per single bin in alpha-delta grid:', alpha_delta_grid.max())
			print('Minimum value per single bin in alpha-delta grid:', alpha_delta_grid.min())
			print('Maximum value per single bin in r grid:', r_grid.max())
			print('Minimum value per single bin in r grid:', r_grid.min())
			print('N_tot_observed = N_tot_alpha_delta = N_tot_r:', 
				N_tot == np.sum(alpha_delta_grid) == np.sum(r_grid))

		start = time.time()
		i = np.arange(N_bins_x)[:,None,None]
		j = np.arange(N_bins_y)[None,:,None]
		k = np.arange(N_bins_z)[None,None,:]
		randoms_grid = alpha_delta_grid[sky_coords_grid[i,j,k,0], sky_coords_grid[i,j,k,1]] \
						* r_grid[sky_coords_grid[i,j,k,2]]
		end = time.time()
		print(f'\nrandoms_grid took time: {end-start} seconds\n')

		randoms_grid /= N_tot  # normalization
		randoms_grid *= vol_adjust_ratio_grid  # volume adjustment
		if self.printout:
			print('Randoms grid shape:', randoms_grid.shape)
			print('Maximum value in randoms grid bin:', randoms_grid.max())
			print('Minimum value in randoms grid bin:', randoms_grid.min())

		if self.save:
			np.save(self.savename + "_randoms_grid.npy", randoms_grid)

		return randoms_grid



	def _volume_adjustment(self, bin_centers_radii: np.array, bin_centers_ra: np.array, 
		bin_centers_dec: np.array, observed_grid_shape: tuple):
		
		# radius
		mid_r = (bin_centers_radii.max() + bin_centers_radii.min()) / 2
		delta_r = bin_centers_radii.max() - bin_centers_radii.min()
		N_bins_r = int(np.ceil(delta_r / self.grid_spacing))
		d_r = self.grid_spacing
		r_sqr = bin_centers_radii ** 2

		# alpha
		delta_alpha = np.deg2rad(bin_centers_ra.max() - bin_centers_ra.min())
		N_bins_alpha = int(np.ceil((delta_alpha * mid_r) / self.grid_spacing))
		d_alpha = delta_alpha / N_bins_alpha

		# delta
		delta_delta = np.deg2rad(bin_centers_dec.max() - bin_centers_dec.min())
		N_bins_delta = int(np.ceil((delta_delta * mid_r) / self.grid_spacing))
		d_delta = delta_delta / N_bins_delta
		cos_delta = np.cos(np.deg2rad(bin_centers_dec))

		# angular volume differential
		dV_ang = d_alpha * cos_delta * d_delta * r_sqr * d_r
		# euclidean volume differential
		dV_xyz = self.grid_spacing ** 3
		# volume adjustment ratio grid; contains the volume adjustment ratio 
		# per each bin in the expected grid
		vol_adjust_ratio_grid = (dV_xyz / dV_ang).reshape(observed_grid_shape)

		if self.printout:
			print('Number of bins in r:', N_bins_r)
			print('Number of bins in alpha:', N_bins_alpha)
			print('Number of bins in delta:', N_bins_delta)
			print('Volume adjustment ratio grid shape:', vol_adjust_ratio_grid.shape)

		return vol_adjust_ratio_grid, d_r, d_alpha, d_delta, N_bins_r, N_bins_alpha, N_bins_delta



	def _alpha_delta_r_projections_from_grid(self, 
		density_grid: np.ndarray, 
		N_bins_x: int, N_bins_y: int, N_bins_z: int, 
		sky_coords_grid: np.ndarray, 
		N_bins_alpha: int, N_bins_delta: int, N_bins_r: int) \
		-> (np.ndarray, np.ndarray):

		alpha_delta_grid, _, _ = np.histogram2d(sky_coords_grid[:,:,:,0].ravel(), 
												sky_coords_grid[:,:,:,1].ravel(), 
												bins=(N_bins_alpha, N_bins_delta), 
												weights=density_grid.ravel())

		r_grid, _ = np.histogram(sky_coords_grid[:,:,:,2].ravel(), 
								bins=N_bins_r, 
								weights=density_grid.ravel())

		return alpha_delta_grid, r_grid



	def cleanup(self):
		del self.centers_grid
		del self.background_grid
		del self.kernel_radius


	def make_grids(self, 
		dencon: (bool,bool) = (False, False), 
		overden: bool = False
		):
		"""Creates density and randoms grids. """
		self.make_histo_grid()
		self.make_randoms_grid()
		self.make_density_grid(dencon, overden)


	def make_convolved_grids(self):
		"""Creates convolved density (signal) and background grids. """
		self._convolve_density_grid()
		self._convolve_randoms_grid()


	def find_centers(self, dencon: bool, overden: bool, cleanup: bool = True):
		"""
		Identifies BAO centers by applying the voting procedure.
		Applies vote threshold and saves the found centers list as .fits catalog.
		"""

		if self.printout:
			print(self)

		self._make_grids(dencon=dencon, overden=overden)
		self._convolve_density_grid()
		
		# TODO: clean this up
		# self.centers_grid[self.centers_grid < self.vote_threshold] = 0
		self.centers_indices = np.asarray(self.centers_grid >= self.vote_threshold).nonzero()
		self.C_weights = self.centers_grid[self.centers_indices]
		if self.printout:
			precut = len(self.centers_grid[self.centers_grid!=0])
			postcut = len(self.C_weights)
			print('Number of found centers before vote cut:', precut)
			print('Number of found centers after vote cut:', postcut)
		if cleanup:
			delattr(self, 'centers_grid')

		# calculates center coords to be exactly at the center of their respective bins
		centers_bin_coords = np.array([(self.density_grid_edges[i][:-1] \
										+ self.density_grid_edges[i][1:]) / 2 
			for i in range(len(self.density_grid_edges))])
		if cleanup:
			delattr(self, 'density_grid_edges')
		C_xyzs = np.array([centers_bin_coords[i][self.centers_indices[i]] \
							for i in range(len(self.centers_indices))])
		self.C_ra, self.C_dec, self.C_redshift, _ = cartesian2sky(*C_xyzs, self.LUT_redshifts, 
																self.G_ra.min(), self.G_ra.max())
		
		# outputs centers catalog in skycoords+weights to out folder
		if self.kernel_type=='ball':
			savename = self.savename + f'ball_found_centers_r_{self.kernel_radius}_cut_{self.vote_threshold}.fits'
		else:
			savename = self.savename + f'found_centers_r_{self.kernel_radius}_cut_{self.vote_threshold}.fits'
		save_data_weighted(savename, self.C_ra, self.C_dec, self.C_redshift, self.C_weights)












