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

import sys, time
import numpy as np
import matplotlib.pyplot as plt
from .centerfinder import CenterFinder
from .utils import remove_ext



class Correlator:
	

	def __init__(
		self, 
		order: int,
		file1: str, 
		file0: str = None,
		fileR: str = None,
		file_gridR: str = None,
		nondiag: bool = False,
		params_file: str = 'params.json',
		save: bool = False,
		save_randoms = False,
		savename: str = None,
		printout: bool = False
		):

		try:
			assert order >= 2
		except AssertionError:
			print('AssertionError: '\
				f'Correlation order has to be >= 2, was given {order}')
			sys.exit(1)

		try:
			assert not ((order == 2) and nondiag)
		except AssertionError:
			print('AssertionError: '\
				f'Can\'t have a nondiagonal correlation of order 2')
			sys.exit(1)

		self.order = order
		self.file1 = file1
		# auto-correlation on default
		self.type = 'auto'
		# cross-correlation if file0 provided
		if file0:
			self.type = 'cross'
			self.file0 = file0
		# makes background from randoms catalog
		self.fileR = fileR
		# loads randoms grid from file
		self.file_gridR = file_gridR

		self.func_type = 'diag'
		if nondiag:
			self.func_type = 'nondiag'

		# default calibration reference
		# this just serves as a default
		# and may be overridden in load_calib()
		self.calibfile = 'calib_ref_random0_DR12v5_CMASS_South_gs_5.npy'
		self.calib: float = None

		# # correlation order has to match radii given
		# assert order == len(r1_list) + 1
		# self.r1_list = r1_list

		self.params_file = params_file
		self.save = save
		self.save_randoms = save_randoms
		self.savename = savename
		self.printout = printout

		self.cfR: CenterFinder = None
		self.cf0: CenterFinder = None
		self.cf1: CenterFinder = None
		self.dencon_args0: tuple = None
		self.dencon_args1: tuple = None
		self.density_grid: np.ndarray = None
		self.randoms_grid: np.ndarray = None
		self.W0: np.ndarray = None
		self.B0: np.ndarray = None
		self.W1_prod: np.ndarray = None
		self.B1_prod: np.ndarray = None
		self.corrfunc: dict = None
		self.ND: float = None
		self.NR: float = None


	def set_cfR(self, cfR):
		self.cfR = cfR


	def get_cfR(self):
		return self.cfR


	def set_cf0(self, cf0):
		self.cf0 = cf0


	def get_cf0(self):
		return self.cf0


	def set_cf1(self, cf1):
		self.cf1 = cf1


	def get_cf1(self):
		return cf1


	def make_cfR(self, dge):
		if self.file_gridR:
			self.cfR.set_randoms_grid(np.load('data/'+self.file_gridR))
			self.randoms_grid = self.cfR.get_randoms_grid()
		else:
			# cfR only needs to build the randoms grid
			self.cfR.set_density_grid_edges(dge) # from cf0
			self.NR = self.cfR.make_histo_grid()
			self.cfR.make_randoms_grid()
			# normalizes the randoms grid
			# sets randoms grid for use with cf0 and cf1
			data2rand_ratio = self.ND / self.NR
			self.randoms_grid = self.cfR.get_randoms_grid()\
								* data2rand_ratio
			if self.save_randoms:
				np.save(self.savename + 'gridR_{}_gs_{}.npy'\
					.format(remove_ext(self.fileR), self.cf0.grid_spacing), 
					self.randoms_grid)


	def _customize_cf0(self, args):
		# defaults kernel_radius to 1/2 grid_spacing for 0th centerfinder
		self.cf0.set_kernel_radius(self.cf0.grid_spacing / 2)
		if args.kernel_radius0:
			self.cf0.set_kernel_radius(args.kernel_radius0)
		if args.show_kernel0:
			self.cf0.set_show_kernel(args.show_kernel0)
		if args.step_kernel0:
			self.cf0.set_kernel_type('step', args.step_kernel0)
		elif args.gaussian_kernel0:
			self.cf0.set_kernel_type('gaussian', args.gaussian_kernel0)
		elif args.wavelet_kernel0:
			self.cf0.set_kernel_type('wavelet', args.wavelet_kernel0)
		elif args.custom_kernel0:
			self.cf0.set_kernel_type('custom', args.custom_kernel0)
		if args.vote_threshold0:
			self.cf0.set_vote_threshold(args.vote_threshold0)
		if args.density_contrast0:
			do_dencon0 = True
			if len(args.density_contrast0)==0:
				keep_neg_wts0 = True
			else:
				keep_neg_wts0 = False
		else:
			do_dencon0 = True
			keep_neg_wts0 = True
		self.dencon_args0 = (do_dencon0, keep_neg_wts0)


	def prep_cf0(self, args):
		# customizes the cf object and the bckgrnd subtraction
		self._customize_cf0(args)
		self.ND = self.cf0.make_histo_grid()


	def make_cf0(self, args):
		self.cf0.set_randoms_grid(self.randoms_grid)
		self.cf0.make_density_grid(dencon=self.dencon_args0,
			overden=args.overdensity0)

		# if auto-correlation, keeps density so it calculates only once
		if self.type == 'auto':
			self.density_grid = self.cf0.get_density_grid()
		# there's no need for convolving 
		# if the whole kernel is just one cell
		if self.cf0.kernel_radius == self.cf0.grid_spacing/2:
			self.W0 = self.cf0.get_density_grid()
			self.B0 = self.cf0.get_randoms_grid()
		# runs the centerfinding algorithm
		else:
			self.cf0.make_convolved_grids()
			self.W0 = self.cf0.get_centers_grid()
			self.B0 = self.cf0.get_background_grid()


	def _customize_cf1(self, args):
		if args.kernel_radius1:
			self.cf1.set_kernel_radius(args.kernel_radius1)
		if args.show_kernel1:
			self.cf1.set_show_kernel(args.show_kernel1)
		if args.step_kernel1:
			self.cf1.set_kernel_type('step', args.step_kernel1)
		elif args.gaussian_kernel1:
			self.cf1.set_kernel_type('gaussian', args.gaussian_kernel1)
		elif args.wavelet_kernel1:
			self.cf1.set_kernel_type('wavelet', args.wavelet_kernel1)
		elif args.custom_kernel1:
			self.cf1.set_kernel_type('custom', args.custom_kernel1)
		if args.vote_threshold1:
			self.cf1.set_vote_threshold(args.vote_threshold1)
		if args.density_contrast1:
			do_dencon1 = True
			if len(args.density_contrast1)==0:
				keep_neg_wts1 = True
			else:
				keep_neg_wts1 = False
		else:
			do_dencon1 = True
			keep_neg_wts1 = True
		self.dencon_args1 = (do_dencon1, keep_neg_wts1)


	def make_cf1(self, args):
		# customizes the cf object and the bckgrnd subtraction
		self._customize_cf1(args)
		self.cf1.set_randoms_grid(self.randoms_grid) # from cfR
		# doesn't recalculate density grid if auto-cor
		if self.type == 'auto':
			self.cf1.set_density_grid(self.density_grid)
		# runs the centerfinding algorithm again if x-cor
		# note cf0 boundaries override cf1 if different in x-corr
		else:
			self.cf1.set_density_grid_edges(self.cf0.get_density_grid_edges())
			self.cf1.make_grids(dencon=self.dencon_args1, 
				overden=args.overdensity1)


	@staticmethod
	def _reduce_mult_till_order(cf1: CenterFinder, order: int) \
		-> (np.ndarray, np.ndarray):
		"""DEPRECATED METHOD
		Reduction that maps high order corrfunc over same object.
		Note this acts like a reduction but the args list is
		created dynamically, i.e. W1, W1-1, W1-2 ... W1-(n-2)
		Strictly for use on Correlator instance's cf1 object.
		"""
		W1_prod = cf1.get_centers_grid()
		B1_prod = cf1.get_background_grid()
		for i in range(1, order-1):
			W1_prod *= (W1_prod - i)
			B1_prod *= (B1_prod - i)

		return W1_prod, B1_prod


	def _make_W1_B1(self):
		self.cf1.make_convolved_grids()
		W1 = self.cf1.get_centers_grid()
		B1 = self.cf1.get_background_grid()

		return W1, B1


	def _correlate_diag(self):
		W1, B1 = self._make_W1_B1()
		W = self.W0 * W1 ** (self.order - 1)
		B = self.B0 * B1 ** (self.order - 1)
		if self.save:
			r1, r0 = self.cf1.kernel_radius, self.cf0.kernel_radius
			np.save(self.savename + '{}pcf_W_r1_{}_r0_{}.npy'\
					.format(self.order, r1, r0), W)
			np.save(self.savename + '{}pcf_B_r1_{}_r0_{}.npy'\
					.format(self.order, r1, r0), B)
		self.cf1.cleanup()

		W = np.sum(W)
		B = np.sum(B)

		return W, B


	@staticmethod
	def _multistep_product_from_dict(dict, steps_list):
		product = 1
		for s in steps_list:
			product *= dict[s]
		return product


	@staticmethod
	def _steps_list_to_idx_tuple(steps, steps_list):
		return tuple((np.argwhere(steps==s)[0,0] for s in steps_list))


	def _recurse_dimensions_correlate(self, 
		steps,
		running_steps_list,
		maxdims,
		currdim, 
		loopidx_start):

		# ready to fill cell in W_hyperarr and B_hyperarr
		if currdim == maxdims:
			steps_idx_list = Correlator._steps_list_to_idx_tuple(steps, running_steps_list)
			print('Current separations and their idxs:\t{}\t{}'\
				.format(running_steps_list, steps_idx_list))
			W = np.sum(self.W0 * 
				Correlator._multistep_product_from_dict(
					self.W1_dict, running_steps_list))
			B = np.sum(self.B0 * 
				Correlator._multistep_product_from_dict(
					self.B1_dict, running_steps_list))
			self.corrfunc_hyperarr[steps_idx_list] = W / B

		# recurses to further dimensions
		else:
			for i, s in enumerate(steps[loopidx_start:]):
				new_running_steps_list = running_steps_list + [s]
				self._recurse_dimensions_correlate(
					steps,
					new_running_steps_list,
					maxdims,
					currdim+1,
					i)


	def save_single(self):
		separation = np.array(list(self.corrfunc.keys()))
		correlation = np.array(list(self.corrfunc.values()))
		if self.printout:
			print('Separation value:\n', separation)
			print('Correlation value:\n', correlation)
		np.save(self.savename + '{}pcf_sep_{}.npy'\
			.format(self.order, separation[0]), separation)
		np.save(self.savename + '{}pcf_corr_{}.npy'\
			.format(self.order, separation[0]), correlation)


	def single_correlate(self):
		s = self.cf1.kernel_radius
		self.corrfunc = {s: None}
		W, B = self._correlate_diag()
		self.corrfunc[s - self.calib] = W / B
		self.save_single()


	def _save_scan(self, scan_args):
		losep, hisep = scan_args
		separation = np.array(list(self.corrfunc.keys()))
		correlation = np.array(list(self.corrfunc.values()))
		if self.printout:
			print('Separation array:\n', separation)
			print('Correlation array:\n', correlation)
		np.save(self.savename + '{}pcf_sep_range_{}_{}.npy'\
			.format(self.order, losep, hisep), separation)
		np.save(self.savename + '{}pcf_corr_range_{}_{}.npy'\
			.format(self.order, losep, hisep), correlation)


	def _save_plot(self, scan_args):
		losep, hisep = scan_args
		separation = list(self.corrfunc.keys())
		correlation = list(self.corrfunc.values())
		plt.plot(separation, correlation, 
			label='conker r0={} gs={}'\
			.format(self.cf0.kernel_radius, self.cf0.grid_spacing))

		plt.title('ConKer: {} {} {}pcf'\
			.format(
				remove_ext(self.file1), 
				'weighted' if self.cf1.weighted else 'unweighted',
				self.order
				)
			)
		plt.xlabel(r'$s$ $[h^{-1}Mpc]$')
		plt.ylabel(r'$\xi(s)$')
		plt.grid(linestyle=':')
		plt.legend()

		plt.savefig('plots/' + '{}pcf_scan_{}_{}_{}.png'\
			.format(self.order, losep, hisep, remove_ext(self.file1)), 
			dpi=300)


	def scan_correlate_diag(self, scan_args):
		start, end = scan_args
		step = self.cf1.grid_spacing
		steps = np.arange(start, end+step, step)
		self.corrfunc = {}

		for s in steps:
			self.cf1.set_kernel_radius(s)
			W, B = self._correlate_diag()
			self.corrfunc[s - self.calib] = W / B
			np.save(self.savename + '{}pcf_W_r1_{}.npy'\
				.format(self.order, s), W)
			np.save(self.savename + '{}pcf_B_r1_{}.npy'\
				.format(self.order, s), B)

		self._save_scan(scan_args)
		if self.func_type == 'diag':
			self._save_plot(scan_args)


	def scan_correlate_nondiag(self, scan_args):
		losep, hisep = scan_args
		step = self.cf1.grid_spacing
		steps = np.arange(losep, hisep+step, step)
		steps_calib = steps - self.calib
		self.W1_dict, self.B1_dict = {}, {}	# TODO: instantiate these in constructor

		if self.printout:
			print('\n\n\nStarting correlation along main diagonal...\n')
		start = time.time()
		# TODO: chunk this up to store just neighborhoods as it fills upwards in s
		# fills up W1_dict and B1_dict for each step
		for s in steps_calib:
			self.cf1.set_kernel_radius(s)
			self.W1_dict[s], self.B1_dict[s] = self._make_W1_B1()
		end = time.time()
		print(f'\nCorrelation on main diagonal took time: {end-start} seconds\n')

		if self.printout:
			print('\n\n\nCalculating correlation everywhere...\n')
		start = time.time()
		# fills up hyperarray sparsely
		hyperarr_idx_levels = self.order - 1
		hyperarr_shape = steps.shape * hyperarr_idx_levels
		self.corrfunc_hyperarr = np.empty(hyperarr_shape)
		self._recurse_dimensions_correlate(steps_calib, [], hyperarr_idx_levels, 0, 0)
		end = time.time()
		print(f'\nCorrelation everywhere took time: {end-start} seconds\n')
		if self.printout:
			print('Shape of complete correlation function:', self.corrfunc_hyperarr.shape)

		np.round(self.corrfunc_hyperarr, decimals=self.cf1.rounding_precision)

		# TODO: flip the hyperarrays across all symmetry lines 
		# to fill up the whole corrfunc for arbitrary order
		if self.order==3:
			full = self.corrfunc_hyperarr + self.corrfunc_hyperarr.T
			for i in range(len(steps_calib)):
				full[i,i] = self.corrfunc_hyperarr[i,i]
			self.corrfunc_hyperarr = full

		np.save(self.savename + '{}pcf_nondiag_1dbins_range_{}_{}.npy'\
			.format(self.order, losep, hisep), steps_calib)
		np.save(self.savename + '{}pcf_nondiag_corrfunc_range_{}_{}.npy'\
			.format(self.order, losep, hisep), self.corrfunc_hyperarr)




	def load_calib(self):
		if self.printout:
			print('\n\n\nLoading calibration data...\n')

		try:
			if self.fileR:
				calibfile = 'calibration/' + 'calib_ref_{}_gs_{}.npy'\
					.format(remove_ext(self.fileR), self.cf1.grid_spacing)
			elif self.file_gridR:
				# reconstructs expected calib filename from gridR filename
				# changes expected file_gridR model: 'gridR_random0_DR12v5_CMASS_South_gs_5.npy'
				# to make: 'calib_ref_random0_DR12v5_CMASS_South_gs_5.npy'
				fname = '_'.join(self.file_gridR.split('_')[1:])
				calibfile = 'calibration/calib_ref_' + fname				

			self.calib = np.load(calibfile)
			self.calibfile = calibfile # success
			if self.printout:
				print('Loaded CUSTOM calibration data from file:\n', self.calibfile)	

		except:
			self.calib = np.load(self.calibfile)
			if self.printout:
				print('Loaded DEFAULT calibration data from file:\n', self.calibfile)









