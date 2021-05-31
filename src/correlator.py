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

import sys
import numpy as np
from .centerfinder import CenterFinder



class Correlator:
	

	def __init__(
		self, 
		order: int,
		file1: str, 
		# wtd1: bool = False,
		# r1_list: list = None,
		file0: str = None,
		# r0: float = 5.,
		# wtd0: bool = False,
		params_file: str = 'params.json',
		save: bool = False,
		savename: str = None,
		printout: bool = False
		):

		try:
			assert order >= 2
		except AssertionError:
			print('AssertionError: '\
				f'Correlation order has to be >= 2, was given {order}')
			sys.exit(1)

		self.order = order
		self.file1 = file1
		# auto-correlation on default
		self.type = 'auto'
		# cross-correlation if 2nd file provided
		if file0:
			self.type = 'cross'
			self.file0 = file0
		
		# # correlation order has to match radii given
		# assert order == len(r1_list) + 1
		# self.r0 = r0
		# self.r1_list = r1_list

		self.params_file = params_file
		self.save = save
		self.savename = savename
		self.printout = printout

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



	def set_cf0(self, cf0):
		self.cf0 = cf0


	def get_cf0(self):
		return self.cf0


	def set_cf1(self, cf1):
		self.cf1 = cf1


	def get_cf1(self):
		return cf1


	def _customize_cf0(self, args):
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


	def make_cf0(self, args):
		# defaults kernel_radius to 1/2 grid_spacing for 0th centerfinder
		self.cf0.set_kernel_radius(self.cf0.grid_spacing / 2)
		# legacy, customizes the cf object and the bckgrnd subtraction
		self._customize_cf0(args)
		# runs the centerfinding algorithm
		self.cf0.make_grids(dencon=self.dencon_args0, 
			overden=args.overdensity0)
		# if auto-correlation, keeps density and randoms so it calculates only once
		if self.type == 'auto':
			self.density_grid = self.cf0.get_density_grid()
			self.randoms_grid = self.cf0.get_randoms_grid()
		# there's no need for convolving 
		# if the whole kernel is just one cell
		if self.cf0.kernel_radius == self.cf0.grid_spacing/2:
			self.W0 = self.cf0.get_density_grid()
			self.B0 = self.cf0.get_randoms_grid()
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
		# legacy, customizes the cf object and the bckgrnd subtraction
		self._customize_cf1(args)
		# doesn't recalculate randoms and density grid if auto-cor
		if self.type == 'auto':
			self.cf1.set_density_grid(self.density_grid)
			self.cf1.set_randoms_grid(self.randoms_grid)
		# runs the centerfinding algorithm again if x-cor
		else:
			self.cf1.make_grids(dencon=self.dencon_args1, 
				overden=args.overdensity1)


	@staticmethod
	def _reduce_mult_till_order(cf1: CenterFinder, order: int) \
		-> (np.ndarray, np.ndarray):
		"""
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


	def _correlate(self):
		self.cf1.make_convolved_grids()
		self.W1_prod, self.B1_prod = Correlator\
			._reduce_mult_till_order(self.cf1, self.order)
		self.cf1.cleanup()
		W = self.W0 * self.W1_prod
		B = self.B0 * self.B1_prod
		if self.save:
			np.save(self.savename + '{}pcf_W_r1_{}_r0_{}.npy'\
				.format(self.order, r, args.kernel_radius0), W)
			np.save(self.savename + '{}pcf_B_r1_{}_r0_{}.npy'\
				.format(self.order, r, args.kernel_radius0), B)
		W = np.sum(W)
		B = np.sum(B)

		return W, B


	def save_single(self):
		separation = list(self.corrfunc.keys())
		correlation = list(self.corrfunc.values())
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
		W, B = self._correlate()
		self.corrfunc[s] = W / B
		self.save_single()


	def save_scan(self):
		separation = list(self.corrfunc.keys())
		correlation = list(self.corrfunc.values())
		if self.printout:
			print('Separation array:\n', separation)
			print('Correlation array:\n', correlation)
		np.save(self.savename + '{}pcf_sep_range_{}_{}.npy'\
			.format(self.order, separation[0], separation[-1]), separation)
		np.save(self.savename + '{}pcf_corr_range_{}_{}.npy'\
			.format(self.order, separation[0], separation[-1]), correlation)
		import matplotlib.pyplot as plt
		plt.plot(separation, correlation)
		plt.savefig(self.savename + '{}pcf_conker_scan_{}_{}.png'\
			.format(self.order, separation[0], separation[-1]))


	def scan_correlate(self, scan_args):
		start, end = scan_args
		self.corrfunc = {s: None for s in range(start, end, self.cf1.grid_spacing)}

		for s in self.corrfunc.keys():
			self.cf1.set_kernel_radius(s)
			W, B = self._correlate()
			self.corrfunc[s] = W / B

		self.save_scan()












