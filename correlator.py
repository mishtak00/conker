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

import os
import subprocess
from argparse import ArgumentParser
from src.centerfinder import CenterFinder
import numpy as np


# class XCorrelator(CenterFinder):

# 	def __init__(self, galaxy_file: str, wtd: bool, params_file: str, save: bool, printout: bool,
# 		kernel_radius: float = 110., kernel_type: str = 'step', 
# 		kernel_args: list = [], vote_threshold: float = -inf):

# 		CenterFinder.__init__(self, )

# 	def xcorrelate():
# 		self.M1


def main():

	print('\n\n\nPlease reference original publication arXiv:XXXXXXXXXXXXX '\
		'when using this software for publishing/redistribution.\n\n\n')
	
	parser = ArgumentParser(description=
		'~~~~~~~~~~~~~~~~~ ( X ) ConKer ( X ) ~~~~~~~~~~~~~~~~~')

	# these read input and parameter files
	parser.add_argument('file1', metavar='INPUT_FILE_1', type=str, 
		help='Name of .fits file with the input data.')
	parser.add_argument('-f0', '--input_file_0', type=str, default=None,
		help='Name of .fits file with the input data.')
	parser.add_argument('-p', '--params_file', type=str, default='params.json', 
		help='Sets custom hyperparameters file.')
	
	# these define 1st kernel behavior
	parser.add_argument('-r1', '--kernel_radius1', type=float, help='Sets kernel radius.')
	parser.add_argument('--show_kernel1', action='store_true', help='Shows 1D kernel plot.')
	kernel_types = parser.add_mutually_exclusive_group()
	kernel_types.add_argument('-e1', '--step_kernel1', nargs='*', type=float,
		help='Fits a step function to the kernel at kernel radius.')
	kernel_types.add_argument('-g1', '--gaussian_kernel1', nargs=1, type=float,
		help='Fits a gaussian function to the kernel at kernel radius.')
	kernel_types.add_argument('-a1', '--wavelet_kernel1', nargs=1, type=float, 
		help='Fits a wavelet function to the kernel at kernel radius.')
	kernel_types.add_argument('-u1', '--custom_kernel1', nargs=1, type=str, 
		help='Fits given custom array to kernel radially.')
	kernel_types.add_argument('-b1', '--ball_kernel1', action='store_true',
		help='Makes a filled sphere of radius kernel_radius.')

	# these define behavior of 1st density grid
	parser.add_argument('-t1', '--vote_threshold1', type=float, 
		help='Centers with number of votes smaller than given argument '\
		'will be discarded from .fits output.')
	parser.add_argument('-w1', '--weighted_input1', action='store_true',
		help='CenterFinder will try to read a fourth column from input data '\
		'and interpret said values as weights.')
	con_or_over = parser.add_mutually_exclusive_group()
	con_or_over.add_argument('-c1', '--density_contrast1', nargs='*',
		help='CenterFinder will subtract the background from the galaxy '\
		'density grid before voting. It will set negative weights to 0 if '\
		'anything is entered after -c.')
	con_or_over.add_argument('-o1', '--overdensity1', action='store_true',
		help='CenterFinder will subtract average density from the galaxy '\
		'density grid before voting.')

	# these define 1st kernel behavior
	parser.add_argument('-r0', '--kernel_radius0', type=float, help='Sets kernel radius.')
	parser.add_argument('--show_kernel0', action='store_true', help='Shows 1D kernel plot.')
	kernel_types = parser.add_mutually_exclusive_group()
	kernel_types.add_argument('-e0', '--step_kernel0', nargs='*', type=float,
		help='Fits a step function to the kernel at kernel radius.')
	kernel_types.add_argument('-g0', '--gaussian_kernel0', nargs=1, type=float,
		help='Fits a gaussian function to the kernel at kernel radius.')
	kernel_types.add_argument('-a0', '--wavelet_kernel0', nargs=1, type=float, 
		help='Fits a wavelet function to the kernel at kernel radius.')
	kernel_types.add_argument('-u0', '--custom_kernel0', nargs=1, type=str, 
		help='Fits given custom array to kernel radially.')
	kernel_types.add_argument('-b0', '--ball_kernel0', action='store_true',
		help='Makes a filled sphere of radius kernel_radius.')

	# these define behavior of density grid
	parser.add_argument('-t0', '--vote_threshold0', type=float, 
		help='Centers with number of votes smaller than given argument '\
		'will be discarded from .fits output.')
	parser.add_argument('-w0', '--weighted_input0', action='store_true',
		help='CenterFinder will try to read a fourth column from input data '\
		'and interpret said values as weights.')
	con_or_over = parser.add_mutually_exclusive_group()
	con_or_over.add_argument('-c0', '--density_contrast0', nargs='*',
		help='CenterFinder will subtract the background from the galaxy '\
		'density grid before voting. It will set negative weights to 0 if '\
		'anything is entered after -c.')
	con_or_over.add_argument('-o0', '--overdensity0', action='store_true',
		help='CenterFinder will subtract average density from the galaxy '\
		'density grid before voting.')
	
	# ancillary behaviors
	parser.add_argument('-s', '--save', action='store_true', 
		help='Grids and .fits output will be automatically saved to an \'out\' folder.')
	parser.add_argument('-v', '--verbose', action='store_true', 
		help='The progress of CenterFinder will be printed out to standard output.')
	parser.add_argument('--scan', nargs=2, type=int,
		help='Calculates correlation function from 1st arg (iclusive) '
				'to 2nd arg (exclusive) by step of grid_spacing.')

	args = parser.parse_args()

	# deletes the .fits extension and
	# allows for other '.'s in the args.file string
	filename1 = '.'.join(args.file1.split('.')[:-1])
	if args.input_file_0:
		filename0 = '.'.join(args.input_file_0.split('.')[:-1])
	else:
		filename0 = filename1


	# creates and customizes instance of 0th CenterFinder object
	file0 = args.file1 if not args.input_file_0 else args.input_file_0
	cf0 = CenterFinder(file0, args.weighted_input0, 
		args.params_file, args.save, args.verbose, kernel_type='ball')
	if args.kernel_radius0:
		cf0.set_kernel_radius(args.kernel_radius0)
	if args.show_kernel0:
		cf0.set_show_kernel(args.show_kernel0)
	if args.step_kernel0:
		cf0.set_kernel_type('step', args.step_kernel0)
	elif args.gaussian_kernel0:
		cf0.set_kernel_type('gaussian', args.gaussian_kernel0)
	elif args.wavelet_kernel0:
		cf0.set_kernel_type('wavelet', args.wavelet_kernel0)
	elif args.custom_kernel0:
		cf0.set_kernel_type('custom', args.custom_kernel0)
	if args.vote_threshold0:
		cf0.set_vote_threshold(args.vote_threshold0)

	# runs the centerfinding algorithm
	if args.density_contrast0 is not None:
		do_dencon0 = True
		if len(args.density_contrast0)==0:
			keep_neg_wts0 = True
		else:
			keep_neg_wts0 = False
	else:
		do_dencon0 = False
		keep_neg_wts0 = False
	dencon_args0 = (do_dencon0, keep_neg_wts0)

	cf0.make_grids(dencon=dencon_args0, overden=args.overdensity0)
	# if self-correlation, keeps density and randoms so it calculates only once
	if not args.input_file_0:
		density_grid = cf0.get_density_grid()
		randoms_grid = cf0.get_randoms_grid()
	cf0.make_convolved_grids()
	W0 = cf0.get_centers_grid()
	B0 = cf0.get_background_grid()
	del cf0


	# creates and customizes instance of 1st CenterFinder object
	cf1 = CenterFinder(args.file1, args.weighted_input1, 
		args.params_file, args.save, args.verbose)
	if args.kernel_radius1:
		cf1.set_kernel_radius(args.kernel_radius1)
	if args.show_kernel1:
		cf1.set_show_kernel(args.show_kernel1)
	if args.step_kernel1:
		cf1.set_kernel_type('step', args.step_kernel1)
	elif args.gaussian_kernel1:
		cf1.set_kernel_type('gaussian', args.gaussian_kernel1)
	elif args.wavelet_kernel1:
		cf1.set_kernel_type('wavelet', args.wavelet_kernel1)
	elif args.custom_kernel1:
		cf1.set_kernel_type('custom', args.custom_kernel1)
	if args.vote_threshold1:
		cf1.set_vote_threshold(args.vote_threshold1)

	# runs the centerfinding algorithm
	if args.density_contrast1:
		do_dencon1 = True
		if len(args.density_contrast1)==0:
			keep_neg_wts1 = True
		else:
			keep_neg_wts1 = False
	else:
		do_dencon1 = False
		keep_neg_wts1 = False
	dencon_args1 = (do_dencon1, keep_neg_wts1)
	# doesn't recalculate randoms and density grid if self-cor
	if not args.input_file_0:
		cf1.set_density_grid(density_grid)
		cf1.set_randoms_grid(randoms_grid)
	else:
		cf1.make_grids(dencon=dencon_args1, overden=args.overdensity1)


	savename = 'out_conker_{}_{}/'.format(filename1, filename0)
	try:
		os.mkdir(savename)
	except FileExistsError:
		pass

	if args.scan:
		start, end = args.scan
		corrfunc = {r:None for r in range(start,end,5)}

		for r in corrfunc.keys():
			cf1.set_kernel_radius(r)
			cf1.make_convolved_grids()
			W1 = cf1.get_centers_grid()
			B1 = cf1.get_background_grid()
			cf1.cleanup()

			W = W0 * W1
			B = B0 * B1
			if args.save:
				np.save(savename + f'W_r1_{r}_r0_{args.kernel_radius0}.npy', W)
				np.save(savename + f'B_r1_{r}_r0_{args.kernel_radius0}.npy', B)

			corrfunc[r] = np.sum(W) / np.sum(B)

		del cf1

		separation, correlation = list(corrfunc.keys()), list(corrfunc.values())
		print(separation)
		print(correlation)
		np.save(savename + f'separation_range_{separation[0]}_{separation[-1]}.npy', separation)
		np.save(savename + f'correlation_range_{separation[0]}_{separation[-1]}.npy', correlation)
		import matplotlib.pyplot as plt
		plt.plot(separation, correlation)
		plt.savefig(savename + f'conker_scan_{separation[0]}_{separation[-1]}.png')

	# else:
	# 	cf1.make_convolves_grids()
	# 	W1 = cf1.get_centers_grid()
	# 	del cf1

	# 	W = W1 * W0
	# 	if args.save:
	# 		np.save(savename + f'W_r1_{args.kernel_radius1}_r0_{args.kernel_radius0}.npy', W)

	# 	correlation = np.sum(W)
	# 	np.save(savename + f'separation_{args.kernel_radius1}.npy', args.kernel_radius1)
	# 	np.save(savename + f'correlation_{args.kernel_radius1}.npy', correlation)


	# subprocess.run(['python','cfdriver.py',args.file1,'-r1', args.kernel_radius1])


if __name__ == '__main__':
	main()




