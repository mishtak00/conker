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
from argparse import ArgumentParser
from src.correlator import Correlator
import numpy as np




def main():

	print('\n\n\nPlease reference original publication arXiv:XXXXXXXXXXXXX '\
		'when using this software for publishing/redistribution.\n\n\n')
	
	parser = ArgumentParser(description=
		'~~~~~~~~~~~~~~~~~ ( X ) ConKer ( X ) ~~~~~~~~~~~~~~~~~')

	# core conker variables
	parser.add_argument('file1', metavar='INPUT_FILE_1', type=str, 
		help='Name of .fits file with the input data.')
	parser.add_argument('-f0', '--file0', type=str, default=None,
		help='Name of .fits file with the input data.'\
				'Self-correlation if ommitted. Cross-correlation if present.')
	parser.add_argument('-n', '--order', type=int, default=2,
		help='Correlation order wanted.')
	parser.add_argument('-p', '--params_file', type=str, default='params.json', 
		help='Sets custom hyperparameters file.')

	# ancillary behaviors
	parser.add_argument('-s', '--save', action='store_true', 
		help='Grids and .fits output will be automatically saved to an \'out\' folder.')
	parser.add_argument('-v', '--verbose', action='store_true', 
		help='The progress of CenterFinder will be printed out to standard output.')
	parser.add_argument('--scan', nargs=2, type=int,
		help='Calculates correlation function from 1st arg (iclusive) '
				'to 2nd arg (exclusive) by step of grid_spacing.')
	
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

	args = parser.parse_args()

	# deletes the .fits extension and
	# allows for other '.'s in the args.file string
	filename1 = '.'.join(args.file1.split('.')[:-1])
	if args.file0:
		filename0 = '.'.join(args.file0.split('.')[:-1])
	else:
		filename0 = filename1
	# sets up output folder
	savename = 'out_conker_{}_{}/'.format(filename1, filename0)
	try:
		os.mkdir(savename)
	except FileExistsError:
		pass

	# sets up correlator object for run
	corr = Correlator(args.order, args.file1, file0=args.file0,
		params_file=args.params_file, printout=args.verbose,
		save=args.save, savename=savename)


	# creates and customizes instance of CenterFinder object 0
	file0 = args.file1 if not args.file0 else args.file0
	cf0 = CenterFinder(file0, args.weighted_input0, 
		args.params_file, args.save, args.verbose, kernel_type='ball')
	corr.set_cf0(cf0)
	corr.make_cf0(args)
	# del cf0

	# creates and customizes instance of CenterFinder object 1
	cf1 = CenterFinder(args.file1, args.weighted_input1, 
		args.params_file, args.save, args.verbose)
	corr.set_cf1(cf1)
	corr.make_cf1(args)
	# del cf1

	# runs requested correlation
	if args.scan:
		corr.scan_correlate(args)

	else:
		pass
	# 	cf1.make_convolved_grids()
	# 	W1 = cf1.get_centers_grid()
	# 	del cf1

	# 	W = W1 * W0
	# 	if args.save:
	# 		np.save(savename + f'W_r1_{args.kernel_radius1}_r0_{args.kernel_radius0}.npy', W)

	# 	correlation = np.sum(W)
	# 	np.save(savename + f'separation_{args.kernel_radius1}.npy', args.kernel_radius1)
	# 	np.save(savename + f'correlation_{args.kernel_radius1}.npy', correlation)



if __name__ == '__main__':
	main()




