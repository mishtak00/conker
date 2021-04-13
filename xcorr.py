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
from scipy.signal import fftconvolve
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
	parser.add_argument('file1', metavar='MATTER_TRACER_FILE_1', type=str, 
		help='Name of .fits file with the input data.')
	parser.add_argument('-f2', '--matter_tracer_file_2', type=str, default=None,
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
	parser.add_argument('-r2', '--kernel_radius2', type=float, help='Sets kernel radius.')
	parser.add_argument('--show_kernel2', action='store_true', help='Shows 1D kernel plot.')
	kernel_types = parser.add_mutually_exclusive_group()
	kernel_types.add_argument('-e2', '--step_kernel2', nargs='*', type=float,
		help='Fits a step function to the kernel at kernel radius.')
	kernel_types.add_argument('-g2', '--gaussian_kernel2', nargs=1, type=float,
		help='Fits a gaussian function to the kernel at kernel radius.')
	kernel_types.add_argument('-a2', '--wavelet_kernel2', nargs=1, type=float, 
		help='Fits a wavelet function to the kernel at kernel radius.')
	kernel_types.add_argument('-u2', '--custom_kernel2', nargs=1, type=str, 
		help='Fits given custom array to kernel radially.')
	kernel_types.add_argument('-b2', '--ball_kernel2', action='store_true',
		help='Makes a filled sphere of radius kernel_radius.')

	# these define behavior of density grid
	parser.add_argument('-t2', '--vote_threshold2', type=float, 
		help='Centers with number of votes smaller than given argument '\
		'will be discarded from .fits output.')
	parser.add_argument('-w2', '--weighted_input2', action='store_true',
		help='CenterFinder will try to read a fourth column from input data '\
		'and interpret said values as weights.')
	con_or_over = parser.add_mutually_exclusive_group()
	con_or_over.add_argument('-c2', '--density_contrast2', nargs='*',
		help='CenterFinder will subtract the background from the galaxy '\
		'density grid before voting. It will set negative weights to 0 if '\
		'anything is entered after -c.')
	con_or_over.add_argument('-o2', '--overdensity2', action='store_true',
		help='CenterFinder will subtract average density from the galaxy '\
		'density grid before voting.')
	
	# ancillary behaviors
	parser.add_argument('-s', '--save', action='store_true', 
		help='Grids and .fits output will be automatically saved to an \'out\' folder.')
	parser.add_argument('-v', '--verbose', action='store_true', 
		help='The progress of CenterFinder will be printed out to standard output.')

	args = parser.parse_args()

	# deletes the .fits extension and
	# allows for other '.'s in the args.file string
	filename1 = '.'.join(args.file1.split('.')[:-1])
	try:
		os.mkdir('out_{}'.format(filename1))
	except FileExistsError:
		pass
	if args.matter_tracer_file_2:
		filename2 = '.'.join(args.matter_tracer_file_2.split('.')[:-1])
	else:
		filename2 = filename1
	# try:
	# 	os.mkdir('out_ball_{}'.format(filename2))
	# except FileExistsError:
	# 	pass

	# creates and customizes instance of 1st CenterFinder object
	cf1 = CenterFinder(args.file1, args.weighted_input1, 
		args.params_file, args.save, args.verbose)
	if args.kernel_radius1 is not None:
		cf1.set_kernel_radius(args.kernel_radius1)
	if args.show_kernel1:
		cf1.set_show_kernel(args.show_kernel1)
	if args.step_kernel1 is not None:
		cf1.set_kernel_type('step', args.step_kernel1)
	elif args.gaussian_kernel1 is not None:
		cf1.set_kernel_type('gaussian', args.gaussian_kernel1)
	elif args.wavelet_kernel1 is not None:
		cf1.set_kernel_type('wavelet', args.wavelet_kernel1)
	elif args.custom_kernel1 is not None:
		cf1.set_kernel_type('custom', args.custom_kernel1)
	if args.vote_threshold1 is not None:
		cf1.set_vote_threshold(args.vote_threshold1)

	# runs the centerfinding algorithm
	if args.density_contrast1 is not None:
		do_dencon1 = True
		if len(args.density_contrast1)==0:
			keep_neg_wts1 = True
		else:
			keep_neg_wts1 = False
	else:
		do_dencon1 = False
		keep_neg_wts1 = False
	dencon_args1 = (do_dencon1, keep_neg_wts1)
	cf1.find_centers(dencon=dencon_args1,
					overden=args.overdensity1,
					cleanup=False)
	W1 = cf1.centers_grid
	cf1.cleanup()

	# creates and customizes instance of 2nd CenterFinder object
	file2 = args.matter_tracer_file_2 if args.matter_tracer_file_2 else args.file1
	cf2 = CenterFinder(file2, args.weighted_input2, 
		args.params_file, args.save, args.verbose)
	if args.kernel_radius2 is not None:
		cf2.set_kernel_radius(args.kernel_radius2)
	if args.show_kernel2:
		cf2.set_show_kernel(args.show_kernel2)
	if args.step_kernel2 is not None:
		cf2.set_kernel_type('step', args.step_kernel2)
	elif args.gaussian_kernel2 is not None:
		cf2.set_kernel_type('gaussian', args.gaussian_kernel2)
	elif args.wavelet_kernel2 is not None:
		cf2.set_kernel_type('wavelet', args.wavelet_kernel2)
	elif args.custom_kernel2 is not None:
		cf2.set_kernel_type('custom', args.custom_kernel2)
	if args.vote_threshold2 is not None:
		cf2.set_vote_threshold(args.vote_threshold2)

	# runs the centerfinding algorithm
	if args.density_contrast2 is not None:
		do_dencon2 = True
		if len(args.density_contrast2)==0:
			keep_neg_wts2 = True
		else:
			keep_neg_wts2 = False
	else:
		do_dencon2 = False
		keep_neg_wts2 = False
	dencon_args2 = (do_dencon2, keep_neg_wts2)
	cf2.find_centers(dencon=dencon_args2,
					overden=args.overdensity2,
					cleanup=False)
	W2 = cf2.centers_grid
	cf2.cleanup()

	W = fftconvolve(W1, W2, mode='same')
	savename = 'out_xcorr_{}_{}/'.format(filename1, filename2)
	try:
		os.mkdir(savename)
	except FileExistsError:
		pass
	np.save(savename + f'W_r1_{args.kernel_radius1}_r2_{args.kernel_radius2}.npy', W)

	# subprocess.run(['python','cfdriver.py',args.file1,'-r1', args.kernel_radius1])


if __name__ == '__main__':
	main()




