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

from argparse import ArgumentParser



class Parser(ArgumentParser):

	def __init__(self):
		super().__init__(description=
			'~~~~~~~~~~~~~~~~~ ( X ) ConKer ( X ) ~~~~~~~~~~~~~~~~~')

		# core conker variables
		self.add_argument('file1', metavar='INPUT_FILE_1', type=str, 
			help='Name of .fits file with the input data.')
		self.add_argument('-f0', '--file0', type=str, default=None,
			help='Name of .fits file with the input data.'\
					'Self-correlation if ommitted. Cross-correlation if present.')
		self.add_argument('-n', '--order', type=int, default=2,
			help='Correlation order wanted.')
		self.add_argument('-p', '--params_file', type=str, default='params.json', 
			help='Sets custom hyperparameters file.')

		# ancillary behaviors
		self.add_argument('-s', '--save', action='store_true', 
			help='Grids and .fits output will be automatically saved to an \'out\' folder.')
		self.add_argument('-v', '--verbose', action='store_true', 
			help='The progress of CenterFinder will be printed out to standard output.')
		self.add_argument('--scan', nargs=2, type=int,
			help='Calculates correlation function from 1st arg (iclusive) '
					'to 2nd arg (exclusive) by step of grid_spacing.')
		
		# these define 1st kernel behavior
		self.add_argument('-r1', '--kernel_radius1', type=float, help='Sets kernel radius.')
		self.add_argument('--show_kernel1', action='store_true', help='Shows 1D kernel plot.')
		kernel_types1 = self.add_mutually_exclusive_group()
		kernel_types1.add_argument('-e1', '--step_kernel1', nargs='*', type=float,
			help='Fits a step function to the kernel at kernel radius.')
		kernel_types1.add_argument('-g1', '--gaussian_kernel1', nargs=1, type=float,
			help='Fits a gaussian function to the kernel at kernel radius.')
		kernel_types1.add_argument('-a1', '--wavelet_kernel1', nargs=1, type=float, 
			help='Fits a wavelet function to the kernel at kernel radius.')
		kernel_types1.add_argument('-u1', '--custom_kernel1', nargs=1, type=str, 
			help='Fits given custom array to kernel radially.')
		kernel_types1.add_argument('-b1', '--ball_kernel1', action='store_true',
			help='Makes a filled sphere of radius kernel_radius.')

		# these define behavior of 1st density grid
		self.add_argument('-t1', '--vote_threshold1', type=float, 
			help='Centers with number of votes smaller than given argument '\
			'will be discarded from .fits output.')
		self.add_argument('-w1', '--weighted_input1', action='store_true',
			help='CenterFinder will try to read a fourth column from input data '\
			'and interpret said values as weights.')
		con_or_over1 = self.add_mutually_exclusive_group()
		con_or_over1.add_argument('-c1', '--density_contrast1', nargs='*',
			help='CenterFinder will subtract the background from the galaxy '\
			'density grid before voting. It will set negative weights to 0 if '\
			'anything is entered after -c.')
		con_or_over1.add_argument('-o1', '--overdensity1', action='store_true',
			help='CenterFinder will subtract average density from the galaxy '\
			'density grid before voting.')

		# these define 1st kernel behavior
		self.add_argument('-r0', '--kernel_radius0', type=float, help='Sets kernel radius.')
		self.add_argument('--show_kernel0', action='store_true', help='Shows 1D kernel plot.')
		kernel_types0 = self.add_mutually_exclusive_group()
		kernel_types0.add_argument('-e0', '--step_kernel0', nargs='*', type=float,
			help='Fits a step function to the kernel at kernel radius.')
		kernel_types0.add_argument('-g0', '--gaussian_kernel0', nargs=1, type=float,
			help='Fits a gaussian function to the kernel at kernel radius.')
		kernel_types0.add_argument('-a0', '--wavelet_kernel0', nargs=1, type=float, 
			help='Fits a wavelet function to the kernel at kernel radius.')
		kernel_types0.add_argument('-u0', '--custom_kernel0', nargs=1, type=str, 
			help='Fits given custom array to kernel radially.')
		kernel_types0.add_argument('-b0', '--ball_kernel0', action='store_true',
			help='Makes a filled sphere of radius kernel_radius.')

		# these define behavior of density grid
		self.add_argument('-t0', '--vote_threshold0', type=float, 
			help='Centers with number of votes smaller than given argument '\
			'will be discarded from .fits output.')
		self.add_argument('-w0', '--weighted_input0', action='store_true',
			help='CenterFinder will try to read a fourth column from input data '\
			'and interpret said values as weights.')
		con_or_over0 = self.add_mutually_exclusive_group()
		con_or_over0.add_argument('-c0', '--density_contrast0', nargs='*',
			help='CenterFinder will subtract the background from the galaxy '\
			'density grid before voting. It will set negative weights to 0 if '\
			'anything is entered after -c.')
		con_or_over0.add_argument('-o0', '--overdensity0', action='store_true',
			help='CenterFinder will subtract average density from the galaxy '\
			'density grid before voting.')




				