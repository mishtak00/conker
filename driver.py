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
from src.parser import Parser
from src.utils import remove_ext
from src.centerfinder import CenterFinder
from src.correlator import Correlator



def main():

	print('\n\n\nPlease reference original publication arXiv:XXXXXXXXXXXXX '\
		'when using this software for publishing/redistribution.\n\n\n')

	parser = Parser()
	args = parser.parse_args()

	# preps output folder
	filename1 = remove_ext(args.file1)
	savename = 'out_f1_{}'.format(filename1)

	if args.file0:
		filename0 = remove_ext(args.file0)
		savename += '_f0_{}'.format(filename0)
	else:
		filename0 = filename1
	if args.randoms_file:
		filenameR = remove_ext(args.randoms_file)
		savename += '_fR_{}'.format(filenameR)
	elif args.randoms_grid:
		filenameR = remove_ext(args.randoms_grid)
		savename += '_gR_{}'.format(filenameR)
	else:
		filenameR = filename1

	savename += '/'
	try:
		os.mkdir(savename)
	except FileExistsError:
		pass

	# sets up correlator object for run
	corr = Correlator(args.order, args.file1, file0=args.file0,
		fileR=args.randoms_file, file_gridR=args.randoms_grid, 
		nondiag=args.nondiag, save_randoms=args.save_randoms,
		params_file=args.params_file, printout=args.verbose,
		save=args.save, savename=savename)

	# creates and puts cf object instance for randoms (non-conv background)
	fileR = args.file1 if not args.randoms_file else args.randoms_file
	cfR = CenterFinder(fileR, args.wtd_randoms,
		args.params_file, args.save, args.verbose,
		dont_factorize=args.dont_factorize_randoms)
	corr.set_cfR(cfR)
	# creates and customizes instance of CenterFinder object 0
	file0 = args.file1 if not args.file0 else args.file0
	cf0 = CenterFinder(file0, args.wtd_input0, 
		args.params_file, args.save, args.verbose, kernel_type='ball')
	corr.set_cf0(cf0)
	# creates and customizes instance of CenterFinder object 1
	cf1 = CenterFinder(args.file1, args.wtd_input1, 
		args.params_file, args.save, args.verbose)
	corr.set_cf1(cf1)

	# histograms input catalog 0 to get boundaries
	corr.prep_cf0(args)

	# makes the randoms grid with input data boundaries
	# sets Correlator object's randoms_grid attribute
	# normalization: ((P_r * P_ad) / NR) * (ND / NR)
	corr.make_cfR(corr.get_cf0().get_density_grid_edges())

	# makes cf0 with randoms grid from cfR
	corr.make_cf0(args)

	# makes cf1 with randoms grid from cfR
	corr.make_cf1(args)

	# finds custom calib file or falls back to default
	corr.load_calib()

	# runs requested correlation and saves output
	# scan command overrides single command
	if args.scan:
		if args.order == 2 or (args.order>2 and not args.noniso):
			corr.scan_correlate_diag(args.scan)
		elif args.order > 2 and args.noniso:
			corr.scan_correlate_nondiag(args.scan)
	else:
		# TODO: integrate higher orders with args from -r1
		corr.single_correlate()



if __name__ == '__main__':
	main()




