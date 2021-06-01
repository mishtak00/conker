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

import os, sys
from src.parser import Parser
from src.centerfinder import CenterFinder
from src.correlator import Correlator



def main():

	print('\n\n\nPlease reference original publication arXiv:XXXXXXXXXXXXX '\
		'when using this software for publishing/redistribution.\n\n\n')

	parser = Parser()
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

	# creates and customizes instance of CenterFinder object 1
	cf1 = CenterFinder(args.file1, args.weighted_input1, 
		args.params_file, args.save, args.verbose)
	corr.set_cf1(cf1)
	corr.make_cf1(args)

	# runs requested correlation and saves output
	# scan command overrides single command
	if args.scan:
		corr.scan_correlate(args.scan)

	else:
		corr.single_correlate()



if __name__ == '__main__':
	main()




