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
import matplotlib.pyplot as plt
from progress.bar import ChargingBar
from math import sqrt
from scipy.signal import fftconvolve
from src.parser import CalibrationParser
from src.kernel import Kernel
from src.utils import *
import time



def main():

	parser = CalibrationParser()
	args = parser.parse_args()
	if args.scan:
		try:
			assert args.scan[0] < args.scan[1]
		except AssertionError:
			print('AssertionError: scan args have to be in ascending order')
			sys.exit(1)

	print('Loading input data...')

	indir = 'data/'
	outdir = 'calibration/'
	input_file = args.fileR
	params_file = 'params.json'
	filename = remove_ext(input_file)

	# for the main_calib_loop restart mechanic
	kernel_center_cell_pts_floor = 5
	restart_limit = 30

	loadname = indir + input_file
	G_ra, G_dec, G_redshift = load_data(loadname)
	# weights don't matter for calibration
	G_weights = np.ones(len(G_ra), dtype=float)
	# gets cosmology and other hyperparameters
	cosmology, grid_spacing, rounding_precision \
		= load_hyperparameters(params_file)
	# calculates lookup tables for fast conversion from r to z and vice versa
	LUT_radii, LUT_redshifts = interpolate_r_z(G_redshift.min(), 
		G_redshift.max(), cosmology)
	G_radii = LUT_radii(G_redshift)

	print('Input loaded successfully...')
	print('Initializing calibration...')

	xyzs = np.array(sky2cartesian(G_ra, G_dec, G_redshift, LUT_radii))
	galaxies_cartesian = xyzs.T  # each point is represented by (x, y, z)

	bin_counts_3d = np.array(
		[np.ceil((xyzs[i].max() - xyzs[i].min()) / grid_spacing) 
			for i in range(len(xyzs))], 
		dtype=int)
	histo_grid, histo_grid_edges = np.histogramdd(
		galaxies_cartesian, 
		bins=bin_counts_3d, 
		weights=G_weights
		)
	print('Original grid idx ranges: (0,{}) (0,{}) (0,{})'\
		.format(*histo_grid.shape))

	xs, ys, zs = xyzs
	x_edges, y_edges, z_edges = histo_grid_edges

	rstart, rend, step = 110, 115, 5
	if args.scan:
		rstart, rend = args.scan
		# have to invert scan boundaries
		# because of the optimization for continuous querying
		# see IMPORTANT below
		step = -grid_spacing
		rstart_temp = rend + step
		rend = rstart + step
		rstart = rstart_temp

	# dict will save tuple per assumed distance: (actual,min,max)
	rvsr = {r:(None,None,None) \
		for r in np.arange(rstart, rend, step)}

	def main_calib_loop(move_idxs_by) -> bool:
		nonlocal outdir, filename
		nonlocal grid_spacing, histo_grid
		nonlocal xs, ys, zs
		nonlocal x_edges, y_edges, z_edges
		nonlocal rstart, rend, step
		nonlocal rvsr

		for r in np.arange(rstart, rend, step):
			print('Calibrating for r =', r)
			actual_ds = []

			kernel = Kernel('step', r, grid_spacing, False, False)
			kernel_grid = kernel.get_grid()
			kc_idx = kernel.get_kernel_center()[0] # kernel center idx
			print('Got kernel grid with shape:', kernel_grid.shape)

			ks = kernel_grid.shape
			f = ks[0] // 2
			xmid, ymid, zmid = np.array(histo_grid.shape) // 2 + move_idxs_by
			xmini, xmaxi = xmid-f, xmid+f+1
			ymini, ymaxi = ymid-f, ymid+f+1
			zmini, zmaxi = zmid-f, zmid+f+1
			box = histo_grid[xmini:xmaxi, ymini:ymaxi, zmini:zmaxi]
			print('Selected subgrid idx ranges: ({},{}) ({},{}) ({},{})'\
				.format(xmini,xmaxi,ymini,ymaxi,zmini,zmaxi))
			print('Got central subgrid with shape', box.shape)
			print('Number of points in subgrid:', np.sum(box))

			# locates correct cell idx boundaries in histo grid
			xcmin_all, xcmax_all = x_edges[xmini], x_edges[xmaxi]
			ycmin_all, ycmax_all = y_edges[ymini], y_edges[ymaxi]
			zcmin_all, zcmax_all = z_edges[zmini], z_edges[zmaxi]

			# finds all points within this coordinate range
			pidxs_all = np.asarray(
				(xs >= xcmin_all) & (xs < xcmax_all) \
				& (ys >= ycmin_all) & (ys < ycmax_all) \
				& (zs >= zcmin_all) & (zs < zcmax_all)
				).nonzero()

			# IMPORTANT: only for use in conjunction with descending kernel radii
			# reassigns data coord arrays to just pts within curr bounds 
			# much faster querying in future loop iters
			xs, ys, zs = xs[pidxs_all], ys[pidxs_all], zs[pidxs_all]

			# calculates coordinates of kernel center
			kcxmin, kcxmax = x_edges[xmini+kc_idx], x_edges[xmini+kc_idx+1]
			kcymin, kcymax = y_edges[ymini+kc_idx], y_edges[ymini+kc_idx+1]
			kczmin, kczmax = z_edges[zmini+kc_idx], z_edges[zmini+kc_idx+1]

			# finds all points within the central cell of the kernel
			pidxs_kc = np.asarray(
				(xs >= kcxmin) & (xs < kcxmax) \
				& (ys >= kcymin) & (ys < kcymax) \
				& (zs >= kczmin) & (zs < kczmax)
				).nonzero()
			pts_in_kernel_center_cell = list(zip(xs[pidxs_kc], ys[pidxs_kc], zs[pidxs_kc]))
			
			# need to restart someplace else in the catalog
			# if there's no points in the center cell 
			# of the kernel in current configuration
			# if len(pts_in_kernel_center_cell) == 0:
			# 	# allows only 10 tries to restart
			# 	if (move_idxs_by < restart_limit):
			# 		print('Couldn\'t find enough points in center cell, '\
			# 				'restart nr. {}...'.format(move_idxs_by+1))
			# 		return main_calib_loop(move_idxs_by+1)
			# 	else:
			# 		print('Reached restart limit, exiting...')
			# 		return False
			
			if (len(pts_in_kernel_center_cell) < kernel_center_cell_pts_floor):
				# allows only 10 tries to restart
				if (move_idxs_by < restart_limit):
					print('Couldn\'t find enough points in center cell, '\
							'restart nr. {}...'.format(move_idxs_by+1))
					return main_calib_loop(move_idxs_by+1)
				else:
					print('Reached restart limit, exiting...')
					return False
			else:
				print(f"Found '{len(pts_in_kernel_center_cell)}' points in kernel center cell")


			# traverses inscribed surface in kernel grid
			idxs = np.asarray(kernel_grid!=0).nonzero()
			enough = 10**5 # upper bound for nr of pts for calib

			# # pretty prints progress bar
			progbar = ChargingBar('Calculating actual distances...', 
				max=len(idxs[0]),
				suffix = '%(percent).1f%% - %(elapsed)ds')

			cut = False
			start = time.time()
			for i, j, k in zip(*idxs):

				# skips if respective cell in histo grid empty
				if box[i,j,k] == 0:
					progbar.next()
					continue

				# locates correct cell idx boundaries in histo grid
				xi, yi, zi = xmini+i, ymini+j, zmini+k
				xcmin, xcmax = x_edges[xi], x_edges[xi+1]
				ycmin, ycmax = y_edges[yi], y_edges[yi+1]
				zcmin, zcmax = z_edges[zi], z_edges[zi+1]

				# finds all points within this i,j,k cell
				pidxs = np.asarray(
					(xs >= xcmin) & (xs < xcmax) \
					& (ys >= ycmin) & (ys < ycmax) \
					& (zs >= zcmin) & (zs < zcmax)
					).nonzero()
				pts_in_cell = list(zip(xs[pidxs], ys[pidxs], zs[pidxs]))

				# registers actual distance vs assumed distance
				assumed_d = r
				for x, y, z in pts_in_cell:
					for kcx, kcy, kcz in pts_in_kernel_center_cell:
						actual_d = sqrt((x-kcx)**2 + (y-kcy)**2 + (z-kcz)**2)
						actual_ds.append(actual_d)

				# can end prematurely if enough pts have been considered
				if len(actual_ds) >= enough:
					cut = True
					break

				progbar.next()

			end = time.time()
			progbar.finish()

			if cut:
				print('Enough distances were calculated for calibration...')
			print('Final nr of pts considered for calibration (were on shell):',
				 len(actual_ds))
			print('Distances calculated in {} seconds'.format(end-start))

			print('Calculating mean for calibration...')
			actual_ds = np.array(actual_ds)
			mean = np.mean(actual_ds)
			meanerrs = actual_ds - mean
			lomeanerrs = meanerrs[meanerrs < 0]
			himeanerrs = meanerrs[meanerrs >= 0]
			loerr = sqrt(np.mean(lomeanerrs**2))
			hierr = sqrt(np.mean(himeanerrs**2))
			# rvsr[r] = (mean,min(actual_ds),max(actual_ds))
			rvsr[r] = (mean, loerr, hierr)

			print('Mean actual distance:', mean)
			mean = np.array(mean)
			np.save('{}calib_mean_{}_gs_{}_r1_{}.npy'\
				.format(outdir, filename, grid_spacing, r), 
				mean)

			print('Histogramming and saving...')
			plt.hist(actual_ds, 
				bins=20, 
				edgecolor='black', 
				label=r'$r_1=$'+str(r)+r'$h^{-1}Mpc$')
			plt.title('ConKer: Actual distance vs. discrete distance')
			plt.xlabel(r'$r$')
			plt.ylabel('Count')
			plt.legend()
			plt.savefig('{}calib_hist_{}_gs_{}_r1_{}.png'\
				.format(outdir, filename, grid_spacing, r), 
				dpi=300)
			plt.clf()

		return True

	# this inline func is a facility for restarting the calibration
	# in the case where there's no points in the central cell of the kernel
	# given a particular configuration in space, the default calibration
	# procedure starts in the middle of the survey space
	if(not main_calib_loop(0)):
		print('Calibration procedure ended unsuccessfully: '\
				'couldn\'t find points in kernel center cell')
		return

	rs = np.array(list(rvsr.keys()))
	vals = np.array(list(rvsr.values()))
	means = vals[:,0]
	loerrs = vals[:,1]
	uperrs = vals[:,2]
	meanerr = np.round(np.mean(rs - means), 
		decimals=rounding_precision)
	np.save('{}calib_ref_{}_gs_{}.npy'\
			.format(outdir, filename, grid_spacing), 
			meanerr)

	plt.errorbar(rs, means, yerr=[loerrs,uperrs], label='uncalibrated', c='orange')
	plt.errorbar(rs, means+meanerr, label='calibrated', c='cornflowerblue')
	plt.title('ConKer: Actual distances vs. discrete distances')
	plt.xlabel(r'$r_{discrete}$ $[h^{-1}Mpc]$')
	plt.ylabel(r'$r_{actual}$ $[h^{-1}Mpc]$')
	plt.xticks(ticks=rs[::2])
	plt.yticks(ticks=rs)
	plt.grid(linestyle=':')
	plt.legend()
	lo, hi = min(rvsr.keys()), max(rvsr.keys())
	plt.savefig('{}calib_vsplot_{}_gs_{}_r1_{}_{}.png'\
		.format(outdir, filename, grid_spacing, lo, hi), 
		dpi=300)
	plt.close()

	print('Mean separation error to {} decimal places:'\
		.format(rounding_precision), meanerr)
	print('Calibration ended successfully...')




if __name__ == '__main__':
	main()






