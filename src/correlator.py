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

from centerfinder import CenterFinder

class Correlator:
	
	def __init__(
		self, 
		order: int,
		file1: str, 
		wtd1: bool = False,
		r1_list: list = None,
		file0: str = None,
		r0: float = 5.,
		wtd0: bool = False,
		params_file: str = 'params.json',
		save: bool = False,
		printout: bool = False
		):

			self.file1 = file1
			# self-correlation on default
			self.type = 'self'
			# cross-correlation if 2nd file provided
			if file0:
				self.type = 'cross'
				self.file0 = file0
			
			# correlation order has to match radii given
			assert order == len(r1_list) + 1
			self.order = order
			self.r0 = r0
			self.r1_list = r1_list

			self.params_file = params_file
			self.save = save
			self.printout = printout

			self.cf0: CenterFinder = None
			self.cf1: CenterFinder = None


	


	def self_correlate():
		pass

	def cross_correlate():
		pass


