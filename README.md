# ConKer

*** DEV ***

Get all available cmd line arguments with explanations:
```
python driver.py --help
```

### Rules on input and randoms catalogs
There are 4 arguments for specifying the inputs. All inputs have to be in the 'data' folder.

- Positional (obligatory) argument for input catalog 1
```
python driver.py filename1.fits
```
- Optional argument for input catalog 0. Defining this argument requests a cross-correlation between input catalogs 1 and 0. Omitting this will cause ConKer to calculate an auto-correlation based solely on catalog 1.
```
-f0 filename0.fits OR --file0 filename0.fits
```
- Optional argument for randoms catalog. Defining this argument requests that the randoms grid (pre-convolution background) be calculated from histogramming and factorizing this randoms catalog. Omitting this will cause ConKer to calculate the background grid from input catalog 1.
- It's a good idea to save the randoms grid for a certain catalog at a specific grid spacing so that it's only calculated once. Saved in the output folder and has to be manually put in the 'data' folder if it's to be used as input for future correlations.
```
--randoms_file filenameR.fits --save_randoms
```
Optional argument for specifying a pre-made randoms grid. Putting this argument requests that the randoms grid be loaded from the 'data' folder and used as is. This randoms grid has to have been made from the same catalog as the input catalog and its grid spacing must be the same as the one in the parameters file for current run. Undefined behavior if otherwise.
```
--randoms_grid filename_gridR_gs_5.npy
```

### Outputs
The requested correlation and its related separation distance will be saved in .npy format in a folder named 'out_f1_FILENAME1_...'. If it's a scan of an isotropic correlation, a plot will be saved in the 'plots' folder.

### Rules on boundaries
- All boundaries will be decided from input catalog 0, i.e. they will override the boundaries from randoms catalog and input catalog 1.
- In the case of a separate randoms catalog, only those randoms within file 0 boundaries will be considered in the grid making processes.
- In the case of a cross-correlation, i.e. 2 different input catalogs, the 2 input catalogs will be made to have the same boundaries by overriding file 1 boundaries with those from file 0.


### Various runs

Run a second-order auto-correlation on filename1 using randoms from filenameR, defaults at r1 = 110 Mpch-1 and r0 = grid_spacing (which defaults at 5 Mpch-1):
```
python driver.py filename1.fits --randoms_file filenameR.fits
```

Run 2pcf scan at K0 radius = 10 Mpch-1 and K1 radius from 70 (incusive) to 140 (exclusive) Mpch-1 in steps of grid_spacing, with progress logs on std out (-v), and saves all intermediate grids to output folders:
```
python driver.py filename1.fits --randoms_file filenameR.fits -r0 10 --scan 70 140 -v -s
```

Run 3pcf scan at K0 radius = 10 Mpch-1 and K1 radius from 70 (incusive) to 140 (exclusive) Mpch-1 in steps of grid_spacing, saves just the resulting correlation and separation arrays in the output folder:
```
python driver.py filename1.fits --randoms_file filenameR.fits -v -r0 10 --scan 70 140 -n 3
```

Run single 2pcf calculation at K1 radius = 110 and K0 radius = 10:
```
python driver.py mock_cmassDR9_north_3001.fits -v -r0 10 -r1 110
```

