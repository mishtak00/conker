# ConKer

*** DEV ***

Get all available cmd line arguments with explanations:
```
python driver.py --help
```

### Rules on input and randoms catalogs
There are 3 arguments for specifying the input catalogs:
- Positional (obligatory) argument for input catalog 1
```
python driver.py filename1.fits
```
- Optional argument for input catalog 0. Defining this argument requests a cross-correlation between input catalogs 1 and 0. Omitting this will cause ConKer to calculate an auto-correlation based solely on catalog 1.
```
-f0 filename0.fits OR --file0 filename0.fits
```
- Optional argument for randoms catalog. Defining this argument requests that the randoms grid (pre-convolution background) be calculated from histogramming and factorizing this randoms catalog. Omitting this will cause ConKer to calculate the background grid from input catalog 1.
```
--randoms_file filenameR.fits
```

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

