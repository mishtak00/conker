# ConKer

*** DEV ***

Get all available terminal arguments with explanations:
```
python driver.py --help
```

Run 2pcf scan at K0 radius = 10 Mpch-1 and K1 radius from 75 (incusive) to 140 (exclusive) Mpch-1 in steps of grid_spacing (default 5 Mpch-1):
```
python driver.py mock_cmassDR9_north_3001.fits -v -r0 10 --scan 70 140
```

Run 3pcf scan at K0 radius = 10 Mpch-1 and K1 radius from 75 (incusive) to 140 (exclusive) Mpch-1 in steps of grid_spacing (default 5 Mpch-1):
```
python driver.py mock_cmassDR9_north_3001.fits -v -r0 10 --scan 70 140 -n 3
```

Run single 2pcf calculation at K1 radius = 110 and K0 radius = 10:
```
python driver.py mock_cmassDR9_north_3001.fits -v -r0 10 -r1 110
```

