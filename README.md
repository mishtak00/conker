# ConKer

*** DEV ***

Run scan at K2 radius = 10 and K1 radius from 75 (incusive) to 140 (exclusive) in steps of grid_spacing:
```
python correlator.py mock_cmassDR9_north_3001.fits -v -r0 10 -c0 --scan 70 140
```

Run one correlation at K1 radius = 110 and K2 radius = 10:
```
python correlator.py mock_cmassDR9_north_3001.fits -v -r1 110 -c1 -r0 10 -c0
```
