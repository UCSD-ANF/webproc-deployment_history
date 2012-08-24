Introduction
------------

These are scripts for generating maps and movies of the USArray deployment of seismic and inframet sensors over the lifetime of the network.

Dependencies
------------

GMT - the generic mapping tools

ffmpeg - for movies

BRTT Antelope with Python bindings. 5.1-64 or earlier for the time being

### Map Data ###

The script usarray_deploy_map.py requires a number of GMT data files that are not included in this Git repository:

 * alaska.grad
 * alaska.grd
 * deathvalley.grad
 * deathvalley.grd
 * deathvalley.xy
 * saltonsea.grad
 * saltonsea.grd
 * saltonsea.xy
 * usa.grad
 * usa.grd
 * land_ocean.cpt
 * land_only.cpt

These files are available somewhere in some format or other
