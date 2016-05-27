layerpack
========================

layerpack provides command-line tools for packing/unpacking NetCDF
arrays in slices, which in some circumstances offers a good tradeoff
between compression and reduction in precision. 

Description
-----

These tools provides a lossy compression format for large NetCDF
datasets, and is most applicable to datasets where the range varies
dramatically along one (or several) dimension, but there is only
moderate variation within the remaining dimensions; one example from
the geophysical sciences is gridded mixing ratios of atmopsheric
constituents, which may vary by orders of magnitude vertically or with
time, but have a limited range of variation within vertical layers.

Linear packing is applied to slices of the data, and the scale and
offset parameters are stored as new array variables with the resulting
NetCDF4 file. Considerable compression can be achieved using the
NetCDF4 "deflate" and "shuffle" algorithms. The packed data is
something of an archive format, and is less useful in this state. The
packed data files can be unpacked to their original format using the
tools provided.

Warning
-----
Users should be aware of the loss of accuracy that results from the
packing/unpacking. The individual packed fields are stored with
precision equal to 1/64000 (roughly 0.0016%) of the range within the
packed slice. Users are strongly encouraged to check that the
resulting errors are acceptable for their purposes. Another similarly
effective compression technique is the PPC method (or "bit-grooming")
implemented in the NCO bundle (see
http://nco.sourceforge.net/nco.html#Compression ). We also recommend
that you check which dimensions, or combinations of dimensions, to
pack along give the best compression-accuracy tradeoff. For example,
the NCO tool ncpdq applies a single linear scale-offset to individual
variables, which would result in the most effective compression,
however also the most severe reduction in accuracy. At the other
extreme, if the layer-packing is applied to a large n-dimensional
variable (e.g. n=7), and it is packed along n-1 dimensions, then the
reduction in accuracy may be minimal, but the compression gains may be
small (or even negative, offset by the need to store scale/offset
arrays for each packed variable).


Usage
-----

~~~
## Apply the linear scaling to variables U, V and W at each level of
## dimensions z and time - all other variables are copied unchanged.
$ ncpack -v U,V,W -d z,time orig.nc packed.nc

## Return the packed variables to their original state
$ ncunpack packed.nc unpacked.nc

## Calculate differences between the original and unpacked versions
## of this dataset, presenting errors as the RMSE difference
## normalised by the standard deviation calculated for each level of
## z and time.
$ nccheckdiff -d z,time orig.nc unpacked.nc

## Get help about usage of a command, e.g. ncpack
$ ncpack -h
~~~~

Installation of the command-line utilities
--------

Installation right from the source tree (or via pip from PyPI)::

    $ python setup.py install

The command line utilites are now available.
