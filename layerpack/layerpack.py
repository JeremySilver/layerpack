"""main help file"""

def __main__():

    print """layerpack: command-line tools for packing/unpacking NetCDF arrays in slices

Usage:
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

Options:
  -h --help                         Help
  -V --verbose                      Verbose output [only: ncpack, ncunpack]
  -v --variables v1,v2,...,vn       Apply to variables v1,v2,...,vn [only: ncpack]
  -d --dimensions d1,d2,...,dm      Pack for each slice across dimensions d1,d2,...,dm [only: ncpack]
  -L --deflatelevel=l               Compress data using deflate l=1,...,9 [only: ncpack, ncunpack]
  -O --overwrite                    Overwrite existing files [only: ncpack, ncunpack]
"""

