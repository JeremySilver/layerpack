import numpy
import netCDF4
# import h5netcdf.legacyapi as netCDF4
import re
import itertools
import os
import argparse
import Ncpackutils

"""layerpack.ncpack: command line utility to compress NetCDF files to layer-packed NetCDF format."""

def main():
    parser = argparse.ArgumentParser("ncpack: command line utility to compress NetCDF files to layer-packed NetCDF format")

    ## optional arguments
    parser.add_argument("-V", "--verbose", help="Print extra information",action="store_true")
    parser.add_argument("-O", "--overwrite", help="If this flag is raised, then over-writing existing output file is allowed",action="store_true")
    parser.add_argument("-v", "--variables", help="Variables to split, separated by commas, with no spaces in between. To split all variables, then leave out this argument", default = '__all__')
    parser.add_argument("-L", "--deflatelevel", help="Deflate compression level (should be an integer from 1-9)", default = 4)
    ## non-optional arguments
    parser.add_argument("-d", "--dimensions", help="Dimensions to split along, separated by commas, with no spaces in between",required=True)
    ## positional arguments
    parser.add_argument('input', default='noinput',help = "Input file")
    ## optional positional arguments
    parser.add_argument('output', default='nooutput', nargs='?',help = "Output file - if not provided, then '.nc' subscript will be replaced with '.ncp', or '.ncp' will be appended to the input name")

    ## process arguments
    args = parser.parse_args()
    verbose = args.verbose
    overwrite = args.overwrite
    splitVarsOrig = args.variables
    splitVars = splitVarsOrig.split(',')
    splitDimsOrig = args.dimensions
    splitDims = splitDimsOrig.split(',')
    level = int(args.deflatelevel)
    if not ((1 <= level) and (level <= 9)):
        parser.error('Deflate level should be between 1 and 9 (inclusive)...')

    infile = args.input
    outfile = args.output
    ## if no output filename is given
    if outfile == 'nooutput':
        if re.sub(r"\.[nN][cC]$",infile):
            ## replace .nc or .NC with .ncp
            outfile = re.sub(r"\.[nN][cC]$",".ncp",infile)
        else:
            ## or just append .ncp if filename doesn't end with .nc
            outfile = outfile + '.ncp'
    
    if not os.path.isfile(infile):
        parser.error('Input file does not exist...')

    Cmp = lambda x, y: int(splitDims.index(x) > splitDims.index(y))*2-1

    ## if outfile exists, delete it (only if over-write flag is raised)
    if os.path.isfile(outfile):
        if overwrite:
            os.remove(outfile)
        else:
            raise Exception("Over-write flag not raised, but {} exists".format(outfile))

    ## open the files
    ncin = netCDF4.Dataset(infile, 'r')
    ncout = netCDF4.Dataset(outfile, 'w', format='NETCDF4')

    ## check that all the required dimensions appear
    for d in splitDims:
        if not (d in ncin.dimensions.keys()):
            raise Exception("Dimension {} not found in {}".format(d,infile))

    variableNames = ncin.variables.keys()

    ## if no variables are specified, apply to all
    if len(splitVars) == 1 and splitVars[0] == '__all__':
        splitVars = ncin.variables.keys()
    else:
        ## check that all the required variables appear
        for v in splitVars:
            if not (v in ncin.variables.keys()):
                raise Exception("Variable {} not found in {}".format(v,infile))

    ## create dimensions in outfile
    for k in ncin.dimensions.keys():
        result = ncout.createDimension(k, len(ncin.dimensions[k]))
        if k in splitDims:
            splitName = '%s___split' % k
            result = ncout.createDimension(splitName, 1)

    shuffle_types = ['u','i','c','s','b','h','l']
    compress_types = ['f4','f8']
    if (len(splitVars) == 1) and splitVars[0] == '*':
        splitVars = ncin.variables.keys()

    ## write out the coordinate variables first
    Vars = ncin.variables.keys()
    for dim in ncin.dimensions.keys():
        if dim in Vars:
            ## insert the coordinate variable at the start
            Vars = [dim] + [v for i,v in enumerate(Vars) if v != dim]
            
    # coordVars = list(set(Vars).intersection(set(ncin.dimensions.keys())))
    # coordVars.sort()
    # noncoordVars = list(set(Vars) - set(ncin.dimensions.keys()))
    # noncoordVars.sort()
    # allVars = coordVars + noncoordVars

    for var in Vars:
        ## extract the type
        Type = ncin.variables[var].dtype.str
        Type = re.sub(r"[<>]", "", Type)
        ## decide whether to shuffle or not
        doShuffle = Type in shuffle_types
        doShuffle = True
        ##
        ## split only if:
        ## a) requested,
        ## b) if there are some dimensions to split along
        ## c) it is a float or a double
        ## d) this variable has at least one dimension to *NOT* split along
        if (var in splitVars) and (len(Ncpackutils.list_intersection(ncin.variables[var].dimensions,splitDims)) > 0) and (Type in compress_types) and (len(Ncpackutils.list_difference(ncin.variables[var].dimensions,splitDims)) > 0):
            ## copy over variables that ARE split
            ##        
            ## list the dimensions along which to split
            theseSplitDims = Ncpackutils.list_intersection(ncin.variables[var].dimensions,splitDims)
            theseSplitDims.sort(Cmp)
            ## get the lengths of these dimensions
            dimLens = [len(ncin.dimensions[d]) for d in theseSplitDims]
            ## allow for the case of multiple dimensions to split along
            gridRange = tuple([tuple(range(l)) for l in dimLens])
            ## get all combinations of indices along the split dimensions
            splitDimIdxs = list(itertools.product(*gridRange))
            ## get the number of dimensions
            ndim = len(ncin.variables[var].dimensions)

            if verbose:
                print "pack",var,"along",theseSplitDims

            ## get the indices of which dimensions for this variable are split
            iSplitDim = numpy.array([ list(ncin.variables[var].dimensions).index(d) for d in theseSplitDims ])
            iSplitDim.sort()
            gridRange = [range(l) for l in dimLens]
            splitDimIdxs = list(itertools.product(*gridRange))
            indices = numpy.array([Ellipsis] * ndim)
            chunksizes = numpy.array([len(ncin.dimensions[d]) for d in ncin.variables[var].dimensions])
            chunksizes[iSplitDim] = 1

            ## define new variables (one with the data in packed format, and two
            ## new arrays, with the scale and offset factors)
            shortname = "%s___short" % (var)
            scalename = "%s___scale" % (var)
            offsetname = "%s___offset" % (var)
            ncout.createVariable(varname = shortname, datatype = 'u2', dimensions = ncin.variables[var].dimensions, zlib = True, shuffle = True, chunksizes = chunksizes, complevel = level)
            ncout.createVariable(varname = scalename, datatype = 'f4', dimensions = tuple(theseSplitDims), zlib = True, complevel = level, shuffle = False)
            ncout.createVariable(varname = offsetname, datatype = 'f4', dimensions = tuple(theseSplitDims), zlib = True, complevel = level, shuffle = False)

            nsplits = len(splitDimIdxs)
            ##
            if 'missing_value' in ncin.variables[var].ncattrs():
                missing_value = ncin.variables[var].missing_value
            else:
                missing_value = None
            ##
            if '_FillValue' in ncin.variables[var].ncattrs():
                fill_value = ncin.variables[var]._FillValue
            else:
                fill_value = None
            ##
            ncin.variables[var].set_auto_mask(False)
            Data = ncin.variables[var][:]
            packedDimLens = tuple([ len(ncin.dimensions[d]) for d in list(theseSplitDims)])

            packed = numpy.ones(Data.shape,numpy.uint16)
            add_offset = numpy.ones(packedDimLens,numpy.float32)
            scale_factor = numpy.ones(packedDimLens,numpy.float32)

            for isplit in range(nsplits):
                indices[iSplitDim] = splitDimIdxs[isplit]
                ## compress the slice
                add_offset[tuple(splitDimIdxs[isplit])], scale_factor[tuple(splitDimIdxs[isplit])], packed[tuple(indices)] = Ncpackutils.pack_to_short(Data[tuple(indices)], missing_value = missing_value)

            ## write the data
            ncout.variables[shortname][:] = packed
            ncout.variables[offsetname][:] = add_offset
            ncout.variables[scalename][:] = scale_factor

            del Data, packed, add_offset, scale_factor

            ## copy the attributes
            for att in ncin.variables[var].ncattrs():
                if not (att[0] == '_'):
                    ncout.variables[shortname].setncattr(att,ncin.variables[var].getncattr(att))
            ncout.sync()

        else: ## copy over variables that are not split
            if verbose:
                print "copy",var," without packing"

            ## 
            ## create variable
            chunksizes = ncin.variables[var].shape
            result = ncout.createVariable(varname = var, datatype = Type, dimensions = ncin.variables[var].dimensions, zlib = True, shuffle = doShuffle, complevel = level, chunksizes = chunksizes)
            ## copy the attributes
            for att in ncin.variables[var].ncattrs():
                ncout.variables[var].setncattr(att,ncin.variables[var].getncattr(att))
            ## copy the data
            ncout.variables[var][:] = ncin.variables[var][:]
            ncout.sync()
    ##
    ## copy global attributes
    for att in ncin.ncattrs():
        ncout.setncattr(att,ncin.getncattr(att))

    Ncpackutils.set_pack_history_attribute(ncin, ncout)

    ##
    ncin.close()
    ncout.close()


