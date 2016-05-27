import numpy
import netCDF4
import re
import itertools
import os
import argparse
import Ncpackutils

"""layerpack.ncunpack: command line utility to unpack layer-packed NetCDF files."""

def main():
    parser = argparse.ArgumentParser(description = "ncunpack: command line utility to unpack layer-packed NetCDF files" )

    ## optional arguments
    parser.add_argument("-V", "--verbose", help="Print extra information",action="store_true")
    parser.add_argument("-O", "--overwrite", help="If this flag is raised, then over-writing existing output file is allowed",action="store_true")
    ## positional arguments
    parser.add_argument('input', default='noinput',help = "Input file")
    ## optional positional arguments
    parser.add_argument('output', default='nooutput', nargs='?',help = "Output file - if not provided, then '.ncp' subscript will be replaced with '.nc', or simply append '.nc'")

    ## process arguments
    args = parser.parse_args()
    verbose = args.verbose
    overwrite = args.overwrite
    infile = args.input
    outfile = args.output
    if outfile == 'nooutput':
        if re.sub(r"\.[nN][cC][pP]$",infile):
            ## replace .nc or .NC with .ncp
            outfile = re.sub(r"\.[nN][cC][pP]$",".nc",infile)
        else:
            ## or just append .ncp if filename doesn't end with .nc
            outfile = outfile + '.nc'
    
    if not os.path.isfile(infile):
        parser.error('Input file does not exist...')

    ## if outfile exists, delete it (only if over-write flag is raised)
    if os.path.isfile(outfile):
        if overwrite:
            os.remove(outfile)
        else:
            raise Exception("Over-write flag not raised, but {} exists".format(outfile))

    Cmp = lambda x, y: int(splitDims.index(x) > splitDims.index(y))*2-1

    ## open the files
    ncin = netCDF4.Dataset(infile, 'r')
    ncout = netCDF4.Dataset(outfile, 'w', format='NETCDF4')

    ## patterns for detecting split variables
    short_pattern = re.compile(".*___short$")    
    split_pattern = re.compile(".*___split$")    
    scale_pattern = re.compile(".*___scale$")    
    offset_pattern = re.compile(".*___offset$")

    ## values to check for and replace
    packed_nan_value     = numpy.uint16(65531)
    packed_neginf_value  = numpy.uint16(65532)
    packed_posinf_value  = numpy.uint16(65533)
    packed_fill_value    = numpy.uint16(65534)
    packed_missing_value = numpy.uint16(65535)

    ## create dimensions in outfile
    splitDims = []
    for k in ncin.dimensions.keys():
        if split_pattern.match(k):
            splitDims.append(re.sub(r"___split$","",k))
        else:
            result = ncout.createDimension(k, len(ncin.dimensions[k]))

    for var in ncin.variables.keys():
        ## skip the scale and offset variables
        if scale_pattern.match(var) or offset_pattern.match(var):
            continue
        ## check if a variable is packed or not
        ispacked = False ## assume that it isn't
        Type = ncin.variables[var].dtype.str
        Type = re.sub(r"[<>]", "", Type)
        if short_pattern.match(var):
            varstem = re.sub(r"___short$","",var)
            varscale = varstem + "___scale"
            varoffset = varstem + "___offset"
            ## the variable is deeemed packed if:
            ## a) It has the appopriate ending ("___short")
            ## b) the corresponding scale and offset variables exist
            ## c) the variable is of the correct type (unsigned short int)
            ## ... additional checks could be applied (e.g. dimensions, type of
            ## scale/offset arrays) ...
            if (varscale in ncin.variables.keys()) and (varoffset in ncin.variables.keys()) and (Type == 'u2'):
                ispacked = True

        if ispacked:
            ## list the dimensions along which to split
            theseSplitDims = Ncpackutils.list_intersection(ncin.variables[var].dimensions,splitDims)
            theseSplitDims.sort(Cmp)
            if verbose:
                print "unpack",varstem,"along",theseSplitDims
            ## get the lengths of these dimensions
            dimLens = [len(ncin.dimensions[d]) for d in theseSplitDims]
            gridRange = [range(l) for l in dimLens]
            splitDimIdxs = list(itertools.product(*gridRange))
            nsplits = len(splitDimIdxs)
            ## get the number of dimensions
            ndim = len(ncin.variables[var].dimensions)
            ## get the indices of which dimensions for this variable are split
            iSplitDim = numpy.array([ list(ncin.variables[var].dimensions).index(d) for d in theseSplitDims ])
            iSplitDim.sort()
            ##
            ## chunksizes = numpy.array([len(ncin.dimensions[d]) for d in ncin.variables[var].dimensions])
            ## chunksizes[iSplitDim] = 1
            ##
            packed = ncin.variables[var][:]
            add_offset = ncin.variables[varoffset][:]
            scale_factor = ncin.variables[varscale][:]
            Data = numpy.ones(packed.shape,numpy.float32)
            ##
            indices = numpy.array([Ellipsis] * ndim)
            for isplit in range(nsplits):
                indices[iSplitDim] = splitDimIdxs[isplit]
                Data[tuple(indices)] = add_offset[tuple(splitDimIdxs[isplit])] + packed[tuple(indices)] * scale_factor[tuple(splitDimIdxs[isplit])]
            ##

            Data[numpy.where(Data == packed_nan_value    )] = numpy.nan
            Data[numpy.where(Data == packed_neginf_value )] = -numpy.inf
            Data[numpy.where(Data == packed_posinf_value )] = numpy.inf

            if '_FillValue' in ncin.variables[var].ncattrs():
                if (Data == packed_fill_value).any():
                    Data[numpy.where(Data == packed_fill_value)] = ncin.variables[var].getncattr('_FillValue')

            if 'missing_value' in ncin.variables[var].ncattrs():
                if (Data == packed_missing_value).any():
                    Data[numpy.where(Data == packed_missing_value)] = ncin.variables[var].missing_value

            ## create an output variable
            varout = ncout.createVariable(varname = varstem, datatype = 'f4', dimensions = ncin.variables[var].dimensions, zlib = True, shuffle = True)
            ## write the data
            varout[:] = Data

            ## copy the attributes
            for att in ncin.variables[var].ncattrs():
                varout.setncattr(att,ncin.variables[var].getncattr(att))

        else:
            if verbose:
                print "copy",var," without packing"

            ## 
            ## create variable
            result = ncout.createVariable(varname = var, datatype = Type, dimensions = ncin.variables[var].dimensions,zlib = True)
            ## copy the attributes
            for att in ncin.variables[var].ncattrs():
                ncout.variables[var].setncattr(att,ncin.variables[var].getncattr(att))
            ## copy the data
            ncout.variables[var][:] = ncin.variables[var][:]

    ##
    ## copy global attributes
    for att in ncin.ncattrs():
        ncout.setncattr(att,ncin.getncattr(att))

    Ncpackutils.set_pack_history_attribute(ncin, ncout)

    ##
    ncin.close()
    ncout.close()
