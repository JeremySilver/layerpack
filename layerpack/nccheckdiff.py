import numpy
import netCDF4
import re
import itertools
import os
import argparse
import Ncpackutils

def main():
    parser = argparse.ArgumentParser()

    ## optional
    parser.add_argument("-V", "--verbose", help="Print extra information",action="store_true")
    parser.add_argument("-d", "--dimensions", help="Dimensions to split along, separated by commas, with no spaces in between",default = 'nodims')
    parser.add_argument('--skipnan', dest='skipnan', action='store_true')
    parser.add_argument('--no-skipnan', dest='skipnan', action='store_false')
    parser.add_argument("-v", "--variables", help="Variables to split, separated by commas, with no spaces in between. To split all variables, then leave out this argument", default = '__all__')
    parser.set_defaults(skipnan=True)

    ## non-optional
    parser.add_argument('fileA', help = "file A")
    parser.add_argument('fileB', help = "file B")


    args = parser.parse_args()

    skipnan = args.skipnan

    splitDimsOrig = args.dimensions
    if splitDimsOrig == 'nodims':
        splitDims = []
    else:
        splitDims = splitDimsOrig.split(',')

    splitVarsOrig = args.variables
    variable_list = splitVarsOrig.split(',')

    verbose = args.verbose
    files = [args.fileA, args.fileB]
    for infile in files:
        if not os.path.isfile(infile):
            parser.error('Input file does not exist...')

    Cmp = lambda x, y: int(splitDims.index(x) > splitDims.index(y))*2-1

    def format_value(val):
        # if type(val) == type(''):
        #    return '{:<10}'.format('--')
        # else:
        return '{:15.4e}'.format(numpy.float(val))

    ## open the files
    nc = [netCDF4.Dataset(f, 'r') for f in files]

    ## check that the dimensions are the same
    for k in nc[0].dimensions.keys():
        if (k not in nc[1].dimensions.keys()):
            raise Exception("Dimensions don't match by name")
        elif len(nc[0].dimensions[k]) != len(nc[1].dimensions[k]):
            raise Exception("Dimensions don't match by length")
    for k in nc[1].dimensions.keys():
        if (k not in nc[0].dimensions.keys()):
            raise Exception("Dimensions don't match by name")

    ## check that the variables are the same
    for k in nc[0].variables.keys():
        if (k not in nc[1].variables.keys()):
            raise Exception("Variables don't match by name")
        elif nc[0].variables[k].dimensions != nc[1].variables[k].dimensions:
            raise Exception("Variables don't match by size")
    for k in nc[1].variables.keys():
        if (k not in nc[0].variables.keys()):
            raise Exception("Variables don't match by name")

    ## if no variables are specified, apply to all
    if len(variable_list) == 1 and variable_list[0] == '__all__':
        variable_list = nc[0].variables.keys()
    else:
        ## check that all the required variables appear
        for v in variable_list:
            if not (v in nc[0].variables.keys()):
                raise Exception("Variable {} not found in {}".format(v,infile))

    ## check that all the required dimensions appear
    for d in splitDims:
        if not (d in nc[0].dimensions.keys()):
            raise Exception("Dimension {} not found in {}".format(d,fileA))

    ## print the header
    print '{:<30} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}'.format('var','rmse','mean_abs','sd','ratio_sd','maxratio_sd','ratio_mean','maxratio_mean')
    for var in variable_list:
        Type = nc[0].variables[var].dtype.str
        Type = re.sub(r"[<>]", "", Type)

        if Type in ['f4','f8']:
            var0 = nc[0].variables[var][:]
            var1 = nc[1].variables[var][:]
            if len(splitDims) == 0 or (len(Ncpackutils.list_intersection(nc[0].variables[var].dimensions,splitDims)) == 0):
                rmse = numpy.sqrt(numpy.mean(numpy.square(var0 - var1)))
                mean0 = numpy.mean(numpy.fabs(var0))
                sd0 = numpy.std(numpy.fabs(var0))
                # with numpy.errstate(invalid='ignore'):
                with numpy.errstate(invalid='ignore'):
                    ratio_sd = rmse/sd0
                    maxratio_sd = numpy.max(numpy.abs(var0 - var1))/sd0
                    ratio_mean = rmse/mean0
                    maxratio_mean = numpy.max(numpy.abs(var0 - var1))/mean0
            else:
                theseSplitDims = Ncpackutils.list_intersection(nc[0].variables[var].dimensions,splitDims)
                theseSplitDims.sort(Cmp)
                ## get the lengths of these dimensions
                dimLens = [len(nc[0].dimensions[d]) for d in theseSplitDims]
                ## allow for the case of multiple dimensions to split along
                gridRange = tuple([tuple(range(l)) for l in dimLens])
                ## get all combinations of indices along the split dimensions
                splitDimIdxs = list(itertools.product(*gridRange))
                ## get the number of dimensions
                ndim = len(nc[0].variables[var].dimensions)
                ## get the indices of which dimensions for this variable are split
                iSplitDim = numpy.array([ list(nc[0].variables[var].dimensions).index(d) for d in theseSplitDims ])
                iSplitDim.sort()
                gridRange = [range(l) for l in dimLens]
                splitDimIdxs = list(itertools.product(*gridRange))
                indices = numpy.array([Ellipsis] * ndim)
                nsplits = len(splitDimIdxs)

                ## define arrays for the output summaries
                sliceDimLens = tuple([ len(nc[0].dimensions[d]) for d in list(theseSplitDims)])
                RMSE = numpy.ones(sliceDimLens,numpy.float32)
                MEAN = numpy.ones(sliceDimLens,numpy.float32)
                SD = numpy.ones(sliceDimLens,numpy.float32)

                for isplit in range(nsplits):
                    indices[iSplitDim] = splitDimIdxs[isplit]
                    RMSE[tuple(splitDimIdxs[isplit])] = numpy.sqrt(numpy.nanmean(numpy.square(var0[tuple(indices)] - var1[tuple(indices)])))
                    MEAN[tuple(splitDimIdxs[isplit])] = numpy.nanmean(numpy.fabs(var0[tuple(indices)]))
                    SD[tuple(splitDimIdxs[isplit])] = numpy.nanstd(var0[tuple(indices)])
                rmse = numpy.nanmean(RMSE)
                mean0 = numpy.nanmean(numpy.fabs(var0))
                sd0 = numpy.nanmean(SD)
                with numpy.errstate(invalid='ignore'):
                    RATIO_SD = RMSE/SD
                    RATIO_MEAN = RMSE/MEAN
                ##
                if skipnan:
                    ratio_sd = numpy.nanmean(RATIO_SD)
                    maxratio_sd = numpy.nanmax(RATIO_SD)
                    ratio_mean = numpy.nanmean(RATIO_MEAN)
                    maxratio_mean = numpy.nanmax(RATIO_MEAN)
                else:
                    ratio_sd = RATIO_SD.mean()
                    maxratio_sd = RATIO_SD.max()
                    ratio_mean = RATIO_MEAN.mean()
                    maxratio_mean = RATIO_MEAN.max()

            print '{:<30} {} {} {} {} {} {} {}'.format(var,format_value(rmse),format_value(mean0),format_value(sd0),format_value(ratio_sd),format_value(maxratio_sd),format_value(ratio_mean),format_value(maxratio_mean))

