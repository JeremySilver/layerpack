import numpy
import sys
import netCDF4
import datetime

def list_intersection(x,y):
    return list(set(x).intersection(y))

def list_difference(x,y):
    return list(set(x) - set(y))

def pack_to_short(x, max_packed_value = 64000, missing_value = None, packed_missing_value = numpy.uint16(65535)):
    if numpy.isinf(x).any():
        x[numpy.where(numpy.isinf(x))] = numpy.NaN
    
    if missing_value is not None:
        x[numpy.where(x == missing_value)] = numpy.NaN
    ##
    xmin = numpy.nanmin(x)
    xmax = numpy.nanmax(x)
    packed = numpy.ones(x.shape, dtype = numpy.uint16)
    if not (numpy.isfinite(xmin) and numpy.isfinite(xmax)): ## all values are NaN
        scale_factor = 1.0
        add_offset = 0.0
        packed.fill(packed_missing_value)
    elif xmin == xmax: ## range is trivial
        scale_factor = xmin
        add_offset = 0.0
    else:
        add_offset = xmin
        scale_factor = (xmax - xmin)/max_packed_value
        # print 'xmax,xmin,max_packed_value,scale_factor',xmax,xmin,max_packed_value,scale_factor
        packed[:] = numpy.uint16( (x - add_offset) / scale_factor )
        ##
    return (add_offset, scale_factor, packed)

def listVariables(dataset):
    subgroups = dataset.groups
    if len(subgroups) == 0:
        return dataset.variables.keys()
    else:
        result = dataset.variables.keys()
        for subgroupname in subgroups.keys():
            result = result + [subgroups[subgroupname].path + '/' + var for var in listVariables(subgroups[subgroupname]) ]
        return result

def set_pack_history_attribute(ncin, ncout):
    ##
    Args = ''
    for arg in sys.argv:
        Args += (arg + ' ')
    ##
    pack_hist_att = 'pack_history'
    pack_hist_string = '{}: {}'.format(datetime.datetime.now().strftime('%c'), Args)
    if pack_hist_att in ncin.ncattrs():
        pack_hist_string = ncin.getncattr(pack_hist_att) + "; THEN: " + pack_hist_string
    ##
    ncout.setncattr(pack_hist_att,pack_hist_string)
    
