import numpy
import sys
import netCDF4
import datetime
import warnings

def list_intersection(x,y):
    return list(set(x).intersection(y))

def list_difference(x,y):
    return list(set(x) - set(y))

def pack_to_short(x, max_packed_value = 64000.0, missing_value = None, fill_value = None):
    
    packed_nan_value     = numpy.uint16(65531)
    packed_neginf_value  = numpy.uint16(65532)
    packed_posinf_value  = numpy.uint16(65533)
    packed_fill_value    = numpy.uint16(65534)
    packed_missing_value = numpy.uint16(65535)

    warnings.filterwarnings('ignore') ## ignore odd warnings from numpy.isnan(x).any() and others
    
    nan_mask = None
    if numpy.isnan(x).any():
        nan_mask = numpy.where(numpy.isnan(x))

    neginf_mask = None
    if numpy.isneginf(x).any():
        neginf_mask = numpy.where(numpy.isneginf(x))
        x[neginf_mask] = numpy.nan
    
    posinf_mask = None
    if numpy.isposinf(x).any():
        posinf_mask = numpy.where(numpy.isposinf(x))
        x[posinf_mask] = numpy.nan

    warnings.resetwarnings() ## turn it off

    fill_mask = None
    if fill_value is not None:
        if (x == fill_value).any():
            fill_mask = numpy.where(x == fill_value)
            x[fill_mask] = numpy.nan

    missing_mask = None
    if missing_value is not None:
        if (x == missing_value).any():
            missing_mask = numpy.where(x == missing_value)
            x[missing_mask] = numpy.nan
    ##
    xmin = numpy.float64(numpy.nanmin(x))
    xmax = numpy.float64(numpy.nanmax(x))
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
        packed[:] = numpy.uint16( (numpy.float64(x) - add_offset) / scale_factor )

    if nan_mask is not None:
        packed[nan_mask]     = packed_nan_value
    
    if neginf_mask is not None:
        packed[neginf_mask]  = packed_neginf_value
    
    if posinf_mask is not None:
        packed[posinf_mask]  = packed_posinf_value

    if fill_mask is not None:
        packed[fill_mask]    = packed_fill_value

    if missing_mask is not None:
        packed[missing_mask] = packed_missing_value
    ##
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
    
