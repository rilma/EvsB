
import pylab


def amax(x, axis=None):

 x_masked = pylab.ma.masked_where(pylab.isnan(x), x)

 return pylab.ma.amax(x_masked, axis=axis)
# End of 'amax'


def amin(x, axis=None):

 x_masked = pylab.ma.masked_where(pylab.isnan(x), x)

 return pylab.ma.amin(x_masked, axis=axis)
# End of 'amin'


def mean(x, axis=None, dtype=None):

 x_masked = pylab.ma.masked_where(pylab.isnan(x), x)
 
 return pylab.ma.mean(x_masked, axis=axis, dtype=dtype) 
# End of 'mean'


def ts_binned(data, xstep, mtype):

 import time_utilities as time_util

 time = data['time']; values = data['values']; 
 dtrange = pylab.floor(time[[0, len(time)-1]]) + pylab.array([0, 1.0])    
 dtbins = time_util.dtbin(dtrange, xstep)
 bin_values = pylab.ones(len(dtbins)) * pylab.nan
 for i in range(len(dtbins)):
  curr_dtrange = dtbins[i] + pylab.array([-1, 1]) * 0.5 * xstep / (24.0 * 3600.0)
  ind1 = pylab.where((time >= curr_dtrange[0]) & (time < curr_dtrange[1]))
  if (len(ind1[0]) > 0):
   ind2 = pylab.where(pylab.isfinite(values[ind1[0]]))
   if (len(ind2[0]) > 0):
    if (mtype == 1):     # Mean
     curr_output = pylab.mean((values[ind1[0]])[ind2[0]])
    elif (mtype == 2):   # Weigthed-Mean
     print('Weigthed-Mean')
    else:   # None
     curr_output = pylab.nan                    
   else:
    curr_output = pylab.nan
  else:
   curr_output = pylab.nan
  bin_values[i] = curr_output
        
 output = {'time' : dtbins, 'values' : bin_values}
 return(output)
# End of 'ts_binned'


def corrcoef(var1, var2, trange, detrend=None):

 if detrend is None: detrend = False

 var1['time'] = pylab.squeeze(var1['time']); var1['value'] = pylab.squeeze(var1['value'])
 var2['time'] = pylab.squeeze(var2['time']); var2['value'] = pylab.squeeze(var2['value'])

 ind1 = pylab.where((var1['time'] >= trange[0]) & (var1['time'] <= trange[1]))
 if len(ind1[0]) > 0:
  var1['time'] = var1['time'][ind1[0]]; var1['value'] = var1['value'][ind1[0]]

 ind2 = pylab.where((var2['time'] >= trange[0]) & (var2['time'] <= trange[1]))
 if len(ind2[0]) > 0:
  var2['time'] = var2['time'][ind2[0]]; var2['value'] = var2['value'][ind2[0]]


 var1_masked = pylab.ma.masked_where(pylab.isnan(var1['value']), var1['value'])
 var2_masked = pylab.ma.masked_where(pylab.isnan(var2['value']), var2['value'])

 #
 # Detrending
 #

 if detrend:

  var1_detrended = (var1_masked - pylab.ma.mean(var1_masked)) / \
   pylab.ma.std(var1_masked)
  var1_value = var1_detrended

  var2_detrended = (var2_masked - pylab.ma.mean(var2_masked)) / \
   pylab.ma.std(var2_masked)
  var2_value = var2_detrended

 else:

  var1_value = var1_masked; var2_value = var2_masked;

 #

 # Correlation coefficients
 corr_coef = pylab.ma.corrcoef(var1_value, var2_value)

 return(corr_coef)
# End of 'corrcoef'


