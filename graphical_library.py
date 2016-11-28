
import scipy.io
import numpy
import numpy.random as numran
import bisect

import matplotlib
# set rendering
#matplotlib.use('Agg')
import matplotlib.pylab
#import matplotlib.numerix
import matplotlib.colors
import matplotlib.cm
import pylab
import matplotlib.pyplot as plt

import time_utilities as timelib
import matplotlib.dates as matdates
import matplotlib.ticker as matticker
#import mpl_toolkits.axis_artist as mpltaa
import mpl_toolkits.axes_grid.axislines as mpltagal


def plot_stack(x2d, y2d, x2drange=None, y2drange=None, xlabel=None, \
 ylabel=None, xmajtick=None, nxmintick=None, ymajtick=None, ymintick=None, \
 ybaseline=None, title=None, nrows=None, x2drange_offset=None):

 if x2d is None:
  nt = 24; nd = 35
  dstart = timelib.totime(2003,11,1)
  x2d = numpy.tile(numpy.nan, (nt, nd))
  for j in range(nd):
   x2d[:, j] = numpy.arange(dstart + 2 * j, dstart + 2 * j + 1, 1. / 24.)
   
 if y2d is None: y2d = numran.rand(nt, nd) - 0.5

 nplots, dummy = numpy.transpose(x2d).shape
 if nrows is None: nrows = 15; 
 #ncols = numpy.int(numpy.ceil(nplots / nrows))
 ncols = numpy.int(numpy.ceil(nplots / float(nrows)))
 #nrows, dummy = numpy.transpose(x2d).shape; ncols = 1;
 
 #print nplots, ncols, nrows
 #print x2d.shape, y2d.shape

 nx = ncols; ny = nrows; #fdx = 0.95; fdy = 0.85;
 fdx = .9; fdy = 1.;
 x0, dx, yf, dy = axis_coord(nx, ny)

 if x2drange is None: x2drange = numpy.array([x2d[0, 0], x2d[-1, 0]])
 if y2drange is None: 
  y2drange = numpy.array([numpy.floor(numpy.min(y2d)), numpy.ceil(numpy.max(y2d))])
 if x2drange_offset is None: x2drange_offset = numpy.array([0., 23.99])

 if xlabel is None: xlabel = 'LT'
 if ylabel is None: ylabel = 'Unit'
 if title is None: title = 'Parameter'
 if ybaseline is None: ybaseline = 0.

 if xmajtick is None: xmajtick = 10
 if nxmintick is None: nxmintick = 2
 if ymajtick is None: ymajtick = 10

 #xmajor_tbin = timelib.dtbin([numpy.floor(x2drange[0]), numpy.ceil(x2drange[1])], xmajtick*60*60)
 #ind = numpy.where((xmajor_tbin >= x2drange[0]) & (xmajor_tbin <= x2drange[1]))
 #xmajor_ticks = xmajor_tbin[ind]

 ylim = 10. * numpy.array([numpy.floor(y2drange[0] / 10.), numpy.ceil(y2drange[1] / 10.)])
 ytick = numpy.arange(ylim[0], ylim[1] + ymajtick, ymajtick)
 indyt = numpy.where((ytick >= y2drange[0]) & (ytick <= y2drange[1]))
 if len(indyt[0] > 0): ymajor_ticks = ytick[indyt[0]]

 # Graphical routines
 xsize = (15 if ncols > 2 else 10 * ncols)
 ysize = (15 if nrows > 3 else 3 * nrows)
 figid = plt.figure(num=20, figsize=(xsize, ysize))
 figid.clf()

 counter = 0; ispanel = 0
  
 ######
 #          
 # From left to right, then up to down, ...
 #
 #for j in range(ny):
 #               
 # for i in range(nx):
 #
 #######
 #
 # From up to down, then left to right
 #
 for i in range(nx):
 #
  for j in range(ny):               
 #                   
   if (counter <= (nx * ny) - 1):
    
    ispanel = 1

    #if ispanel == 1:
    if counter < nplots:

     # in normal coordinates
     pn = plt.axes([x0 + i * dx, yf - (j + 1) * dy, fdx * dx, fdy * dy])
 
     # Plot commands here ....
     #

     # Defining "box" limits
     x2d_range_dummy = [numpy.floor(numpy.min(x2d[:, counter])), numpy.ceil(numpy.max(x2d[:, counter]))]
     toffset = numpy.array(x2drange_offset) / 24.
     x2d_range_dummy = [x2d_range_dummy[0] + toffset[0], x2d_range_dummy[1] - (1. - toffset[1])]
     x2d_range = x2d_range_dummy
     #
     plt.plot(x2d[:, counter], y2d[:, counter], linestyle='-', color='k')
     plt.plot(x2d_range, numpy.tile(ybaseline, 2), linestyle=':', color='k')
     plt.xlim(x2d_range[0], x2d_range[1])
     plt.ylim(y2drange[0], y2drange[1])

     xmajor_tbin = timelib.dtbin([numpy.floor(x2d_range[0]), numpy.ceil(x2d_range[1])], xmajtick*60*60)
     ind = numpy.where((xmajor_tbin >= x2d_range[0]) & (xmajor_tbin <= x2d_range[1]))
     xmajor_ticks = xmajor_tbin[ind]

     ## minorticks
     #pn.xaxis.set_minor_locator(plt.MultipleLocator(base=xmajtick / 24. / nxmintick))

     plt.xlabel(xlabel); plt.ylabel(ylabel)
     #plt.ylabel(matdates.num2date(x2d_range[0]).strftime('%m/%d/%Y'), rotation='horizontal')
     plt.text(0.8, 0.75, matdates.num2date(x2d_range[0]).strftime('%m/%d/%Y'), \
      horizontalalignment='center', verticalalignment='center', \
      transform = pn.transAxes, bbox=dict(facecolor='white', alpha=0.5))

     pn.axes.set_yticks(ymajor_ticks)

     if i >= 0 and j == 0: # Top panel

      plt.title(title)

      pn.axes.set_xticks(xmajor_ticks)          
      pn.xaxis.set_major_formatter(matdates.DateFormatter('%H:%M'))

      pn.spines['bottom'].set_color('none')

      clear_ticklabels(plt.gca(), 1, 1)

      for ii, line in enumerate(pn.get_xticklines()): # Delete ticks (bottom x-axis)
       if not ii % 2 == 1:   # even indices
        line.set_visible(False)

     elif i >= 0 and (j > 0 and j < (ny - 1)): # Intermediate panels      

      if counter == (nplots - 1):
       clear_ticklabels(plt.gca(), 0, 1)
       pn.axes.set_xticks(xmajor_ticks)          
       pn.xaxis.set_major_formatter(matdates.DateFormatter('%H:%M'))
       for ii, line in enumerate(pn.get_xticklines()): # Delete ticks (top x-axis)
        if ii % 2 == 1:   # odd indices
         line.set_visible(False)
      else:
       pn.spines['bottom'].set_color('none')
       clear_ticklabels(plt.gca(), 1, 1)       
       pn.xaxis.set_major_locator(plt.NullLocator())
       pn.xaxis.set_minor_locator(plt.NullLocator())

      pn.spines['top'].set_color('none')

     elif i >= 0 and j == (ny - 1):	# Bottom panel

      pn.axes.set_xticks(xmajor_ticks)          
      pn.xaxis.set_major_formatter(matdates.DateFormatter('%H:%M'))

      pn.spines['top'].set_color('none')	

      #for ii, line in enumerate(pn.get_xticklines() + pn.get_yticklines()): # Delete ticks (top x-axis and right y-axis)
      for ii, line in enumerate(pn.get_xticklines()): # Delete ticks (top x-axis)
       if ii % 2 == 1:   # odd indices
        line.set_visible(False)
                        
      if i > 0:
       clear_ticklabels(plt.gca(), 0, 1)
      else:
       clear_ticklabels(plt.gca(), 0, 0)

     #elif i > 0 and j == (ny - 1):
     #elif i > 0 and j < (ny - 1):
      
     # print 'Ronald'
          

   counter += 1
   ispanel = 0


def plot(x, y, fignum=0, figsize=(8, 6), graph=0, title='Title', xlabel='X-axis', xlog=0, \
         ylabel='Y-axis', ylog=0):
    
    if (fignum > 0):
        pylab.figure(fignum, figsize=figsize)
        pylab.add_subplot(1, 1, 1)
    #else:
    #    pylab.axes([0 0 0 0])
    
    if ((xlog == 1) & (ylog == 0)):
        pylab.semilogx(x, y, color=color, linestyle=linestyle)
    elif ((xlog == 0) & (ylog == 1)):
        pylab.semilogy(x, y, color=color, linestyle=linestyle)
    elif ((xlog == 1) & (ylog == 1)):
        pylab.loglog(x, y, color=color, linestyle=linestyle)
    else:
        pylab.plot(x, y, color=color, linestyle=linestyle)
        
    pylab.title(title); pylab.xlabel(xlabel); pylab.ylabel(ylabel)
        
    
def read_maghcomp_mat(dpath, year, month):
    filename = dpath + ('jic-mag-%04i-%02i.mat') % (year, month)
    data = scipy.io.loadmat(filename)
    nx, ny = data.get('time').shape
    matdata = {'time':data.get('time'), 'H':data.get('H')}
    vecdata = {'time':data.get('time').reshape(nx * ny).shape, \
               'H':data.get('H').reshape(nx * ny).shape}
    return({'matdata':matdata, 'vecdata':vecdata})


def add_ticklabels(pn):

 curr_xticklabels = pn.axes.get_xticklabels()
 nxticks = len(curr_xticklabels)
 xticklabels = []

 for i in range(nxticks):
  xticklabels.append('dotPE')                

 pn.xaxis.set_ticklabels(xticklabels)

#
# End of 'add_ticklabels'
#####


def clear_ticklabels(pn, in_x, in_y):
    
    if in_x > 0:
        curr_xticklabels = pn.axes.get_xticklabels()
        nxticks = len(curr_xticklabels)
        xticklabels = [];
        for i in range(nxticks):
            xticklabels.append('')                
            pn.xaxis.set_ticklabels(xticklabels)
        pn.set_xlabel('')
    
    if in_y > 0:
        curr_yticklabels = pn.axes.get_yticklabels()
        nyticks = len(curr_yticklabels)
        yticklabels = [];
        for i in range(nyticks):
            yticklabels.append('')                
            pn.yaxis.set_ticklabels(yticklabels)
        pn.set_ylabel('')
    
    return(0)

def clear_spines(pn, in_x, in_y):

 if in_x[0] > 0: 
  #pn.xaxis["bottom"].set_visible(False)
  #pn.spines['bottom'].set_visible(False)
  pn.spines['bottom'].set_color('none')
  #pn.xaxis.set_ticks_position('bottom')
  #pn.xaxis.tick_bottom()
 if in_x[1] > 0: pn.spines['top'].set_color('none')
 #if in_x[0] > 0 or in_x[1] > 0:
  #pn.xaxis.set_ticks([])
  #pn.xaxis.set_minor_locator(matticker.NullLocator())

 if in_y[0] > 0: pn.spines['left'].set_color('none')
 if in_y[1] > 0: pn.spines['right'].set_color('none')
 #if in_y[0] > 0 or in_y[1] > 0:
  #pn.yaxis.set_ticks([])
  #pn.yaxis.set_minor_locator(matticker.NullLocator()) 

#def clear_ticks(pn, in_x, in_y):
##def clear_axis(pn, in_x, in_y):

## if in_x[0] > 0: pn.axis['bottom'].set_visible(False)
## if in_x[1] > 0: pn.axis['top'].set_visible(False)

## if in_y[0] > 0: pn.axis['left'].set_visible(False)
## if in_y[1] > 0: pn.axis['right'].set_visible(False)

#
# End of ''
#####


def insert_gap(z, x):

 #
 # Inserts "gaps" on 2-D array data if time interval between two consecutive
 # data points is greater than the time step used on the processing to #
 # generate the original data set
 #

 print('Inserting gap ...')

 # first find the time interval greater than 90% of others
 #
 #timeIntervalList = []
 #for i in range(len(x) - 1):
 # timeIntervalList.append(x[i + 1] - x[i])
 dx = numpy.diff(x); timeIntervalList = numpy.append(dx, dx[-1]).tolist()

 timeIntervalList.sort()
 index = int(0.9 * len(timeIntervalList))
 maxInterval = timeIntervalList[index]

 nx, ny = z.shape	# get dimensional size of 2-D array data

 # Filling up the new time series vector and 2-D array
 #
 for i in range(len(x) - 1):

  if i == 0:

   zout = numpy.resize(z[i, :], (1, ny)); xout = x[i];

  else:

   if x[i + 1] - x[i] > maxInterval:	# if condition is satisfied, add a row of nan values

    zout = numpy.append(zout, numpy.tile(numpy.nan, (1, ny)), axis=0)
    xout = numpy.append(xout, numpy.mean([x[i], x[i + 1]]))	# arithmetic mean of two closest time series values

   else:

    zout = numpy.append(zout, numpy.resize(z[i, :], (1, ny)), axis=0)
    xout = numpy.append(xout, x[i])

 zout = numpy.append(zout, numpy.resize(z[nx - 1, :], (1, ny)), axis=0)
 xout = numpy.append(xout, x[nx - 1])

 return(zout, xout)

#
# End of 'insert_gap'
#####    


def pcolorplot(x, y, z, minColormap = None, maxColormap = None, \
               startTime = None, endTime = None, yMinimum = None, yMaximum = None, \
               maxNumAlt = None, smoothAltitude = True, \
               insertDataGap = 5, colorMap = matplotlib.cm.jet, xLabelStr = '', \
               yLabelStr = '', titleStr = '', useAbsoluteTime = False, \
               altYTitle = None, altYLabels = None, \
               fullFilename = '/home/rilma/Desktop/mad-pcolor.png'):

#    minColormap = None; maxColormap = None;
#    startTime = None; endTime = None;
#    yMinimum = None; yMaximum = None;
#    maxNumAlt = None; smoothAltitude = True;
#    insertDataGap = 5; colorMap = matplotlib.cm.jet;
#    xLabelStr = ''; yLabelStr = ''; titleStr = '';
#    useAbsoluteTime = False; altYTitle = None; altYLabels = None;
#    fullFilename = '/home/rilma/Desktop/mad-pcolor.png';

    nx, ny = z.shape;
    array_data = numpy.ones((nx * ny, 3)) * numpy.nan;
    counter = 0;
    
    xmat = x; ymat = y;

    if len(x.shape) == 1: xmat = x.reshape((len(x), 1)) * numpy.ones((1, ny))
    xvec = xmat.reshape(nx * ny);    

    if len(y.shape) == 1: ymat = numpy.ones((nx, 1)) * y.reshape((1, len(y)))
    yvec = ymat.reshape(nx * ny);

    zvec = z.reshape(nx * ny);

    array_data[:, 0] = xvec; array_data[:, 1] = yvec; array_data[:, 2] = zvec;

    missing = 1.0E30 # special value, used to create a masked array
    parameter_count = 3

    # The first pass through is to obtain the number and range of the  x and y variables
        
    xList = []; yList = []; zList = [];
    zMin = None; zMax = None;

    # loop over all the lines of data in the array
    for j in range(len(array_data)):

        try:
                
            x = array_data[j][0]; y = array_data[j][1]; z = array_data[j][2]

            if x not in xList:
                xList.append(x)
        
            if y not in yList:
                yList.append(y)

            if z != missing:
                zList.append(z)

            if zMin != None:
                if z < zMin and z != missing:
                    zMin = z
            elif z != missing:
                zMin = z

            if zMax != None:
                if z > zMax and z != missing:
                    zMax = z
            elif z != missing:
                zMax = z

        except:
            continue
            

    
    #print 'Ronald!'
    

    if zMin == None:
        raise ValueError('No valid z data found')

    # if both minColormap and maxColormap == None, use autoscaling
    if minColormap == None and maxColormap == None:
        zList.sort()
        d10 = zList[int(len(zList)*0.10)]
        d90 = zList[int(len(zList)*0.90)]

        zMin = d10 - (d90-d10) * 0.75
        zMax = d90 + (d90-d10) * 0.75
    
    # now sort the X and Y axis lists and pull their length
        
    xList.sort()
    if startTime == None:
        xMin = xList[0]
    else:
        xMin = startTime
    if endTime == None:
        xMax = xList[-1]
    else:
        xMax = endTime
    yList.sort()
    max_x_dimension = len(xList)
    max_y_dimension = len(yList)

    if yMinimum == None:
        yMinimum = yList[0]
    if yMaximum == None:
        yMaximum = yList[-1]
    
    truncateAlt = False
    if maxNumAlt != None:
        if max_y_dimension > maxNumAlt:
            truncateAlt = True

    # build dictonary of indexes into xList
    xListDict = {}
    for i in range(len(xList)):
        xListDict[xList[i]] = i

    # if self.truncateAlt == False, build dictonary of indexes into yList,
    # else truncate y values by builing a list of maxNumAlt ranges
    if truncateAlt == False:
        yListDict = {}
        for i in range(len(yList)):
            yListDict[yList[i]] = i
    else:
        yListRanges = []
        for i in range(maxNumAlt):
            yListRanges.append(yList[int(i*(len(yList)/maxNumAlt))])
        max_y_dimension = maxNumAlt
    
    # now build arrays to handle the X axis label, Y axis label, and the Z data

    #X = matplotlib.numerix.zeros((max_x_dimension, max_y_dimension), matplotlib.numerix.Float32)
    #Y = matplotlib.numerix.zeros((max_x_dimension, max_y_dimension), matplotlib.numerix.Float32)
    #Z = matplotlib.numerix.ones((max_x_dimension, max_y_dimension), matplotlib.numerix.Float32)
    
    X = numpy.zeros((max_x_dimension, max_y_dimension), dtype=numpy.float32)
    Y = numpy.zeros((max_x_dimension, max_y_dimension), dtype=numpy.float32)
    Z = numpy.ones((max_x_dimension, max_y_dimension), dtype=numpy.float32)

    # all parameter values default to missing
    Z = Z * missing

    # fill the X and Y arrays

    for i in range(max_x_dimension):
        for j in range(max_y_dimension):
            X[i][j] = float(xList[i])
            if truncateAlt:
                Y[i][j] = float(yList[int(j*(len(yList)/maxNumAlt))])
            else:
                Y[i][j] = float(yList[j])    
    
        # Now load up the data array Z with the array_data measurements as a function of x and y
        previousIndex = None
        previousValue = None
        presentTime = None
        newTimeFound = True

    for k in range(len(array_data)):
        try:

            xdata = array_data[k][0]
            ydata = array_data[k][1]
            zdata = array_data[k][2]

            if zdata == missing:
                continue

            if xdata != presentTime:
                newTimeFound = True
            else:
                newTimeFound = False
            presentTime = xdata

            # now find the right place in the array for this data point
            i = xListDict[xdata]
            #j = getYIndex(ydata)
            #j = getYIndex(truncateAlt, yListDict, yListRanges, ydata)
            j = getYIndex(truncateAlt, yListDict, [], ydata)

            Z[i][j] = zdata
		
            # now see if we need to fill in any gaps
            if (not newTimeFound) and smoothAltitude:
                if previousIndex < j - 1:
                    # fill in all missing points
                    for k in range(previousIndex + 1, j):
                        # simply average between the points based on index
                        thisValue = previousValue + (zdata - previousValue)*(float(k-previousIndex)/float(j-previousIndex))
                        Z[i][k] = thisValue
			   
            previousIndex = j
            previousValue = zdata
		
                
        except:
            continue
    
    # insert missing data to represent gaps if needed
    if insertDataGap != None:
        # first find the time interval greater than 90% of others
        timeIntervalList = []
        for i in range(len(xList) - 1):
            timeIntervalList.append(xList[i+1] - xList[i])
        timeIntervalList.sort()
        index = int(len(timeIntervalList)*0.9)
        maxInterval = timeIntervalList[index]
	    
        for i in range(len(xList) - 1):
            if xList[i+1] - xList[i] > maxInterval * insertDataGap:
                Z[i,:] = missing


    # set up plotting parameters
    if minColormap == None:
        minColormap = zMin
    if maxColormap == None:
        maxColormap = zMax

    matplotlib.pylab.pcolor(X,Y,Z, shading='flat', vmin=minColormap, \
     vmax=maxColormap, cmap = colorMap, norm = matplotlib.pylab.normalize())
    matplotlib.pylab.colorbar()

    matplotlib.pylab.xlabel(xLabelStr)
    matplotlib.pylab.ylabel(yLabelStr)
    matplotlib.pylab.xlim(xMin, xMax)
    matplotlib.pylab.ylim(yMinimum, yMaximum)
    matplotlib.pylab.yticks()
    matplotlib.pylab.title(titleStr)

    if useAbsoluteTime:
        locs, labels = matplotlib.pylab.xticks()
        if len(locs) > 5:
            # truncate by len(locs) / 5
            scale = 1 + int(len(locs) / 5)
            new_locs = []
            for i in range(len(locs)):
                if i % scale == 0:
                    new_locs.append(locs[i])
            locs = new_locs
        newXTicks = convertToAbsoluteTimeStr(locs)
        matplotlib.pylab.xticks(locs, newXTicks, rotation=15)

    matplotlib.pylab.xticks()

    # add second y-axis if desired
    if altYTitle != None and altYLabels != None:
        ax2 = matplotlib.pylab.twinx()
        matplotlib.pylab.ylabel(altYTitle)
        ax2.yaxis.tick_right()
        matplotlib.pylab.yticks(range(len(altYLabels)), altYLabels)
        matplotlib.pylab.yticks()

    # get the handle to the figure now that it has been created so we can manipulate on subsequent calls.
        
    #self.__figure = matplotlib.pylab.gcf()
          
    matplotlib.pylab.savefig(fullFilename)
    #matplotlib.pylab.show()

    #matplotlib.pylab.clf()
        
    #print 'Stop!'

#
# End of ''
#####


def getYIndex(truncateAlt, yListDict, yListRanges, yvalue):
    """__getYIndex returns the correct index into the y dimension for a given y value.

    Input: yvalue - value of y parameter

    Returns: the correct index into the y dimension

    Algorithm: if self.truncateAlt == False, simple use the dictionary self.yListDict.  Else
    loop through self.yListRanges and return the first greater than the requested value
    """
    if truncateAlt == False:
        return yListDict[yvalue]
    else:
        i = bisect.bisect_left(yListRanges, yvalue)
        if i >= len(yListRanges):
            i = len(yListRanges) - 1
        return i
        
    
def axis_coord(nx, ny):
    
    xmin = 0.0; xmax = 1.0; ix = 0.05 * 1.5;
    x0 = xmin + ix; xf = xmax - ix; dx = (xf - x0) / nx
    
    ymin = 0.0; ymax = 1.0; iyb = 0.05 * 2; iyt = 0.10;
    y0 = ymin + iyb; yf = ymax - iyt; dy = (yf - y0) / ny    
    
    return(x0, dx, yf, dy)

def save_figure(svfig, orientation_value=None):
 #
 # Saving panel figure content on a file
 #
 figfname = svfig.get('figname')
 format_type = svfig.get('format')
 dpi_value = svfig.get('dpi')
 if orientation_value is None: orientation_value = 'portrait'

 matplotlib.pylab.savefig(figfname, format=format_type, dpi=dpi_value, \
  bbox_inches='tight', orientation=orientation_value)

 print('Saving figure file --> %s' % figfname)
 #



def get_majtick(valrange, majtick):

 #
 # calculate the values written as "major ticks" in plot axes.
 #

 lim = 100. * pylab.array([pylab.floor(valrange[0] / 100.), \
  pylab.ceil(valrange[1] / 100.)])
 tick = pylab.arange(lim[0], lim[1] + majtick, majtick)
 indt = pylab.where((tick >= valrange[0]) & (tick <= valrange[1]))
 if len(indt[0] > 0): major_ticks = tick[indt[0]]

 return(major_ticks)

# End of 'get_majtick'


#def reverse_cmap(cmap, name='colormap', N=256):
# #
# # Reversing a "colormap"
# #
# # INPUTS:
# #	cmap:
# #	name:
# #	N: Number of rgb quatization levels
# #
# cmap_dict = cmap._segmentdata
# cmap_output = {}
# dummy1 = []
# for key in ('red','green','blue'):
#  cmap_key = numpy.array(cmap_dict[key])
#  nlev, ncomp = cmap_key.shape
#  dummy = cmap_key[range(nlev - 1, 0 - 1, -1), :]
#  for j in range(nlev):
#   dummy1.append(tuple(dummy[j, :]))
#  cmap_output[key] = list(dummy1)
#  dummy1 = []
#
# return matplotlib.colors.LinearSegmentedColormap(name, cmap_output, N)

def main():
    
    #year = 1997; month = 03;
    ##dpath = '/media/sda5/database/magnetometer/jic97/jro-hcomp/'
    #dpath = '/home/rilma/Desktop/elnino-1997/jro-hcomp/'
    #data = read_maghcomp_mat(dpath, year, month)
    
    #dtime = data.get('matdata').get('time');
    #x = dtime[0, :]; y = dtime[:, 0]; 
    #z = numpy.transpose(data.get('matdata').get('H'));
    
 x = numpy.linspace(0, 254, 255); y = x;
 z = numpy.clip(numran.randn(255, 255), -1, 1)
  
 #print x.shape, y.shape, z.shape  
 #print len(x.shape)
 ##pcolorplot(x, y, z)
 pcolorplot(x, y, z, -1, 1)
 pylab.show()
    
def test_stack_plot():
 x2d = None; y2d = None;
 xmajtick = 4; nxmintick = 8;
 plot_stack(x2d, y2d, x2drange=None, y2drange=None, xmajtick=xmajtick, nxmintick=nxmintick)
 figname = '/home/rilma/Desktop/test-stack.png'
 svfig = {'save' : 1, 'show' : 1, 'dpi' : 300, 'figname' : figname, 'format' : 'png'}
 save_figure(svfig)
 plt.show()

if __name__ == '__main__':
 main()
 #test_stack_plot()
    
