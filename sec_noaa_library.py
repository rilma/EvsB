
import pylab
import time_utilities as time_util
import graphical_library as graph_lib

pylab.rcParams['xtick.direction'] = 'out'
pylab.rcParams['ytick.direction'] = 'out'

## Tuning the Figure setup (Tex)
##
pylab.rc('font',**{'family':'serif', 'serif':['New Century Schoolbook']})
pylab.rc('text', usetex=True)
##   


def calc_hrows(filename):

 rows_counter = 0

 fid = open(filename, 'r')

 while 1:

  strline = fid.readline()

  # Checking EOF
  if not strline: break

  strmark = '#---------'
  if strline[0:len(strmark)] == strmark:
   fid.close()
   return(rows_counter)

  rows_counter += 1
 
 fid.close()

 return(rows_counter)


def read_lists(date, dpath=None, satname=None, inst=None):

 import os

 if satname is None: satname = 'ace'
 if inst is None: inst = ['mag', 'swepam'] 
 if dpath is None: dpath = '/home/rilma/database/sec-noaa/lists/'
 curr_dpath = dpath + satname + os.sep

 for i in range(len(inst)):

  if inst[i] == 'mag':

   filename = ('%04i%02i%02i' % (date[2], date[0], date[1])) + '_ace_mag_1m.txt' 
   fname = curr_dpath + filename
   print( 'Reading: %s' % fname )
   skiprows = calc_hrows(fname)
   mag_data = pylab.loadtxt(fname, skiprows=skiprows)

   #yr, mo, da, hhmm, mjday, secday, S, bx, by, bz, bt, lat, lon = \
   # pylab.loadtxt(fname, skiprows=skiprows, unpack=True)
   #hh = 100 * pylab.int32(pylab.floor(0.01 * hhmm))
   #mm = hmm - 100 * hh
   #time = time_util.totime(yr, mo, da, hh, mm, 0, 0)

  #else: mag_data = []

  if inst[i] == 'swepam':

   filename = ('%04i%02i%02i' % (date[2], date[0], date[1])) + '_ace_swepam_1m.txt'
   fname = curr_dpath + filename
   print( 'Reading: %s' % fname )
   skiprows = calc_hrows(fname)
   swepam_data = pylab.loadtxt(fname, skiprows=skiprows)

  #else: swepam_data = []   

 return({'mag' : mag_data, 'swepam' : swepam_data})


def read_all_lists(trange, dpath=None, inst=None, satname=None):

 ndays = int(trange[1] - trange[0]) + 1

 for i in range(ndays):

  curr_trange = trange[0] + pylab.array([0., 0.999]) + i
  curr_date = pylab.floor(curr_trange[0])
  dt = pylab.num2date(curr_date)
  date = [dt.month, dt.day, dt.year]
  currdt_data = read_lists(date, dpath=dpath, satname=satname, inst=inst)

  if i == 0:
   mag_data = currdt_data['mag']; swepam_data = currdt_data['swepam']
  else:
   if len(mag_data) > 0:
    if len(currdt_data['mag']) > 0:
     mag_data = pylab.append(mag_data, currdt_data['mag'], axis=0)
   if len(swepam_data) > 0:
    if len(currdt_data['swepam']) > 0:
     swepam_data = pylab.append(swepam_data, currdt_data['swepam'], axis=0)

 yr, mo, da, hhmm, mjday, secday, S, bx, by, bz, bt, lat, lon = mag_data.T
 hh = pylab.int32(pylab.floor(0.01 * hhmm))
 mm = hhmm - 100 * hh
 mag_time = time_util.totime(yr, mo, da, hh, mm, 0, 0)
 mag = {'time' : mag_time, 'S' : S, 'bx' : clean_1d_list(bx, S), \
 'by' : clean_1d_list(by, S), 'bz' : clean_1d_list(bz, S), \
  'lat' : lat, 'lon' : lon}

 yr, mo, da, hhmm, mjday, secday, S, swpd, swbs, swit = swepam_data.T
 hh = pylab.int32(pylab.floor(0.01 * hhmm))
 mm = hhmm - 100 * hh
 swepam_time = time_util.totime(yr, mo, da, hh, mm, 0, 0)
 swepam = {'time' : swepam_time, 'S' : S, 'swpd' : clean_1d_list(swpd, S), \
  'swbs' : clean_1d_list(swbs, S), 'swit' : clean_1d_list(swit, S)}

 return({'mag' : mag, 'swepam' : swepam})


def clean_1d_list(input_data, status):
 ind = pylab.where(status > 0)
 if len(ind[0]) > 0: input_data[ind[0]] = pylab.nan
 return(input_data)


def plot_lists(data, ipar=None, ntmintick=None, trange=None, \
 tmajtick=None, tlabel=None):
 
 if ipar is None: ipar = pylab.array([[1, 11], [2, 12], [3, 13]])
 nrows, ncols = ipar.shape
 nx = ncols; ny = nrows; fdx = 0.85; fdy = 0.85;
 x0, dx, yf, dy = graph_lib.axis_coord(nx, ny)

 time = data['mag']['time']
 
 tinfo = time;

 if trange == None: trange = pylab.array([tinfo[0], max(tinfo)])

 indt = pylab.where((tinfo >= trange[0]) & (tinfo <= trange[1]))
 indt = indt[0] if len(indt[0]) > 0 else pylab.array(range(len(tinfo)))
 
 dtmin_obj = min(trange); dtmax_obj = max(trange);
    
 xmajor_tbin = time_util.dtbin([pylab.floor(dtmin_obj), pylab.ceil(dtmax_obj)], tmajtick*60*60)
 ind = pylab.where((xmajor_tbin >= dtmin_obj) & (xmajor_tbin <= dtmax_obj))
 xmajor_ticks = xmajor_tbin[ind]
 
 # Graphical routines
 figid = pylab.figure(num=41, figsize=(10 * ncols, 3 * nrows))
 figid.clf()

 counter = 0; ispanel = 0

 for j in range(ny):
                
  for i in range(nx):
                    
   if (counter <= (nx * ny) - 1):
    if ipar[j, i] == 1: # Bx
     xvalues = data['mag']['time']; yvalues = data['mag']['bx'];
     vmin = -40.; vmax = 40.; title = '$B_x$';
     title_units = 'nT'; ispanel = 1
    if ipar[j, i] == 2: # By
     xvalues = data['mag']['time']; yvalues = data['mag']['by'];
     vmin = -40.; vmax = 40.; title = '$B_y$';
     title_units = 'nT'; ispanel = 1
    if ipar[j, i] == 3: # Bz
     xvalues = data['mag']['time']; yvalues = data['mag']['bz'];
     vmin = -40.; vmax = 40.; title = '$B_z$';
     title_units = 'nT'; ispanel = 1
    if ipar[j, i] == 11: # Np
     xvalues = data['swepam']['time']; yvalues = data['swepam']['swpd'];
     vmin = 0.; vmax = 35.; title = 'Proton Density';
     title_units = 'cm$^{-3}$'; ispanel = 1
    if ipar[j, i] == 12: # Vsw
     xvalues = data['swepam']['time']; yvalues = data['swepam']['swbs'];
     vmin = 300.; vmax = 700.; title = 'Bulk Speed';
     title_units = 'km/s'; ispanel = 1
    if ipar[j, i] == 13: # Ti
     xvalues = data['swepam']['time']; yvalues = data['swepam']['swit'];
     vmin = 0.; vmax = 1e6; title = 'Ion Temperature';
     title_units = '$^\circ$K'; ispanel = 1

    if ispanel == 1:

     # in normal coordinates
     pn = pylab.axes([x0 + i * dx, yf - (j + 1) * dy, fdx * dx, fdy * dy])     

     #
     # Using "plot"
     #
     x = xvalues[indt]; y = yvalues[indt];
     ipn = pylab.plot(x, y, linestyle='-', clip_on=True, color='k')
     ipn0 = pylab.plot(trange, pylab.zeros(2), linestyle=':', color='k')
     # Defining "box" limits
     pylab.xlim(trange[0], trange[1])
     pylab.ylim(vmin, vmax)
     #

     pn.axes.set_xticks(xmajor_ticks);
     pn.axes.minorticks_on()

     pn0 = pylab.gca()
     pn0.xaxis.set_minor_locator(pylab.MultipleLocator(base=tmajtick / 24. / ntmintick))
     #pn0.yaxis.set_minor_locator(pylab.MultipleLocator(base=hmintick))

     # "Add minorticks" here ...        
     pn.xaxis.set_major_formatter(pylab.DateFormatter("%H:%M"))        

  ##plt.tick_params(axis='both', direction='out', which='both')

     pylab.title(title); pylab.xlabel(tlabel); pylab.ylabel(title_units);
                    
     if ((i > 0) & (j < (ny - 1))):
      graph_lib.clear_ticklabels(pylab.gca(), 1, 1-1)
     elif ((i == 0) & (j < (ny - 1))):
      graph_lib.clear_ticklabels(pylab.gca(), 1, 0)
     elif ((i > 0) & (j == (ny - 1))):
      graph_lib.clear_ticklabels(pylab.gca(), 0, 1-1)

   counter += 1
   ispanel = 0


def batch_sec_noaa_library():

 dpath = '/home/rilma/database/sec-noaa/lists/'
 trange = time_util.dtlim(2011, 8, [4, 7], [0, 23], [0, 59], [0, 59])

 satname = 'ace'; inst = ['mag', 'swepam'];
 data = read_all_lists(trange, dpath=dpath, inst=inst, satname=satname)

# ipar = pylab.array([[3], [12]])
 ipar = None
 tmajtick = 4*3; ntmintick = 4*3

 curr_trange = trange
 time_label = time_util.time_label(curr_trange, strutc='UT')
 tlabel = time_label.get('tlabel');

 plot_lists(data, ipar=ipar, ntmintick=ntmintick, trange=trange, \
  tmajtick=tmajtick, tlabel=tlabel)

 fnameprefix = 'sec-noaa_ace-'
 figformat = 'png'
 figsave = 1
 figshow = 0
 figdpi = 300
 gpath = '/home/rilma/tmp/results/ace/'

 figname = gpath + fnameprefix + time_label.get('tfname') + '.' + figformat
 svfig = {'save' : figsave, 'show' : figshow, 'dpi' : figdpi, \
   'figname' : figname, 'format' : figformat}
 if (svfig.get('save') > 0): graph_lib.save_figure(svfig)
 if (svfig.get('show') > 0): pylab.show()


if __name__ == '__main__':
 batch_sec_noaa_library()

