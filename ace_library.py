
import time_utilities as timelib
#import numpy
import numpy as np
import scipy.io as scio
import matplotlib.dates as matdates
import os
import file_library as filelib
import graphical_library as graphlib
import matplotlib.pyplot as plt
from matplotlib import rc

plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'

# Tuning the Figure setup (Tex)
#
rc('font',**{'family':'serif', 'serif':['New Century Schoolbook']})
rc('text', usetex=True)
#   

def get_bartels(trange):

 dimension = 5000
 start = timelib.totime(np.atleast_1d(1832), np.atleast_1d(2), \
  np.atleast_1d(8), np.atleast_1d(0), np.atleast_1d(0), np.atleast_1d(0))

 return timelib.bartels_calendar(dimension, trange, start)

def val_ranges(trange):

 imf_Bx_lim = [-25, 20]; imf_By_lim = [-30, 50]; imf_Bz_lim = [-60, 60];
 sw_lim = [200, 1000]; ief_Ey_lim = [-50, 50];
 time_shift = 0. / 60.

 if trange[0] >= timelib.totime(1991,8,20) and trange[1] < timelib.totime(1991,8,30):
  imf_Bx_lim = [-10, 10]; imf_By_lim = [-10, 10]; imf_Bz_lim = [-10, 10];

 if trange[0] >= timelib.totime(1991,8,30) and trange[1] < timelib.totime(1991,9,5):
  imf_Bx_lim = [-15, 15]; imf_By_lim = [-20, 20]; imf_Bz_lim = [-20, 20];

 if trange[0] >= timelib.totime(1993,9,1) and trange[1] <= timelib.totime(1993,9,15):
  sw_lim = [300., 450.];
  imf_Bx_lim = [-5, 5]; imf_By_lim = [-5, 5]; imf_Bz_lim = [-5, 5];

 if trange[0] >= timelib.totime(1993,9,8) and trange[1] <= timelib.totime(1993,9,9):
  sw_lim = [340., 400.]; ief_Ey_lim = [-2, 2];
  #imf_Bx_lim = [-5, 5]; imf_By_lim = [-5, 5]; imf_Bz_lim = [-5, 5];

 if trange[0] >= timelib.totime(1993,9,9) and trange[1] <= timelib.totime(1993,9,10):
  sw_lim = [320., 400.]; ief_Ey_lim = [-2, 2];

 if trange[0] >= timelib.totime(1993,9,10) and trange[1] <= timelib.totime(1993,9,11):
  sw_lim = [290., 360.]; ief_Ey_lim = [-2, 2];

 if trange[0] >= timelib.totime(2004,11,1) and trange[1] <= timelib.totime(2004,11,15):
  sw_lim = [615., 875.]; ief_Ey_lim = [-35., 35.]; imf_Bz_lim = [-45., 45.];
  time_shift = 45. / 60.

 if trange[0] >= timelib.totime(2005,3,20) and trange[1] <= timelib.totime(2005,3,25):
  imf_Bz_lim = [-6, 6]; sw_lim = [275, 375]; ief_Ey_lim = [-2., 2.];
  
 return{'SW' : sw_lim, 'E_y' : ief_Ey_lim, \
  'B_x' : imf_By_lim, 'B_y' : imf_By_lim, 'B_z' : imf_Bz_lim, \
  'time_shift' : time_shift}


def plot_ace_magswepam(data, ipar=None, ntmintick=None, \
    trange=None, tmajtick=None, tlabel=None):

 nrows, ncols = ipar.shape
 nx = ncols; ny = nrows; fdx = 0.95; fdy = 0.85;
 x0, dx, yf, dy = graphlib.axis_coord(nx, ny)

 tinfo = data.get('time');

 if trange == None: trange = np.array([tinfo[0], np.max(tinfo)])

 indt = np.where((tinfo >= trange[0]) & (tinfo <= trange[1]))
 indt = indt[0] if len(indt[0]) > 0 else np.array(range(len(tinfo)))
 
 dtmin_obj = np.min(trange); dtmax_obj = np.max(trange);
    
 xmajor_tbin = timelib.dtbin([np.floor(dtmin_obj), np.ceil(dtmax_obj)], tmajtick*60*60)
 ind  = np.where((xmajor_tbin >= dtmin_obj) & (xmajor_tbin <= dtmax_obj))
 xmajor_ticks = xmajor_tbin[ind]

 # Graphical routines
 figid = plt.figure(num=3, figsize=(8 * ncols, 4 * nrows))
 figid.clf()

 counter = 0; ispanel = 0

 imf_Bx_lim = [-25, 20]; imf_By_lim = [-30, 50]; imf_Bz_lim = [-60, 60];
 sw_x_lim = [-900., -600]; sw_y_lim = [-150., 75]; sw_z_lim = [-150., 50];
 vr = val_ranges(trange)

 for j in range(ny):
                
  for i in range(nx):
                    
   if (counter <= (nx * ny) - 1):
                        
    if ipar[j, i] == 1:
     yvalues = data['addpar']['Vsw'];
     vmin = vr['SW'][0]; vmax = vr['SW'][1]; 
     title = r'$V_{sw}$ (GSM)'; title_units = 'km/s'; ispanel = 1  
    if ipar[j, i] == 2:
     yvalues = data['addpar']['E_gsm_y'];
     vmin = vr['E_y'][0]; vmax = vr['E_y'][1]; 
     title = r'IEF, $E_y$ (GSM)'; title_units = 'mV/m'; ispanel = 1  
    if ipar[j, i] == 3:
     yvalues = data['addpar']['E_gsm_y'] / 10.;
     vmin = -3.5; vmax = 3.5; 
     title = r'IEF, $E_y/10$ (GSM)'; title_units = 'mV/m'; ispanel = 1  
    if ipar[j, i] == 11:
     yvalues = data['MAGSWE_data_64sec'].V_gsm_x[0]
     vmin = sw_x_lim[0]; vmax = sw_x_lim[1];
     title = r'$V_{sw}$ ($x$-Component, GSM)'; title_units = 'km/s'; ispanel = 1  
    if ipar[j, i] == 12:
     yvalues = data['MAGSWE_data_64sec'].V_gsm_y[0]
     vmin = sw_y_lim[0]; vmax = sw_y_lim[1];
     title = r'$Vy_{sw}$ (GSM)'; title_units = 'km/s'; ispanel = 1  
    if ipar[j, i] == 13:
     yvalues = data['MAGSWE_data_64sec'].V_gsm_z[0]
     vmin = sw_z_lim[0]; vmax = sw_z_lim[1];
     title = r'$Vz_{sw}$ (GSM)'; title_units = 'km/s'; ispanel = 1  
    if ipar[j, i] == 51:
     yvalues = data['MAGSWE_data_64sec'].B_gsm_x[0];
     vmin = imf_Bx_lim[0]; vmax = imf_Bx_lim[1]; 
     title = r'IMF, $B_x$ (GSM)'; title_units = 'nT'; ispanel = 1
    if ipar[j, i] == 52:
     yvalues = data['MAGSWE_data_64sec'].B_gsm_y[0];     
     vmin = imf_By_lim[0]; vmax = imf_By_lim[1]; 
     title = r'IMF, $B_y$ (GSM)'; title_units = 'nT'; ispanel = 1
    if ipar[j, i] == 53:
     yvalues = data['MAGSWE_data_64sec'].B_gsm_z[0];
     vmin = vr['B_z'][0]; vmax = vr['B_z'][1]; 
     title = r'IMF, $B_z$ (GSM)'; title_units = 'nT'; ispanel = 1
            
    if ispanel == 1:

     # in normal coordinates
     pn = plt.axes([x0 + i * dx, yf - (j + 1) * dy, fdx * dx, fdy * dy])
     
     xvalues = np.squeeze(data['time']) + vr['time_shift'] / 24.

     # Adding an additional "time-shift" to Nov 2004 data
     #
     if trange[0] >= timelib.totime(2004,11,6) and trange[1] <= timelib.totime(2004,11,14):
      ind_0 = np.where(xvalues >= timelib.totime(2004,11,9,23,0,0))
      if len(ind_0[0]) > 0: 
       #xvalues[ind_0[0]] += ((20. / 60.) / 24.)
       ind_1 = np.where(xvalues < timelib.totime(2004,11,9,23,0,0))
       if len(ind_1[0]) > 0:
        xvalues = np.append(np.append(xvalues[ind_1[0]], \
         np.mean([xvalues[ind_1[0][-1]], xvalues[ind_0[0][0]]])), \
         xvalues[ind_0[0]] + ((20. / 60.) / 24.))
        yvalues = np.append(np.append(yvalues[ind_1[0]], np.nan), yvalues[ind_0[0]])
     #

     plt.plot(xvalues, yvalues, color='k', linestyle='-', marker='None', markersize=3.)
     plt.plot(trange, np.zeros(2), color='k', linestyle='--')

     # Defining "box" limits
     plt.xlim(trange[0], trange[1])
     plt.ylim(vmin, vmax)

     pn.axes.set_xticks(xmajor_ticks)
     pn.axes.minorticks_on()

     # "Add minorticks" here ...        
     pn.xaxis.set_major_formatter(matdates.DateFormatter("%H:%M"))        
     pn.xaxis.set_minor_locator(plt.MultipleLocator(base=tmajtick / 24. / ntmintick))

     plt.title(title); plt.xlabel(tlabel); plt.ylabel(title_units);
                    
     if ((i > 0) & (j < (ny - 1))):
      graphlib.clear_ticklabels(plt.gca(), 1, 1)
     elif ((i == 0) & (j < (ny - 1))):
      graphlib.clear_ticklabels(plt.gca(), 1, 0)
     elif ((i > 0) & (j == (ny - 1))):
      graphlib.clear_ticklabels(plt.gca(), 0, 1)
   
   counter = counter + 1
   ispanel = 0


def read_ace_magswepam(filename):

 import datetime

 hdf_data = scio.loadmat(filename, struct_as_record=False)
        
 #time = hdf_data.get('time')[0]
 magswe_data = hdf_data.get('MAGSWE_data_64sec')[0][0]
 l2_rec = hdf_data.get('level2_record_doy')[0]

 #print hdf_data.get('MAGSWE_data_64sec')

 year = magswe_data.year[0]; day = magswe_data.day[0];
 hr = magswe_data.hr[0]; minute = magswe_data.min[0]; #second = magswe_data.sec[0]
 second = np.floor(magswe_data.sec[0]); microsec = 1e6*(magswe_data.sec[0] - second)

 #
 delta = datetime.timedelta(seconds=64)
 
 start_month, start_dom = timelib.calc_mdom(day[0], year[0])
 dstart = matdates.num2date(np.atleast_1d(timelib.totime(year[0], start_month, \
  start_dom, hr[0], minute[0], second[0], microsec[0]))[0])

 end_month, end_dom = timelib.calc_mdom(day[-1], year[-1])
 dend = matdates.num2date(np.atleast_1d(timelib.totime(year[-1], end_month, \
  end_dom, hr[-1], minute[-1], second[-1], microsec[-1]))[0] + \
  2*64. / (3600. * 24.))

 time = matdates.drange(dstart, dend, delta)
 #

 #for i in range(4):
 # print year[i], day[i], hr[i], minute[i], second[i], microsec[i]
 # print matdates.num2date(dtbins[i])

 #
 # Derived parameters
 #
 # Vsw (GSM), in unit of km/s 
 Vsw = np.sqrt(magswe_data.V_gsm_x[0]**2 + magswe_data.V_gsm_y[0]**2 + \
  magswe_data.V_gsm_z[0]**2)
 #                
 # IMF, Bz (GSM), in units of nT        
 imf_Bz = magswe_data.B_gsm_z[0]
 #                
 # IEF, Ey (GSM), in units of mV/m        
 ief_Ey = -Vsw * imf_Bz * 1e-9 * 1e6
 #
 addpar = {'Vsw' : Vsw, 'E_gsm_y' : ief_Ey}
 #

 return({'time' : time, 'MAGSWE_data_64sec' : magswe_data, 'addpar' : addpar, \
  'level2_record_doy' : l2_rec})


def main_proc_ace_magswepam(dpath=None, figformat=None, figsave=None, figshow=None, \
  gpath=None, ipar=None, ppath=None, trange=None, tmajtick=None, ntmintick=None):

 # INPUTS
 ###
 # Directory of JULIA ESF data (NCAR format)
 if dpath is None: dpath = '/home/rilma/database/ace/bartels/'
 #
 # Folder which plots are saved
 if gpath is None: gpath = '/home/rilma/tmp/results/ace/'
 filelib.create_directory(gpath)

 # Folder which reduced files are saved
 if ppath is None: ppath = '/home/rilma/tmp/results/ace/'
 filelib.create_directory(ppath)

 # Specifies the time period to analyze, in LT
 #
 if trange is None: trange = timelib.dtlim(2006, 4, 2, [0, 23], [0, 59], [0, 59])

 tzone = -5. / 24.

 if tmajtick is None: tmajtick = 3
 if ntmintick is None: ntmintick = 3

 if figformat is None: figformat = 'png'
 if figsave is None: figsave = 1
 if figshow is None: figshow = 0
 figdpi = 300

 #
 ###
 # End of INPUTS
 ###

 time_label = timelib.time_label(trange, strutc='UT')
 tlabel = time_label.get('tlabel')
 
 bc = get_bartels(trange)
 brotation = bc.get('rotation')

 nbr = brotation[1] - brotation[0] + 1
 for i in range(nbr):

  # Reading ACE data with 64 seconds of time resolution
  curr_br = brotation[0] + i        
  curr_fname = dpath + 'magswe_data_64sec_%04i.mat' % curr_br

  if os.path.isfile(curr_fname):
   print( "Reading: " + curr_fname )
   magswepam = read_ace_magswepam(curr_fname)
   
   # Time-serie plot
   figname = gpath + 'magswe-plot-64sec-' + time_label.get('tfname') + '.' + figformat
   svfig = {'save' : figsave, 'show' : figshow, 'dpi' : figdpi, \
    'figname' : figname, 'format' : figformat}
   plot_ace_magswepam(magswepam, ipar=ipar, ntmintick=ntmintick, \
    trange=trange, tmajtick=tmajtick, tlabel=tlabel)
   if (svfig.get('save') > 0): graphlib.save_figure(svfig)
   if (svfig.get('show') > 0): plt.show()   
   #dummy = magswepam['time']
   #print matdates.num2date(dummy[0])


def batch_proc_ace_magswepam():

 dpath = '/home/rilma/database/satellites/ace/bartels/'
 gpath = None; ppath = None;
 
 figformat = 'png'
 figsave = 1
 figshow = 1
 avg = 0

 tmajtick = 2; ntmintick = 4;

 if 1 == 1:

  #ipar = np.array([[51], [52], [53]])

  #trange = timelib.dtlim(2005, 3, 23, [0, 23], [0, 59], [0, 59])
  #ipar = np.array([[1], [53], [2]])

  trange = timelib.dtlim(2004, 11, [9, 10], [15, 5], [0, 59], [0, 59])
  #ipar = np.array([[2], [53], [1]])
  #ipar = np.array([[11], [12], [13]])
  ipar = np.array([[3], [53], [1]])


 main_proc_ace_magswepam(dpath=dpath, figformat=figformat, figsave=figsave, \
  figshow=figshow, gpath=gpath, ipar=ipar, ntmintick=ntmintick, ppath=ppath, \
  tmajtick=tmajtick, trange=trange)


def read_ace_data(drange):
    
    dpath = '/home/rilma/database/satellites/ace/bartels/'
    
    dimension = 5000
    start = timelib.totime(numpy.atleast_1d(1832), numpy.atleast_1d(2), numpy.atleast_1d(8), \
                           numpy.atleast_1d(0), numpy.atleast_1d(0), numpy.atleast_1d(0))
    
    bc = timelib.bartels_calendar(dimension, drange, start)
    brotation = bc.get('rotation')
    nbr = brotation[1] - brotation[0] + 1
    for i in range(nbr):
        
        # Reading ACE data with 64 seconds of time resolution
        curr_br = brotation[0] + i        
        file = dpath + 'magswe_data_64sec_%04i.mat' % curr_br
        print( "Reading: " + file )
        hdf_data = scio.loadmat(file)
        
        time = hdf_data.get('time')[0]
        magswe_data = hdf_data.get('MAGSWE_data_64sec')[0][0]
        l2_rec = hdf_data.get('level2_record_doy')[0]
        
        year = magswe_data.year[0]; day = magswe_data.day[0];
        hr = magswe_data.hr[0]; min = magswe_data.min[0];
        sec = magswe_data.sec[0]
        
        #ace_dtime = []
        #for j in range(len(day)):
        #    ace_dtime.append(timelib.totime2(year[j], day[j], hr[j], min[j], sec[j]))

        #print day, year
        #print day.shape, year.shape
        #month, dom = timelib.calc_mdom(day, year)

        #ace_dtime = timelib.totime(year, month, dom, hr, min, sec)
        print( time )

        print( 'pause' )

    
    return(0)
    
def main():
    
    #drange = timelib.dtlim(2004, [11, 12], [9, 13], [0, 23], [0, 59], [0, 59])
    #drange = timelib.dtlim(2005, 3, 23, [0, 23], [0, 59], [0, 59])
    drange = timelib.dtlim(2004, 11, [7, 13], [0, 23], [0, 59], [0, 59])
    ace_data = read_ace_data(drange)
    
    print( 'pause!' )
    
if __name__ == '__main__':
    #main()
    batch_proc_ace_magswepam()

