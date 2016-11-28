
import ace_library as acelib
import datetime
import file_library as file_lib
import goes_library as goeslib
import graphical_library as graph_lib
import iaga_library as iagalib
import magnetometer_library as maglib
import matplotlib.collections as collections
import os
import pylab
import sec_noaa_library as secnoaalib
import statistic_library as statlib
import time_utilities as time_util

pylab.rcParams['xtick.direction'] = 'out'
pylab.rcParams['ytick.direction'] = 'out'

# Tuning the Figure setup (Tex)
#
font = {'family' : 'serif', 'serif' : 'New Century Schoolbook', 'size' : 24}
pylab.rc('font', **font)
pylab.rc('text', usetex=True)


def plot_pef_1d(ipar=None, ntmintick=None, trange=None, tmajtick=None, \
 tlabel=None):

 # Turning to date/time object
 trange_obj = pylab.num2date(trange)

 is_ef = (ipar == 1) or (ipar == 7) or (ipar == 17) or (ipar == 18) or (ipar == 20) or (ipar == 28) 

 ace_ind = pylab.where(is_ef)

 ace_type = 2

 if len(ace_ind[0]) > 0:

  if ace_type == 1: # data from SEC NOAA

   dpath = '/home/rilma/database/sec-noaa/lists/'  
   satname = 'ace'; inst = ['mag', 'swepam'];
   ace_data = secnoaalib.read_all_lists(trange, dpath=dpath, inst=inst, \
    satname=satname)
   #
   # IEF, Ey, in units of mV/m        
   Vsw = ace_data['swepam']['swbs']; imf_Bz = ace_data['mag']['bz'];
   ief_Ey = -Vsw * imf_Bz * 1e-9 * 1e6
   ace_time = ace_data['mag']['time']

  if ace_type == 2: # data from Caltech

   # merged MAGSWEPAM .... ACE Satellite
   #
   dpath = '/home/rilma/database/satellites/ace/bartels/'
   bc = acelib.get_bartels(trange); brotation = bc.get('rotation')
   nbr = brotation[1] - brotation[0] + 1
   for i in range(nbr):
    # Reading ACE data with 64 seconds of time resolution
    curr_br = brotation[0] + i        
    curr_fname = dpath + 'magswe_data_64sec_%04i.mat' % curr_br
    if os.path.isfile(curr_fname):
     print( "Reading: " + curr_fname )
     magswepam = acelib.read_ace_magswepam(curr_fname)
    #
   Vsw = magswepam['addpar']['Vsw']
   imf_Bz = magswepam['MAGSWE_data_64sec'].B_gsm_z[0]
   ief_Ey = magswepam['addpar']['E_gsm_y']
   ace_time = magswepam['time']

  ief_Ey = ief_Ey * 1e-1


 jro_ind = pylab.where(is_ef)
 if len(jro_ind[0]) > 0:
  dpath = '/home/rilma/tmp/results/isr-perp/drifts-ncar-data/'
  #matdata['time'] = matdata.get('time') + tzone	# from UT to LT  
  jro_data = read_all_avg1d_drifts(trange, dpath=dpath)


 ###############################################################################
 #
 # Inputs for the JIC magnetic data filter
 #
 tstep = 5 * 60		# in seconds
 mtype = 1		# Simple mean
 #

 # Time offset to be added after data is binned
 toffset = 0.5 * tstep / (24 * 3600.0)
 #

 is_jic_mag = ((ipar == 2) | (ipar == 3) | (ipar == 6) | \
  (ipar == 7) | (ipar == 17) | (ipar == 18) | (ipar == 20) | (ipar == 28))

 jromag_ind = pylab.where(is_jic_mag)
 if len(jromag_ind[0]) > 0:

  path_mag_jro = '/home/rilma/database/jro/magnetometer/'
  # Loading pre-defined configuration values
  jro_ds = maglib.jromag_setup(trange)

  # Jicamarca magnetometer
  #
  jic_magdata = maglib.load_jromag_ascii(trange + pylab.array([0, 1.]), \
   path_mag_jro, 'jic')
  # Compute dB
  jic_magdata['dB'] = jic_magdata['H'] - jro_ds['jic']['hcomp_zero']
  # 1st derivative to H-comp (dB/dt), in nT/min
  jic_magdata['dB_by_dt'] = pylab.append(pylab.diff(jic_magdata['dB']), \
   pylab.nan) / 1.
  #

  # Smoothing the H-Comp
  raw_jic_magdata = {'time' : jic_magdata['time'], 'values' : jic_magdata['H']}
  bin_jic_magdata = statlib.ts_binned(raw_jic_magdata, tstep, mtype)
  jic_magdata['binned_time'] = bin_jic_magdata['time'] + toffset
  # 1st derivative of smoothed H-comp (in nT/min)
  jic_magdata['binned_dB_by_dt'] = pylab.append(pylab.diff( \
   bin_jic_magdata['values']), pylab.nan) / (tstep / 60.0)


  # Piura magnetometer data
  #
  piu_magdata = maglib.load_jromag_ascii(trange + pylab.array([0, 1.]), \
   path_mag_jro, 'piu')
  # Compute dB
  piu_magdata['dB'] = piu_magdata['H'] - jro_ds['piu']['hcomp_zero']
  # 1st derivative to H-comp (dB/dt)
  piu_magdata['dB_by_dt'] = pylab.append(pylab.diff(piu_magdata['dB']), \
   pylab.nan)
  #


  # Reading magnetometer data (IAGA ascii format)
  #
  dobj0 = datetime.date(trange_obj[0].year, trange_obj[0].month, \
   trange_obj[0].day)
  dobjf = datetime.date(trange_obj[1].year, trange_obj[1].month, \
   trange_obj[1].day)
  dates = []; dates.append(dobj0); dates.append(dobjf);

  path_mag_iaga = '/home/rilma/database/intermagnet/definitive/'

  # Tirunelveli (India)
  #
  tir_magdata = iagalib.load_iaga_ascii(dates, path_mag_iaga, 'tir', fntype=2)
  tir_magdata['time'] = pylab.date2num(tir_magdata['time'])
#    
#    tir_h = numpy.array(tir_magdata.get('H'))
#    # Median Filter to H-component
#    tir_h = filtlib.median(tir_h, median_width)
#
  # Default setup values (TIR)
  tir_dsetup = iagalib.iaga_dsetup(trange, 'tir')
#
  # Compute dB
  tir_magdata['dB'] = pylab.array(tir_magdata['H']) - tir_dsetup['hcomp_zero']
  #tir_dB = tir_h - tir_dsetup.get('hcomp_zero')    

  # 1st derivative to H-comp
  tir_magdata['dB_by_dt'] = pylab.append(pylab.diff(tir_magdata['dB']), pylab.nan) / 1.
  #tir_diff_H = numpy.append(numpy.diff(tir_dB), numpy.nan)    

  # Smoothing the H-Comp
  raw_tir_magdata = {'time' : tir_magdata['time'], \
   'values' : pylab.array(tir_magdata['H'])}
  bin_tir_magdata = statlib.ts_binned(raw_tir_magdata, tstep, mtype)
  tir_magdata['binned_time'] = bin_tir_magdata['time'] + toffset
  # 1st derivative to the smoothed H-comp
  tir_magdata['binned_dB_by_dt'] = pylab.append(pylab.diff( \
   bin_tir_magdata['values']), pylab.nan) / (tstep / 60.0)
#
#    # Compute dB (smooth)
#    tir_dB_sm = bin_tir_magdata.get('values') - tir_dsetup.get('hcomp_zero')
#
#    # Preparing E-Field data to be used in cross-correlation and additional processing    
#    tir_diffBdata = {'xvalues':(bin_tir_magdata.get('time') + toffset), 'yvalues':smooth_tir_diff_H}
#    tir_diffBdata = nov04lib.seldata(tir_diffBdata, dtrange)        
#
#

  # Barrow (Alaska)
  #
  brw_magdata = iagalib.load_iaga_ascii(dates, path_mag_iaga, 'brw', fntype=1)
  brw_magdata['time'] = pylab.date2num(brw_magdata['time']) 

  # Default setup values (BRW)
  brw_dsetup = iagalib.iaga_dsetup(trange, 'brw')

  # Compute dB
  brw_magdata['dB'] = pylab.array(brw_magdata['H']) - brw_dsetup['hcomp_zero']

  # 1st derivative to H-comp
  brw_magdata['dB_by_dt'] = pylab.append(pylab.diff(brw_magdata['dB']), pylab.nan) / 1.

  # Smoothing the H-Comp
  raw_brw_magdata = {'time' : brw_magdata['time'], \
   'values' : pylab.array(brw_magdata['H'])}
  bin_brw_magdata = statlib.ts_binned(raw_brw_magdata, tstep, mtype)
  brw_magdata['binned_time'] = bin_brw_magdata['time'] + toffset

  # 1st derivative to the smoothed H-comp
  brw_magdata['binned_dB_by_dt'] = pylab.append(pylab.diff( \
   bin_brw_magdata['values']), pylab.nan) / (tstep / 60.0)


  # College (Alaska)
  #
  cmo_magdata = iagalib.load_iaga_ascii(dates, path_mag_iaga, 'cmo', fntype=1)
  cmo_magdata['time'] = pylab.date2num(cmo_magdata['time']) 

  # Default setup values (CMO)
  cmo_dsetup = iagalib.iaga_dsetup(trange, 'cmo')

  # Compute dB
  cmo_magdata['dB'] = pylab.array(cmo_magdata['H']) - cmo_dsetup['hcomp_zero']

  # 1st derivative to H-comp
  cmo_magdata['dB_by_dt'] = pylab.append(pylab.diff(cmo_magdata['dB']), pylab.nan) / 1.

  # Smoothing the H-Comp
  raw_cmo_magdata = {'time' : cmo_magdata['time'], \
   'values' : pylab.array(cmo_magdata['H'])}
  bin_cmo_magdata = statlib.ts_binned(raw_cmo_magdata, tstep, mtype)
  cmo_magdata['binned_time'] = bin_cmo_magdata['time'] + toffset

  # 1st derivative to the smoothed H-comp
  cmo_magdata['binned_dB_by_dt'] = pylab.append(pylab.diff( \
   bin_cmo_magdata['values']), pylab.nan) / (tstep / 60.0)


  #

 is_ef = (True if len(pylab.where(is_ef)[0]) > 0 else False)
 is_jic_mag = (True if len(pylab.where(is_jic_mag)[0]) > 0 else False)
 if is_ef & is_jic_mag:
  dtrange = trange
  dtrange = time_util.dtlim(2004, 11, 10, [7, 8], [0, 29], [0, 59])
  var1 = {'time' : jro_data['time'], 'value' : jro_data['Ex']}
  var2 = {'time' : jic_magdata['binned_time'], 'value' : jic_magdata['binned_dB_by_dt']}
  dummy = statlib.corrcoef(var1, var2, dtrange, detrend=None)


 #is_goes = ((ipar == 4) | (ipar == 5) | (ipar == 7) | (ipar == 40))
 is_goes = ((ipar == 4) | (ipar == 5) | (ipar == 40))
 goes_ind = pylab.where(is_goes)
 if len(goes_ind[0]) > 0:
  path_goes = '/home/rilma/database/spidr/goes/'
  # Set the GOES data folder
  path_goes = path_goes + 'mag-%s' % trange_obj[0].strftime('%Y%m%d') + os.sep
  # Read GOES-12 data (75 W)
  goes12 = goeslib.load_spidr_all1sd(path_goes, trange[0], 'goes12', 'hp')
  # Loading pre-defined configuration values
  gms12 = goeslib.goes_mag_setup(goes12['station'], trange)
#  # Median Filter to Hp-component
#  goes12_hp = filtlib.median(goes12['hp'], median_width)
#  # 1st derivative to H-comp
#  goes12_diff_hp = numpy.append(numpy.diff(goes12_hp), numpy.nan)    
#    
  # Smoothing the hp-Comp
  raw_goes12_hp = {'time' : goes12['time'], 'values' : goes12['hp']}
  bin_goes12_hp = statlib.ts_binned(raw_goes12_hp, tstep, mtype)
  goes12['binned_hp'] = bin_goes12_hp['values']
  goes12['binned_time'] = bin_goes12_hp['time'] + toffset
  # 1st derivative to the smoothed hp-comp
  goes12['binned_dB_by_dt'] = \
   pylab.append(pylab.diff(bin_goes12_hp['values']), pylab.nan) / (tstep / 60.0)    

  if len(pylab.where(ipar == 40)[0]) > 0:

   import proc_goes_sat as pgoess

   dpath = '/home/rilma/database/satellites/goes/'
   isat = 12

   date = [trange_obj[0].year, trange_obj[0].month, trange_obj[0].day]
   goes = pgoess.read_goes_data(date=date, dpath=dpath, isat=isat)
   goes12['sc_pos_gse'] = goes['sc_pos_gse']

#   for i in range(len(goes['time'])):
#    print pylab.num2date(goes12['time'][i]), pylab.num2date(goes['time'][i])
#   return

##


 nrows, ncols = ipar.shape
 nx = ncols; ny = nrows; fdx = 0.85; fdy = 0.85;
 x0, dx, yf, dy = graph_lib.axis_coord(nx, ny)

 if len(pylab.where(ipar == 40)[0]) > 0:
  #fdy = 0.65; 
  dy -= 0.1

 #if trange == None: trange = pylab.array([tinfo[0], max(tinfo)])

 #indt = pylab.where((tinfo >= trange[0]) & (tinfo <= trange[1]))
 #indt = indt[0] if len(indt[0]) > 0 else pylab.array(range(len(tinfo)))
 
 dtmin_obj = min(trange); dtmax_obj = max(trange);
    
 xmajor_tbin = time_util.dtbin([pylab.floor(dtmin_obj), pylab.ceil(dtmax_obj)], tmajtick*60*60)
 ind = pylab.where((xmajor_tbin >= dtmin_obj) & (xmajor_tbin <= dtmax_obj))
 xmajor_ticks = xmajor_tbin[ind]

 # Graphical routines
 figid = pylab.figure(num=42, figsize=(10 * ncols, 6 * nrows))
 figid.clf()

 counter = 0; ispanel = 0

 for j in range(ny):
                
  for i in range(nx):
                    
   if (counter <= (nx * ny) - 1):

    if ipar[j, i] == 1: # Electric Field
     xvalues = ace_time; yvalues = ief_Ey;
     vmin = valrange(trange)['jic_ex'][0]; vmax = valrange(trange)['jic_ex'][1];
     title = 'Electric Field'; label = '$IEF_y (ACE)$';
     color = 'r'; ymintick = 0.25
     title_units = 'mV/m'; linestyle = '-'; ispanel = 1

    if (ipar[j, i] == 2) or (ipar[j, i] == 6) or (ipar[j, i] == 20):
     xvalues = jic_magdata['time']; yvalues = jic_magdata['dB'];
     vmin = valrange(trange)['jic_hcomp'][0]; vmax = valrange(trange)['jic_hcomp'][1];
     title = r'$H_{comp}$'; label = 'JIC';
     color = 'k'; ymintick = 10.
     title_units = 'nT'; linestyle = '-'; ispanel = 1;

    if (ipar[j, i] == 3) or (ipar[j, i] == 7):
     xvalues = jic_magdata['binned_time']; #yvalues = jic_magdata['dB_by_dt'];
     yvalues = jic_magdata['binned_dB_by_dt']
     #vmin = valrange(trange)['jic_dB_by_dt'][0]; vmax = valrange(trange)['jic_dB_by_dt'][1];
     vmin = valrange(trange)['tir_dB_by_dt'][0]; vmax = valrange(trange)['tir_dB_by_dt'][1];
     #title = r'$\frac{\partial B}{\partial t}$'; label = 'JIC';
     title = r'$\partial B/\partial t$'; label = 'JIC';
     color = 'k'; ymintick = 1.
     title_units = 'nT/min'; linestyle = '-'; ispanel = 1;

    if (ipar[j, i] == 4) or (ipar[j, i] == 40): 
     xvalues = goes12['time']; yvalues = goes12['hp'];
     #yvalues = goes12_hp;
     vmin = valrange(trange)['goes_hp'][0]; vmax = valrange(trange)['goes_hp'][1];
     #ylim = gms12.get('mag_lim'); yminorticks = gms12.get('mag_minorticks');
     title = goes12['station'] + r': Long. W 75$^\circ$';
     label = r'$h_p$'
     color = 'k'; ymintick = 5.#gms12['mag_minorticks'];
     title_units = 'nT'; linestyle = '-'; ispanel = 1;

    if ipar[j, i] == 5:
     xvalues = goes12['binned_time']; yvalues = goes12['binned_dB_by_dt'];
     vmin = valrange(trange)['goes_dB_by_dt'][0]; vmax = valrange(trange)['goes_dB_by_dt'][1];
     #ylim = gms12.get('mag_lim'); yminorticks = gms12.get('mag_minorticks');
     title = goes12['station'] + r': Long. W 75$^\circ$';
     label = r'$\partial B/\partial t, h_p$'
     color = 'k'; ymintick = 0.5#gms12['mag_minorticks'];
     title_units = 'nT/min'; linestyle = '-'; ispanel = 1;

    if ipar[j, i] == 17:
     #xvalues = brw_magdata['time']; #yvalues = jic_magdata['dB_by_dt'];
     ##yvalues = pylab.array(brw_magdata['H'])
     #yvalues = pylab.array(brw_magdata['dB'])
     xvalues = brw_magdata['binned_time']; yvalues = brw_magdata['binned_dB_by_dt'] / 20.;
     #vmin = valrange(trange)['jic_dB_by_dt'][0]; vmax = valrange(trange)['jic_dB_by_dt'][1];
     #vmin = 9150; vmax = 9300;
     #vmin = -1750; vmax = 250;
     vmin = -10; vmax = 10;
     #vmin = -200; vmax = 200;
     #title = r'$H_{comp}$'; label = 'BRW';
     #title = r'$\frac{\partial B}{\partial t}$'; label = 'JIC';
     title = r'$\partial B/\partial t$'; label = 'BRW (Scaled)';
     #color = 'k'; ymintick = 50.
     #color = 'k'; ymintick = 1.
     color = 'k'; ymintick = 10.
     #title_units = 'nT'; linestyle = '-'; ispanel = 1;
     title_units = 'nT/min'; linestyle = '-'; ispanel = 1;

    if ipar[j, i] == 18:
     xvalues = cmo_magdata['time']; yvalues = cmo_magdata['dB'];
     #yvalues = pylab.array(cmo_magdata['H'])
     vmin = -1500; vmax = 250;
     #vmin = cmo_dsetup['hcomp_range'][0]; vmax = cmo_dsetup['hcomp_range'][1]
     #title = r'$H_{comp}$'; label = 'CMO';
     title = ''; label = r'H$_{comp}$ (CMO)';
     color = 'k'; ymintick = 50.
     title_units = 'nT'; linestyle = '-'; ispanel = 1;

    if ipar[j, i] == 28:
     xvalues = cmo_magdata['binned_time']; yvalues = cmo_magdata['binned_dB_by_dt'] / 8;
     vmin = valrange(trange)['tir_dB_by_dt'][0]; vmax = valrange(trange)['tir_dB_by_dt'][1];
     title = r'$\partial B/\partial t$'; label = 'CMO / 8';
     color = 'k'; ymintick = 1.;
     title_units = 'nT/min'; linestyle = '-'; ispanel = 1;


    if ispanel == 1:

     # in normal coordinates
     pn = pylab.axes([x0 + i * dx, yf - (j + 1) * dy, fdx * dx, fdy * dy])     

     #
     # Using "plot"
     #

     indt = pylab.where((xvalues >= trange[0]) & (xvalues <= trange[1]))
     indt = indt[0] if len(indt[0]) > 0 else pylab.array(range(len(tinfo)))

     x = xvalues[indt]; y = yvalues[indt];
#     ipn = pylab.plot(x, y, clip_on=True, color=color, label=label, \
#      linestyle=linestyle)
     p0, = pylab.plot(x, y, clip_on=True, color=color, label=label, \
      linestyle=linestyle)
     if ipar[j, i] == 1: pylab.plot(jro_data['time'], jro_data['Ex'], \
      color='k', label=r'E$_X$ (JIC)', linestyle='-')
     if ipar[j, i] == 2: pylab.plot(piu_magdata['time'], piu_magdata['dB'], \
      color='k', label='PIU', linestyle='--')
     #if ipar[j, i] == 3: pylab.plot(piu_magdata['time'], piu_magdata['dB_by_dt'], \
     # color='k', label='PIU', linestyle='--')
     if ipar[j, i] == 6: pylab.plot(tir_magdata['time'], tir_magdata['dB'], \
      color='k', label='TIR', linestyle='--')
     if ipar[j, i] == 7: p1, = pylab.plot(tir_magdata['binned_time'], \
      tir_magdata['binned_dB_by_dt'], color='k', label='TIR', linestyle='--')     
     #if ipar[j, i] == 7: p2, = pylab.plot(goes12['binned_time'], \
     # goes12['binned_dB_by_dt'], color='k', \
     # label=goes12['station']+r': W 75$^\circ$', linestyle='-.')
     if (ipar[j, i] == 17) | (ipar[j, i] == 28): p1, = pylab.plot(jic_magdata['binned_time'], \
      jic_magdata['binned_dB_by_dt'], color='k', label='JIC', linestyle='--')

     ipn0 = pylab.plot(trange, pylab.zeros(2), linestyle=':', color='k')

     if (ipar[j, i] == 7) or (ipar[j, i] == 17) or (ipar[j, i] == 18) or (ipar[j, i] == 20) or (ipar[j, i] == 28) :


      ax3 = pn.twinx()
      p3, = ax3.plot(jro_data['time'], jro_data['Ex'], color='r', \
       label=r'E$_X$ (JIC)', linestyle='-')
      ax3.set_ylabel('mV/m')

#      ax3_attr = dir(pylab.axes)
#      for k in range(len(ax3_attr)):       
#       print k, ax3_attr[k]
#       if k == 388: return

      ax3.yaxis.label.set_color(p3.get_color())      
      #ax3.tick_params(axis='y', colors=p3.get_color())
      ax3.set_xlim(trange[0], trange[1]); ax3.set_ylim(-2.5, 2.5);

      if (ipar[j, i] != 18) & (ipar[j, i] != 28) :
       if (ipar[j, i] != 20):
        barh_xrange = time_util.dtlim(2004, 11, 10, [7, 8], [38, 28], [0, 59])
        y1 = pylab.zeros(len(y))
        iy1 = pylab.where((x >= barh_xrange[0]) & (x <= barh_xrange[1]))
        y1[iy1[0]] = 1.
        collection = collections.BrokenBarHCollection.span_where(x, \
         ymin=-10., ymax=10., where=y1>0, facecolor='yellow', alpha=0.5)
        pn.add_collection(collection)

     if (ipar[j, i] == 7) or (ipar[j, i] == 17) or (ipar[j, i] == 18) or (ipar[j, i] == 20) or (ipar[j, i] == 28) :
      #lines = [p0, p1, p3] if ipar[j, i] == 17 else [p0, p1, p2, p3]       
      if ipar[j, i] == 17:
       lines = [p0, p1, p3]
      elif (ipar[j, i] == 18) or (ipar[j, i] == 20) :
       lines = [p0, p3]
      else:
       #lines = [p0, p1, p2, p3]
       lines = [p0, p1, p3]
      
      pn.legend(lines, [l.get_label() for l in lines], fancybox=True, \
       loc='upper right', prop={'size':12})
      #attr = dir(pn)
      #for i in range(len(attr)):
      # print attr[i]
     else:
      if ipar[j, i] == 1: pylab.legend(fancybox=True, fontsize=16, loc='upper right')
      print( 'No legend!' )

      
      #ax2.yaxis.set_minor_locator(pylab.MultipleLocator(base=0.25))

     # Defining "box" limits
     pn.set_xlim(trange[0], trange[1])
     pn.set_ylim(vmin, vmax)

#     print vmin, vmax
     #

     pn.axes.set_xticks(xmajor_ticks);
     pn.axes.minorticks_on()

     pn0 = pylab.gca()
     pn0.xaxis.set_minor_locator(pylab.MultipleLocator(base=tmajtick / 24. / ntmintick))
     pn0.yaxis.set_minor_locator(pylab.MultipleLocator(base=ymintick))

     # "Add minorticks" here ...        
     pn.xaxis.set_major_formatter(pylab.DateFormatter("%H:%M"))        


     if ipar[j, i] == 40:

      import scipy.interpolate as scint
 
      # Adding multiple tick labels ...
      #
      xticks_loc = pn.get_xticks()

      sc_time = goes['time']
      sc_pos_x = pylab.squeeze(goes12['sc_pos_gse'][:, 0])
      sc_pos_y = pylab.squeeze(goes12['sc_pos_gse'][:, 1])
      sc_pos_z = pylab.squeeze(goes12['sc_pos_gse'][:, 2])

      sc_time = pylab.append(sc_time[0] - (sc_time[1] - sc_time[0]), sc_time)
      sc_pos_x = pylab.append(sc_pos_x[0] - (sc_pos_x[1] - sc_pos_x[0]), sc_pos_x)
      sc_pos_y = pylab.append(sc_pos_y[0] - (sc_pos_y[1] - sc_pos_y[0]), sc_pos_y)
      sc_pos_z = pylab.append(sc_pos_z[0] - (sc_pos_z[1] - sc_pos_z[0]), sc_pos_z)

      f = scint.interp1d(sc_time, sc_pos_x)
      sc_pos_x_ticks_loc = f(xticks_loc)

      f = scint.interp1d(sc_time, sc_pos_y)
      sc_pos_y_ticks_loc = f(xticks_loc)

      f = scint.interp1d(sc_time, sc_pos_z)
      sc_pos_z_ticks_loc = f(xticks_loc)
   
      xticklabels = []
      for i in range(len(xticks_loc)):       
       hour = pylab.num2date(xticks_loc[i]).hour
       minute = pylab.num2date(xticks_loc[i]).minute
       ut_time_label = '%02i:%02i' % (hour, minute)
       curr_xticks = '%s\n%6.1f\n%6.1f\n%6.1f' % (ut_time_label, \
        sc_pos_x_ticks_loc[i], sc_pos_y_ticks_loc[i], sc_pos_z_ticks_loc[i])
       xticklabels.append(curr_xticks)
      xticklabels[0] = 'UT\nX(GSE)\nY(GSE)\nZ(GSE)'
      pn0.set_xticklabels(xticklabels, fontsize=15.)

      pylab.title(title); pn.set_xlabel(''); pn.set_ylabel(title_units);

     else:
      if ipar[j, i] == 20: title = ''
#     pylab.title(title); pylab.xlabel(tlabel); pylab.ylabel(title_units);
      pylab.title(title); pn.set_xlabel(tlabel); pn.set_ylabel(title_units);


     if ((i > 0) & (j < (ny - 1))):
      graph_lib.clear_ticklabels(pylab.gca(), 1, 1-1)
     elif ((i == 0) & (j < (ny - 1))):
      graph_lib.clear_ticklabels(pylab.gca(), 1, 0)
     elif ((i > 0) & (j == (ny - 1))):
      graph_lib.clear_ticklabels(pylab.gca(), 0, 1-1)

   counter += 1
   ispanel = 0
# End of 'plot_pef_1d'


def read_all_avg1d_drifts(trange, dpath=None):

 import scipy.io as scio

 ndays = int(trange[1] - trange[0]) + 1

 for i in range(ndays):

  curr_trange = pylab.floor(trange[0]) + pylab.array([0., 0.999]) + i
  curr_date = curr_trange[0]
  dt = pylab.num2date(curr_date)
  mat_fname = 'jro%04i%02i%02i_avg1d_drifts.mat' % (dt.year, dt.month, dt.day)
  curr_fname = dpath + mat_fname

  print( 'Reading ... %s' % curr_fname )
  curr_matdata = scio.loadmat(curr_fname, struct_as_record=False)
  
  if i == 0:
   matdata = curr_matdata
  else:
   matdata = {'time' : pylab.append(matdata['time'], curr_matdata['time'], axis=0), \
    'vd' : pylab.append(matdata['vd'], curr_matdata['vd'], axis=0), \
    'errvd' : pylab.append(matdata['errvd'], curr_matdata['errvd'], axis=0), \
    'zd' : pylab.append(matdata['zd'], curr_matdata['zd'], axis=0), \
    'errzd' : pylab.append(matdata['errzd'], curr_matdata['errzd'], axis=0), \
    'nvpoints' : pylab.append(matdata['nvpoints'], curr_matdata['nvpoints'], axis=0), \
    'ntpoints' : pylab.append(matdata['ntpoints'], curr_matdata['ntpoints'], axis=0), \
    'Ex' : pylab.append(matdata['Ex'], curr_matdata['Ex'], axis=0), \
    'errEx' : pylab.append(matdata['errEx'], curr_matdata['errEx'], axis=0)}

 return(matdata)
# End of 'read_all_avg1d_drifts'


def valrange(trange):

 jic_dB_by_dt, tir_dB_by_dt = None, None

 tperiod = time_util.dtlim(2004, 11, [9, 13], [0, 23], [0, 59], [0, 59]) 
 if (trange[0] >= tperiod[0]) & (trange[1] <= tperiod[1]):
  jic_ex_range = [-4., 4.]; jic_hcomp_range = [-400., 200.];
  goes_hp = [-25., 150.]; goes_dB_by_dt = [-100., 100.]

 tperiod = time_util.dtlim(2004, 11, 10, [5, 10], [0, 59], [0, 59]) 
 if (trange[0] >= tperiod[0]) & (trange[1] <= tperiod[1]):
  jic_ex_range = [-2.5, 1.5]; jic_hcomp_range = [-275., 50.];
  jic_dB_by_dt = [-5., 7.]; tir_dB_by_dt = [-10., 10.]; goes_hp = [-10., 100.]
  goes_dB_by_dt = [-5., 10.]

 output = {'jic_ex' : jic_ex_range, 'jic_hcomp' : jic_hcomp_range, \
  'jic_dB_by_dt' : jic_dB_by_dt, 'tir_dB_by_dt' : tir_dB_by_dt, \
  'goes_hp' : goes_hp, 'goes_dB_by_dt' : goes_dB_by_dt}
 return(output)
# End of 'dsetup'


def batch_proc_pef():

## trange = time_util.dtlim(2011, 8, [5, 6], [0, 23], [0, 59], [0, 59])
# trange = time_util.dtlim(2004, 11, [9, 10], [0, 23], [0, 59], [0, 59])
# tmajtick = 6; ntmintick = 6;


## Fig. 3
# trange = time_util.dtlim(2004, 11, 10, [5, 10], [0, 59], [0, 59])    
# tmajtick = 1; ntmintick = 3;
# ipar = pylab.array([[1],[2]])

## Fig. 3b
# trange = time_util.dtlim(2004, 11, 10, [7, 8], [0, 29], [0, 59])    
# tmajtick = 1./3; ntmintick = 4;
# ipar = pylab.array([[20]])

## Fig. 4
# trange = time_util.dtlim(2004, 11, 10, [7, 8], [0, 59], [0, 59])    
# tmajtick = 1./3; ntmintick = 4;
# ipar = pylab.array([[1],[3]])

### Fig.5
## trange = time_util.dtlim(2004, 11, 10, [7, 8], [0, 59], [0, 59])     
## tmajtick = 1./3; ntmintick = 4;
## ipar = pylab.array([[5],[7]])

## another Fig.5 (I like this better)
# trange = time_util.dtlim(2004, 11, 10, [7, 8], [0, 59], [0, 59])     
# tmajtick = 1./3; ntmintick = 4;
# ipar = pylab.array([[7]])

# Fig.5b
 trange = time_util.dtlim(2004, 11, 10, [7, 8], [0, 59], [0, 59])     
 tmajtick = 1./3; ntmintick = 4;
 ipar = pylab.array([[40]])

## Fig. 6
# trange = time_util.dtlim(2004, 11, 10, [7, 11], [0, 59], [0, 59])     
# tmajtick = 1.; ntmintick = 6;
# ipar = pylab.array([[17]])

## Fig. 6 (new)
# tmajtick = 1.; ntmintick = 6;
# trange = time_util.dtlim(2004, 11, 10, [7, 8], [38, 28], [0, 59])
# tmajtick = 1./6; ntmintick = 4;
# ipar = pylab.array([[18]])

## Fig. 6 (new 02)
# tmajtick = 1.; ntmintick = 6;
# trange = time_util.dtlim(2004, 11, 10, [7, 8], [38, 28], [0, 59])
# tmajtick = 1./6; ntmintick = 4;
# ipar = pylab.array([[28]])

 curr_trange = trange
 time_label = time_util.time_label(curr_trange, strutc='UT')
 tlabel = time_label.get('tlabel');

 plot_pef_1d(ipar=ipar, ntmintick=ntmintick, trange=trange, tmajtick=tmajtick, \
  tlabel=tlabel)

 fnameprefix = 'pef-'
 #figformat = 'png'
 figformat = 'eps'
 figsave = 1
 figshow = 0
 figdpi = 300

 # Folder which plots are saved
 gpath = '/home/rilma/tmp/results/pef/'
 file_lib.create_directory(gpath)

 figname = gpath + fnameprefix + time_label.get('tfname') + '.' + figformat
 svfig = {'save' : figsave, 'show' : figshow, 'dpi' : figdpi, \
   'figname' : figname, 'format' : figformat}
 if (svfig.get('save') > 0): graph_lib.save_figure(svfig)
 if (svfig.get('show') > 0): pylab.show()
## End of 'batch_proc_pef'


if __name__ == '__main__':
 batch_proc_pef()


