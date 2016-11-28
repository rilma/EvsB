
import datetime
from matplotlib.dates import date2num, drange, DateFormatter, HourLocator, MinuteLocator, num2date
from matplotlib.ticker import MultipleLocator
#from pylab import axis, figure, savefig, show, title, xlabel, ylabel
from matplotlib.pyplot import axis, figure, savefig, show, title, xlabel, ylabel
import numpy
import os
import time_utilities

def read_iaga_ascii(filename):
    
    print( "Reading: " + filename )
    
    hlock = 0; lock = 0; counter = 0; dscounter = 0;
    
    dtcomp = []; xcomp = []; ycomp = []; zcomp = []; ftotal = []
    hcomp = []; dcomp = [];

    #f = open(filename, 'r')    
    f = open(filename, 'r', errors='ignore')
    
    for line in f:
        
        # Reading header infomation
        if ((hlock == 0) & (line[1:2] != '#')):
            fdname = line[1:24].rstrip(); fdvalue = line[24:(len(line)-2-1)].rstrip()
            if (counter == 0):
                header = { fdname : fdvalue }
            else:
                header[fdname] = fdvalue
        else:
            hlock = 1
        
        # Select the initial row number to read data columns
        if (line[0:4] == "DATE"):
            lock = 1; dscounter = counter;
        
        # Reading data columns
        if ((lock == 1) and ((dscounter - counter) < 0)):
            
            dt_str = line[0:23]
            d = dt_str.split(' ')[0].split('-')
            dint = [int(i) for i in d]
            t = dt_str.split(' ')[1].split(':')
            tint = [float(i) for i in t]
            ms = int((tint[2] - float(int(tint[2])))*1e3)
            dobj = datetime.date(dint[0], dint[1], dint[2])
            tobj = datetime.time(int(tint[0]), int(tint[1]), int(tint[2]), ms)
            dtval = datetime.datetime.combine(dobj, tobj)
            
            dtcomp.append(dtval)
            
            if (header.get('Reported') == 'XYZF'):
                xcomp.append(float(line[30:40]))
                ycomp.append(float(line[40:50]))
            
            if ((header.get('Reported') == 'HDZF') | (header.get('Reported') == 'HDZ')):
                hcomp.append(float(line[30:40]))
                dcomp.append(float(line[40:50]))
            
            zcomp.append(float(line[50:60]))
            ftotal.append(float(line[60:70]))
            
        counter = counter + 1
        
    f.close()
        
    output = {'time' : dtcomp, 'X' : xcomp,  'Y' : ycomp,  'Z' : zcomp, \
              'F' : ftotal, 'H' : hcomp, 'D' : dcomp, 'header' : header} 
    
    return(output)

def load_iaga_ascii(dates, path, station='hua', fntype=1):
 
    time = []; xcomp = []; ycomp = []; zcomp = []; ftotal = [];
    hcomp = []; dcomp = [];
 
    dtrange = date2num(dates)
    timedelta = datetime.timedelta(days=1)

    datetime_date = drange(num2date(dtrange[0]), num2date(dtrange[1] + 1.0), timedelta)

    for i in range(len(datetime_date)):
        dobj = num2date(datetime_date[i])
        if (fntype == 1):
            fname = "%s%04i%02i%02id.min" % (station, dobj.year, dobj.month, dobj.day)
        elif (fntype == 2):
            fname = ('%s.%s') % (dobj.strftime("%b%d%y").lower(),station)
        else:
            print( '?' )
        filename = path + ("%s%04i" % (station, dobj.year)) + os.sep + fname
        if (not os.access(filename, os.R_OK)):
            print( 'The file %s could not be found!' % filename )
        else:
            
            magdata = read_iaga_ascii(filename)            
            
            time = numpy.append(time, magdata.get('time'))
            
            if (magdata.get('header').get('Reported') == 'XYZF'):
                xcomp_temp = magdata.get('X'); ycomp_temp = magdata.get('Y');
                xcomp = numpy.append(xcomp, xcomp_temp); ycomp = numpy.append(ycomp, ycomp_temp);
                hcomp = numpy.append(hcomp, numpy.sqrt(numpy.array(xcomp_temp)**2 + numpy.array(ycomp_temp)**2))
            
            if ((magdata.get('header').get('Reported') == 'HDZF') | (magdata.get('header').get('Reported') == 'HDZ')):
                hcomp = numpy.append(hcomp, magdata.get('H'))
                dcomp = numpy.append(dcomp, magdata.get('D'))
            
            zcomp = numpy.append(zcomp, magdata.get('Z'))            
            ftotal = numpy.append(ftotal, magdata.get('F'))
            
#    hcomp = numpy.sqrt(xcomp**2 + ycomp**2)
        
    output = {'time':list(time), 'X':list(xcomp), 'Y':list(ycomp), 'Z':list(zcomp),\
              'F':list(ftotal), 'H':list(hcomp)} 
    
    # Replacing bad values (99999.0) with NaN
    fdname = list(output.keys())
    for i in range(len(fdname)):
        if (fdname[i] != 'time'):
            fdvalue = numpy.array(output.get(fdname[i]))
            ind = numpy.where((fdvalue == 99999.0) | (fdvalue == 99999.90))
            if (len(ind[0]) > 0):
                fdvalue[ind[0]] = numpy.nan
                output[fdname[i]] = list(fdvalue)    
    
    return(output)

def iaga_dsetup(dtrange, station):

    dt_lim = time_utilities.dtlim(2003, [10, 11], [20, 5], [0, 23], [0, 59], [0, 59])
    if ((dtrange[0] >= dt_lim[0]) & (dtrange[1] <= dt_lim[1])):
        if (station == 'frd'):            
            hcomp_range = [20550, 21550]; hcomp_step = 100; 
            hcomp_zero = 21010.40; # Midnight value on previous quiet day: Oct 23, 2003
            zcomp_range = [48000, 49000]; zcomp_step = 100            
        else:  
            hcomp_range = [3300, 4600]; hcomp_step = 100; hcomp_zero = 0.0;
            zcomp_range = [56000, 57100]; zcomp_step = 100

    dt_lim = time_utilities.dtlim(2004, 11, [6, 15], [0, 23], [0, 59], [0, 59])
    if ((dtrange[0] >= dt_lim[0]) & (dtrange[1] <= dt_lim[1])):
        if (station == 'abg'):
            hcomp_range = [37700, 38200]; hcomp_step = 100;
            #hcomp_zero = 38119.20;   # Midnight value on quietest day on Nov 2004 (6th)
            hcomp_zero = 38064.50;   # Midnight value on 1st quiet day after storm: Nov 15, 2004
            zcomp_range = [18800, 18900]; zcomp_step = 10;
        elif (station == 'asc'):            
            hcomp_range = [20700, 21300]; hcomp_step = 100;
            hcomp_zero = 21032.10;   # Midnight value on 1st quiet day after storm: Nov 15, 2004
            zcomp_range = [390, 480]; zcomp_step = 10            
        elif (station == 'cmo'):	# College, Alaska
            hcomp_range = [10000., 12500.]; hcomp_step = 10;
            hcomp_zero = 12663.		# Midnight value on 1st quiet day after storm: Nov 15, 2004
            zcomp_range = [390, 480]; zcomp_step = 10
        elif (station == 'brw'):	# Barrow, Alaska
            hcomp_range = [9150., 9300.]; hcomp_step = 10;
            hcomp_zero = 9255.0;   # Midnight value on 1st quiet day after storm: Nov 15, 2004
            zcomp_range = [390, 480]; zcomp_step = 10
        elif (station == 'gua'):
            hcomp_range = [35250, 35850]; hcomp_step = 100; 
            hcomp_zero = 35675.44;   # Midnight value on 1st quiet day after storm: Nov 15, 2004
            zcomp_range = [7800, 7875]; zcomp_step = 25            
        elif (station == 'hua'):
            hcomp_range = [25350, 26150]; hcomp_step = 100; hcomp_zero = 0.0;
            zcomp_range = [390, 480]; zcomp_step = 10
        elif (station == 'thl'):
            hcomp_range = [3300, 4600]; hcomp_step = 100; hcomp_zero = 0.0;
            zcomp_range = [56000, 57100]; zcomp_step = 100
        elif (station == 'tir'):            
            hcomp_range = [39750, 40400]; hcomp_step = 100;
            #hcomp_zero = 40242.10;   # Midnight value on 1st quiet day after storm: Nov 15, 2004
            hcomp_zero = 40196.70;   # Midnight value on 1st quiet day after storm: Nov 15, 2004
            zcomp_range = [670, 940]; zcomp_step = 50;
        else:  
            hcomp_range = [3300, 4600]; hcomp_step = 100; hcomp_zero = 0.0;
            zcomp_range = [56000, 57100]; zcomp_step = 100

    dt_lim = time_utilities.dtlim(2005, 3, [1, 31], [0, 23], [0, 59], [0, 59])
    if ((dtrange[0] >= dt_lim[0]) & (dtrange[1] <= dt_lim[1])):
        if (station == 'hua'):            
            hcomp_range = [25650, 25950]; hcomp_step = 50; hcomp_zero = 0.0;
            zcomp_range = [400, 450]; zcomp_step = 10;
        else:  
            hcomp_range = [3300, 4600]; hcomp_step = 100; hcomp_zero = 0.0;
            zcomp_range = [56000, 57100]; zcomp_step = 100
    
    dt_lim = time_utilities.dtlim(2008, 6, [14, 21], [0, 23], [0, 59], [0, 59])
    if ((dtrange[0] >= dt_lim[0]) & (dtrange[1] <= dt_lim[1])):
        if (station == 'hua'):            
            hcomp_range = [25400, 25700]; hcomp_step = 50; hcomp_zero = 0.0;
            zcomp_range = [250, 350]; zcomp_step = 25;
        else:  
            hcomp_range = [3300, 4600]; hcomp_step = 100; hcomp_zero = 0.0;
            zcomp_range = [56000, 57100]; zcomp_step = 100
    
    return({'hcomp_range':hcomp_range, 'hcomp_step':hcomp_step, 'hcomp_zero':hcomp_zero, \
            'zcomp_range':zcomp_range, 'zcomp_step':zcomp_step})
    
def plot_iaga(data, gsetup):
 
    dtime = data.get('time')
    xvalues = date2num(dtime)
    
    iagas = iaga_dsetup(xvalues[[0, len(xvalues) - 1]], gsetup.get('station'))
    
    xlim = date2num(gsetup.get('xlim'))
#    ylim = gsetup.get('ylim')

    xsmajorl = HourLocator(None, gsetup.get('xmajorticks_step'))
    xsminorl = MinuteLocator(None, gsetup.get('xminorticks_step'))
    
#    ysmajorl = MultipleLocator(gsetup.get('ymajorticks_step'))
    ysminorl = MultipleLocator(gsetup.get('yminorticks_step'))
    
    par = [5, 3]; npanels = len(par);
    
    fig = figure(num=1, figsize=(16, 6*npanels))
    
    for i in range(npanels):
        
        pn = fig.add_subplot(npanels, 1, i + 1)    
        
        if (par[i] == 3):
            yvalues = data.get('Z'); ylim = iagas.get('zcomp_range');
            ymajorticks_step = iagas.get('zcomp_step')
            myTitle = 'Z-Comp'; myYLabel = 'nT';            
        elif (par[i] == 5):
            yvalues = data.get('H'); ylim = iagas.get('hcomp_range');
            ymajorticks_step = iagas.get('hcomp_step')
            myTitle = 'H-Comp'; myYLabel = 'nT';
        else:
            break
    
        pn.plot_date(xvalues, yvalues, 'k-')

        pn.xaxis.set_major_locator(xsmajorl)
        pn.xaxis.set_minor_locator(xsminorl)
 
        ysmajorl = MultipleLocator(ymajorticks_step) 
        
        pn.yaxis.set_major_locator(ysmajorl)
        pn.yaxis.set_minor_locator(ysminorl)
    
        pn.xaxis.set_major_formatter(DateFormatter("%H:%M"))
    
        axis([xlim[0], xlim[1], ylim[0], ylim[1]])    
    
        pn.grid(True)
    
        true_dtrange = date2num([dtime[0],dtime[len(dtime)-1]])
        if (true_dtrange[1] - true_dtrange[0] >= 1.0):
            curr_xlabel = "%s-%s (UT)" % (dtime[0].strftime('%m/%d'),\
                                        dtime[len(dtime)-1].strftime('%d/%Y'))
        else :
            curr_xlabel = "%s (UT)" % (dtime[0].strftime('%m/%d/%Y'))
    
        station = gsetup.get('station')
    
        title("%s: %s" % (station.upper(), myTitle)); xlabel(curr_xlabel); ylabel(myYLabel);
    
    if (gsetup.get('graph') > 0):
        format_type = gsetup.get('format')
        fname = ('%s%04i%02i%02i.' + format_type) % (station, dtime[0].year, dtime[0].month, dtime[0].day);
        figfname = gsetup.get('path') + fname;
        savefig(figfname, format=format_type, dpi=300)
    
    show()
    
    return(0)
