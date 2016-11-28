
import datetime
import matplotlib.dates
import matplotlib.ticker
import numpy
import os
import string
import time_utilities
    
def read_jromag_ascii(filename):
    
    print( 'Reading %s' % filename )
    
    tmp = numpy.loadtxt(filename, skiprows=4)
    
    dobj = datetime.date(int(tmp[0, 2]), int(tmp[0, 1]), int(tmp[0, 0]))    
    time = []
    
    for i in range(len(tmp[:, 0])):
        tobj = datetime.time(int(tmp[i, 3]), int(tmp[i, 4]), 0)
        time.append(datetime.datetime.combine(dobj, tobj))

    for i in (numpy.arange(3) + 5):
        tmp2= time_utilities.fillgaps1D(matplotlib.dates.date2num(time), tmp[:, i], 60)
        if i == 5:
            time2 = tmp2.get('xvalues'); dcomp = tmp2.get('yvalues');
        elif i == 6:
            hcomp = tmp2.get('yvalues')
        elif i == 7:
            zcomp = tmp2.get('yvalues')            

    output = {'time':time2, 'D':dcomp, 'H':hcomp, 'Z':zcomp}
    
    return(output)

def load_jromag_ascii(dtrange, path, station='jic'):
    time = []; dcomp = []; hcomp = []; zcomp = [];
    datetime_date = matplotlib.dates.num2date(dtrange)
    dtvalues = matplotlib.dates.drange(datetime_date[0], datetime_date[1],\
                            datetime.timedelta(days=1))                            
    for i in range(len(dtvalues)):
        dobj = matplotlib.dates.num2date(dtvalues[i])
        stryy = "%s" % dobj.year; yy = '%s%s' % (stryy[2],stryy[3]);
        fname = station + dobj.strftime('%d%b.%ym').lower()
        filename = path + station + yy + '/' + fname
        if (not os.access(filename, os.R_OK)):
            print( 'The file %s could not be found!' % filename )
        else:
            jromag = read_jromag_ascii(filename)
            time = numpy.append(time, jromag.get('time'))
            dcomp = numpy.append(dcomp, jromag.get('D'))
            hcomp = numpy.append(hcomp, jromag.get('H'))
            zcomp = numpy.append(zcomp, jromag.get('Z'))

    output = {'station':station, 'time':time, 'D':dcomp, 'H':hcomp, 'Z':zcomp}
    
    return(output)

def jromag_setup(dtrange):

    dt_lim = time_utilities.dtlim(2004, 11, [6, 16], [0, 23], [0, 59], [0, 59])
    if ((dtrange[0] >= dt_lim[0]) & (dtrange[1] <= dt_lim[1])):
        jic_hcomp_lim = [25525, 26150]; jic_hminorticks_step = 10
        jic_hcomp_zero = 25900.0; # Midnight value on previous quiet day: Nov 15, 2004
        piu_hcomp_lim = [27225, 27750]; piu_hminorticks_step = 10
        piu_hcomp_zero = 27517.0 + 106.5; # Midnight value on previous quiet day: Nov 15, 2004        

    dt_lim = time_utilities.dtlim(2004, 11, 10, [0, 23], [0, 59], [0, 59])
    if ((dtrange[0] >= dt_lim[0]) & (dtrange[1] <= dt_lim[1])):
        jic_hcomp_lim = [25600, 26000]; jic_hminorticks_step = 100
        jic_hcomp_zero = 25900.0; # Midnight value on 1st quiet day after storm: Nov 15, 2004
        piu_hcomp_lim = [27350, 27650]; piu_hminorticks_step = 100
        piu_hcomp_zero = 27517.0 + 106.5; # Midnight value on 1st quiet day after storm: Nov 15, 2004        

    dt_lim = time_utilities.dtlim(2005, 3, [1, 31], [0, 23], [0, 59], [0, 59])
    if ((dtrange[0] >= dt_lim[0]) & (dtrange[1] <= dt_lim[1])):
     jic_hcomp_lim = [25750, 26000]; jic_hminorticks_step = 10.
     jic_hcomp_zero = 25810.0 # Midnight value on 2nd quietest day of month: Mar 22, 2005
     piu_hcomp_lim = [27500, 27650]; piu_hminorticks_step = 5.
     piu_hcomp_zero = 27546.1 # Midnight value on 2nd quietest day of month: Mar 22, 2005

    dt_lim = time_utilities.dtlim(2008, 6, [14, 21], [0, 23], [0, 59], [0, 59])
    if ((dtrange[0] >= dt_lim[0]) & (dtrange[1] <= dt_lim[1])):
        jic_hcomp_lim = [25750, 26000]; jic_hminorticks_step = 10
        jic_hcomp_zero = 0.0; # Midnight value on previous quiet day: ??? ??, ????
        piu_hcomp_lim = [27150, 27400]; piu_hminorticks_step = 10        
        piu_hcomp_zero = 0.0; # Midnight value on previous quiet day: ??? ??, ????        
        
    jic = {'hcomp_lim':jic_hcomp_lim, 'hcomp_minorticks_step':jic_hminorticks_step, 'hcomp_zero':jic_hcomp_zero}
    piu = {'hcomp_lim':piu_hcomp_lim, 'hcomp_minorticks_step':piu_hminorticks_step, 'hcomp_zero':piu_hcomp_zero}
    
    output = {'jic':jic, 'piu':piu}
    
    return(output)
