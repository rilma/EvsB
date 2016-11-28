
import numpy
import pylab
import string
import matplotlib.dates
import time_utilities
import scriptutil
import pylab
import scipy.io
import magnetometer_library as maglib

from pylab import datestr2num
from scipy import asarray


def readb_spidr_ascii(filename, comments='#'):


    def GeneratorFunction():

        f = open(filename, mode='r')

        for line in f:
            if line[0] == comments and line[0:3] == (comments + 'dd'): break

        for line in f:
            currLine = line.split()
            dtnum = datestr2num('{:s} {:s}'.format(currLine[0], currLine[1]))
            value = float(currLine[2])
            yield dtnum, value

        f.close()


    value1, value2 = [], []

    g = GeneratorFunction()
    while True:
        try:
            val1, val2 = next(g)
        except:
            break
        value1.append(val1)
        value2.append(val2)

    return asarray(value1), asarray(value2)


# Read an individual SPIDR file (1 day and 1 parameter)
def read_spidr_ascii(filename):
    
    hd = readh_spidr_ascii(filename)
    
    if False:

        tmp = numpy.loadtxt(filename, comments='#', converters={0:pylab.datestr2num}, dtype='str')
        
        time0 = numpy.double(tmp[:, 0]); time1 = pylab.datestr2num(tmp[:, 1]);
        tdelta = time1[0] - time0[0]
        
        time = time1 - tdelta
        data = numpy.double(tmp[:, 2])

    else:

        time, data = readb_spidr_ascii(filename)        
    
    ind = numpy.where(data == hd.get('Missing value'))
    if (len(ind[0]) > 0): data[ind[0]] = numpy.nan
    
#    return({'time':time, hd.get('Element'):data})
    return({'time':time, 'data':data, 'header':hd})
             
def readh_spidr_ascii(filename):
    fid = open(filename, 'r')
    i = 0
    for line in fid:        
        
        if (line[1] == '>'):
            
            for line in fid:
                
                if (line[1] == '>'):
                    break
                
                dummy = line[1:len(line)-1].split(':')                
                
                # Upgrade the output structure
                if (i == 0):
                    header = {dummy[0] : dummy[1].lstrip()}
                else:
                    header[dummy[0]] = dummy[1].lstrip()
                
                i = i + 1
                
            break        
            
    fid.close
    
    header['Missing value'] = numpy.double(header.get('Missing value'))
    
    return(header)
    
# Load all the parameter available in one day
def load_spidr_all1sd(path, dt, station, filetype):

    dtobj = matplotlib.dates.num2date(dt)

    # List the compressed image files available in the "path" folder
    wildcard = 'spidr-' + station + '-' + dtobj.strftime('%Y%m%d') + ('-%s.txt' % filetype)
#    wildcard = 'spidr-' + station + '-' + dtobj.strftime('%Y%m%d') + '-*.txt'
#    wildcard = 'spidr-' + station + '-' + dtobj.strftime('%Y%m%d') + '-he.txt'    
#    wildcard = 'spidr-' + station + '-' + dtobj.strftime('%Y%m%d') + '-hp.txt'
#    wildcard = 'spidr-' + station + '-' + dtobj.strftime('%Y%m%d') + '-hn.txt'
    flist = scriptutil.ffind(path, shellglobs=(wildcard, ))

    for i in range(len(flist)):
        
        print( ('Reading: %s') % flist[i] )
        
        tmp = read_spidr_ascii(flist[i])
        
        varname = tmp.get('header').get('Element')
        
        if (i == 0):
            output = {'time':tmp.get('time'), varname:tmp.get('data'), \
                      'station':tmp.get('header').get('Station name')}
#            print 'Stop!'
        else:
            output[varname] = tmp.get('data')
            
    
    return(output)

def plot_goes_magdata(data1, data2, data3, data4, data5, gsetup):
    
    time_lim = gsetup.get('xlim'); 
    
    gms = goes_mag_setup(data1.get('station'), time_lim)    
    gms_d5 = goes_mag_setup(data5.get('station'), time_lim)                    
        
    tlabels = time_utilities.time_label(time_lim); xlabel = tlabels.get('tlabel');
    xsmajorl = matplotlib.dates.HourLocator(None, gsetup.get('xmajorticks'))
    xsminorl = matplotlib.dates.MinuteLocator(None, gsetup.get('xminorticks'))    
    fig = pylab.figure(num=1, figsize=(12, 15))
#    pn = fig.add_subplot(111)

    partype = gsetup.get('partype'); numpanels = len(partype);

    for i in range(numpanels):
        pn = fig.add_subplot(numpanels, 1, i + 1)
    
        if (partype[i] == 0):    # GOES 10 data
#            xvalues = data1.get('time'); yvalues = data1.get('he');
            xvalues = data1.get('time'); yvalues = data1.get('hp');
            yvalues_he = data1.get('he'); yvalues_hp = data1.get('hp');
            myYLabel = 'nT'; myTitle = data1.get('station') + ' Long. W 135';
            myYLim = gms.get('mag_lim'); myYMinorticks_step = gms.get('mag_minorticks');            
        elif (partype[i] == 1):  # JRO data
            xvalues = data2.get('time'); yvalues = data2.get('efield');
            myYLabel = 'mV/m'; myTitle = 'JRO E-Field';
            myYLim = [-4, 4]; myYMinorticks_step = 0.25
        elif (partype[i] == 2):
            xvalues = data3.get('time'); yvalues = data3.get('H');
            myYLabel = 'nT'; myTitle = 'Jicamarca H-Component';
            myYLim = [25600, 26000]; myYMinorticks_step = 50;
        elif (partype[i] == 3):
            xvalues = data4.get('time'); yvalues = data4.get('H');
            myYLabel = 'nT'; myTitle = 'Piura H-Component';
            myYLim = [27350, 27650]; myYMinorticks_step = 50;
        elif (partype[i] == 4): # GOES-12
#            xvalues = data5.get('time'); yvalues = data5.get('he');
            xvalues = data5.get('time'); yvalues = data5.get('hp');
            yvalues_he = data5.get('he'); yvalues_hp = data5.get('hp');
            myYLabel = 'nT'; myTitle = data5.get('station') + ' Long. W 75';
            myYLim = gms_d5.get('mag_lim'); myYMinorticks_step = gms_d5.get('mag_minorticks');
        elif (partype[i] == 12):
            xvalues = data3.get('time'); yvalues = numpy.diff(data3.get('H'));
            yvalues = numpy.append(yvalues, numpy.nan)
            myYLabel = 'nT/min'; myTitle = 'Jicamarca dH/dt';
            myYLim = [-20, 20]; myYMinorticks_step = 10;
        elif (partype[i] == 13):
            xvalues = data4.get('time'); yvalues = numpy.diff(data4.get('H'));
            yvalues = numpy.append(yvalues, numpy.nan)
            myYLabel = 'nT/min'; myTitle = 'Piura dH/dt';
            myYLim = [-20, 20]; myYMinorticks_step = 10;
        else:
            print( 'Else' )
            
        if ((partype[i] == 0) | (partype[i] == 4)):
            pn.plot_date(xvalues, yvalues_he, color='r', label='he', linestyle='-', marker='None')
            pn.plot_date(xvalues, yvalues_hp, color='g', label='hp', linestyle='-', marker='None')
            pylab.legend(loc='best')
        else:
            pn.plot_date(xvalues, yvalues, color='k', linestyle='-', marker='None')
        
#        pn.plot_date(xvalues, yvalues, color='r', label='he', linestyle='-', marker='None')    
#        pn.plot_date(xvalues, data1.get('hn'), color='b', label='hn', linestyle='-', marker='None')
#        pn.plot_date(xvalues, data1.get('hp'), color='g', label='hp', linestyle='-', marker='None')
#        pn.plot_date(xvalues, data1.get('ht'), color='k', label='ht', linestyle='-', marker='None')
    
        pn.xaxis.set_major_locator(xsmajorl)    
        pn.xaxis.set_minor_locator(xsminorl)    
        pn.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
        ysminorl = matplotlib.ticker.MultipleLocator(myYMinorticks_step)
        pn.yaxis.set_minor_locator(ysminorl)
        pn.grid(True)
    
        curr_axis = [time_lim[0], time_lim[1], myYLim[0], myYLim[1]]
        pylab.axis(curr_axis);

        pylab.ylabel(myYLabel); pylab.title(myTitle);
        
 #       pylab.legend()
    
        if i == (numpanels - 1):
            pylab.xlabel(xlabel)
        else:
            old_xticklabels = pn.axes.get_xticklabels()            
            nxticks = len(old_xticklabels)
            ticklabels = [];
            for ii in range(nxticks):
                ticklabels.append('')                
            pn.xaxis.set_ticklabels(ticklabels)

    if (gsetup.get('graph') > 0):
        format_type = gsetup.get('format')
        fname = ('goes-mag-%s.' + format_type) % tlabels.get('tfname')
        figfname = gsetup.get('path') + fname;
        pylab.savefig(figfname, format=format_type, dpi=300)

    pylab.show()
    
    return(0)

def goes_mag_setup(station, dtrange):
    
    if (station == 'GOES-10'):
        dt_lim = time_utilities.dtlim(2004, 11, [1, 30], [0, 23], [0, 59], [0, 59])
        if ((dtrange[0] >= dt_lim[0]) & (dtrange[1] <= dt_lim[1])):
            mag_lim = [-150, 200]; mag_minorticks_step = 10;
    
    if (station == 'GOES-12'):
        dt_lim = time_utilities.dtlim(2004, 11, [1, 30], [0, 23], [0, 59], [0, 59])
        if ((dtrange[0] >= dt_lim[0]) & (dtrange[1] <= dt_lim[1])):
            mag_lim = [-100, 250]; mag_minorticks_step = 10;
    
    output = {'mag_lim':mag_lim, 'mag_minorticks':mag_minorticks_step}

    return(output)
    
