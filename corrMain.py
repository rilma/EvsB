if __name__ == '__main__':

    import file_library as file_lib
    import graphical_library as graph_lib
    import proc_pef2 as prcpef
    import pylab
    import time_utilities as time_util


    figid = [2,5,6,7]

    for id in figid:

        if id == 2:
            # Fig. 2
            trange = time_util.dtlim(2004, 11, [9, 10], [12, 20], [0, 59], [0, 59])    
            tmajtick = 6; ntmintick = 6;
            ipar = pylab.array([[1]])

        # elif id == 3:
        #     # Fig. 3
        #     trange = time_util.dtlim(2004, 11, 10, [5, 10], [0, 59], [0, 59])    
        #     tmajtick = 1; ntmintick = 3;
        #     ipar = pylab.array([[1],[2]])

        # elif id == 4:
        #     # Fig. 4
        #     trange = time_util.dtlim(2004, 11, 10, [7, 8], [0, 59], [0, 59])    
        #     tmajtick = 1./3; ntmintick = 4;
        #     ipar = pylab.array([[1],[3]])

        elif id == 5:
            # Fig. 5
            trange = time_util.dtlim(2004, 11, 10, [7, 8], [0, 29], [0, 59])    
            tmajtick = 1./3; ntmintick = 4;
            ipar = pylab.array([[20]])

        elif id == 6:
            # Fig.6
            trange = time_util.dtlim(2004, 11, 10, [7, 8], [0, 59], [0, 59])     
            tmajtick = 1./3; ntmintick = 4;
            ipar = pylab.array([[40]])

        elif id == 7:
            # Fig.7
            tmajtick = 1.; ntmintick = 6;
            trange = time_util.dtlim(2004, 11, 10, [7, 8], [38, 28], [0, 59])
            tmajtick = 1./6; ntmintick = 4;
            ipar = pylab.array([[18]])

        #
        ########################################################################
        #

        curr_trange = trange
        time_label = time_util.time_label(curr_trange, strutc='UT')
        tlabel = time_label.get('tlabel');

        prcpef.plot_pef_1d(ipar=ipar, ntmintick=ntmintick, trange=trange, 
            tmajtick=tmajtick, tlabel=tlabel)

        fnameprefix = 'fig%02i-pef-' % id
        figformat = 'png'
        figsave = 1
        figshow = 0
        figdpi = 300

        # Folder which plots are saved
        gpath = '/home/rilma/tmp/results/pef/'
        file_lib.create_directory(gpath)

        figname = gpath + fnameprefix + time_label.get('tfname') + '.' + figformat
        svfig = {'save' : figsave, 'show' : figshow, 'dpi' : figdpi,
            'figname' : figname, 'format' : figformat}
        if (svfig.get('save') > 0): graph_lib.save_figure(svfig)
        if (svfig.get('show') > 0): pylab.show()
